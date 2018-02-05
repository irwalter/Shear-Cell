!Units in SI (m, N, Pa, etc)
program cube_periodic
use data_structures
use ini_posit
use torque
use segment_parameter_calc
use motion
use position_update
use output
use excl_volume
use breakage
use omp_lib
use set_up_grids

implicit none

real(8)                                :: E_Young, min_curv, r_fiber, viscosity, ex_vol_const,&
                                          gamma_dot, epsilon_dot, dt, inertia_moment, fric_coeff, distanceFactor,alpha
integer(8)                             :: nbr_neighbors, flow_case, nbr_intgr, writ_period, break_period,&
                                          i, j, k, n, frame
logical                                :: periodic_boundary, allow_breakage
real(8), dimension(3)                  :: box_size

type(rod)  , allocatable, dimension(:) :: hinges
type(fiber), allocatable, dimension(:) :: fibers
                                           
real(8), parameter                     :: pi=3.141592

integer(8), dimension(:,:), allocatable:: neighbor_list

real(8), dimension(:,:), allocatable   :: distance_neighbors
logical                                :: recover_simulation
type(segment), dimension(:),allocatable:: ghost_segments
type(cell), dimension(:), allocatable  :: cells
real(8)                                :: start, finish, start2, finish2, timex, t
logical                                :: is_fric_wall, printVelocities 
type(simulationParameters)             :: simParameters
!*******************************************************************
! Default values
simParameters%IsPeriodicY =.false. 
print *, "Maximum number of threads" , omp_get_max_threads( )  

#ifdef TENSOR            
print *, "Hydrodynamic representation being used is TENSOR"
#else
print *, "Hydrodynamic representation being used is BEAD"
#endif

!$OMP PARALLEL
print *, "Number of threads being used" , omp_get_num_threads( )  
!$OMP END PARALLEL
!print *,"1"
 call read_data(frame,&
                recover_simulation,&
                fric_coeff,&
                is_fric_wall,&
                E_Young,&
                min_curv,&
                r_fiber,&
                viscosity,&
                ex_vol_const,&
                nbr_neighbors,&
                gamma_dot,&
                epsilon_dot,&
                flow_case,&
                periodic_boundary,&
                box_size,&
                dt,&
                nbr_intgr,&
                writ_period,&
                allow_breakage,&
                break_period,&
                fibers,&
                hinges,&
                printVelocities,&
                distanceFactor,&
                simParameters)


open(3,file='OUTPUT/positions.out')
open(5,file='OUTPUT/vels.out')
open(6,file='OUTPUT/forces.out')

Inertia_Moment=(pi/4.d0)*r_fiber**4d0



if(recover_simulation.eqv..true.) then 
	n=frame*writ_period+1
else 
	n=1
	frame =1
end if
start = OMP_get_wtime()

call intializeVariables(hinges, r_fiber, gamma_dot, epsilon_dot, flow_case, simParameters )

do i=n,  nbr_intgr
!print *,"i ",i 
    t = dt* (i-1)
    
    call cpu_time (start2)

    !print *, "Main 1"
 	if (MODULO(i,break_period)==0 .or. (i.eq.n)) then
 	   print *,"Integration", i
 	   
 	   call hinges_break_config_grids_neighbors (fibers,&
 	                                             hinges,&
 	                                             ghost_segments,&
 	                                             allow_breakage,&
 	                                             min_curv,&
 	                                             r_fiber,&
 	                                             box_size,&
 	                                             cells,&                                  
 	                                             nbr_neighbors,&
                                                 neighbor_list,&
                                                 distance_neighbors,&
                                                 gamma_dot,&
                                                 t,&
                                                 distanceFactor,&
                                                 simParameters)
                                                 !print *,"out config"
        !call cpu_time(finish)
        !print *, "Neighbors", finish-start
         
 	end if
 	
 	!do j=1, ubound (neighbor_list,1)
 	!    
 	!    print *, j, "Vecino", neighbor_list(j,:)
 	!end do
 	!stop
 	!print *,"Main 2"
	!print *,"5"
	
	!call cpu_time(start)
#ifdef TENSOR    
    call fiber_par_calc_tensor(fibers,&
                        hinges,&
                        r_fiber,&
                        viscosity,&
                        gamma_dot,&
                        epsilon_dot,&
                        flow_case,&
                        simparameters)
#else    
    call fiber_par_calc(fibers,&
                        hinges,&
                        r_fiber,&
                        viscosity,&
                        gamma_dot,&
                        epsilon_dot,&
                        flow_case)
#endif
    !call cpu_time(finish)
   ! print *, "Fiber Parameters", finish-start                    
	!print *,"6"
	!do j=1,ubound(neighbor_list,1)
	!    print *,"Neighbor List", neighbor_list(j,:)
	!end do
	
    !call cpu_time(start)
 	call ex_vol_forces_moments_total(fibers,&
                                     hinges,&
                                     ghost_segments,&
                                     r_fiber,&
                                     ex_vol_const,&
                                     nbr_neighbors,&
                                     neighbor_list,&
                                     distance_neighbors,&
                                     fric_coeff,&
                                     box_size,&
                                     distanceFactor,&
                                     gamma_dot,&
                                     t)
    !call cpu_time(finish)
    !print *, "Interactions", finish-start
    !call cpu_time(start)
    if(.NOT. simParameters%IsPeriodicY ) then
    call ex_vol_forces_moments_walls2(fibers,&
                                        hinges,&
                                        r_fiber,&
                                        ex_vol_const,&
                                        box_size,&
                                        fric_coeff,&
                                        is_fric_wall,&
                                        gamma_dot)  
    end if
!    call cpu_time(finish)
    !print *, "Interactions walls", finish-start
    !call cpu_time(start)
 	call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)
 	!call cpu_time(finish)
    
    !print *, "time for bending", finish-start
 	!call cpu_time(start)
    timex=0;
    call mot(fibers, hinges, r_fiber)
    !call cpu_time(finish)
    !print *, "Dealing with matrix ", finish-start
    
    !call cpu_time(start)
 	call update_periodic(fibers, hinges, dt, periodic_boundary, box_size, gamma_dot,dt* (i-1) )
    !call cpu_time(finish)
    !print *, "Update", finish-start
 	if (MODULO(i,writ_period)==0) then
 		call data_out(fibers, hinges, frame,printVelocities)
        !call ComputeOutputFandT(fibers, hinges, r_fiber, viscosity, gamma_dot, epsilon_dot, flow_case)
 		frame=frame+1
 	end if
 	!call cpu_time (finish2)
 	!print *, "total incluye busqueda", finish2-start2	
end do 
!call ComputeOutputFandT(fibers, hinges, r_fiber, viscosity, gamma_dot, epsilon_dot, flow_case)
call cpu_time(finish)
print *, "Time Elpased",  OMP_get_wtime()-start, "s"
 close (3) 
 close (5)
  close (6)
write( *, * ) 'Press Enter to continue' 
read( *, * ) 
end program cube_periodic
