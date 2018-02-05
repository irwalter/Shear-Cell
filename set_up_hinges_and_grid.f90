module sorts
use data_structures
contains
Subroutine QCKSRT(N,ARR, Targ)

Parameter(M=7,NSTACK=50,FM=7875.,FA=211.,FC=1663.,FMI=1./FM)
! Here M is the size of subarrays sorted by straight insertion, 
! NSTACK is the required auxuliary storage, and the remaining
! constants are used by the random generating statements.
   integer, dimension(N):: ARR
   integer              :: A
   type (segment), dimension(N):: Targ
   type (segment):: B
   Dimension ISTACK(NSTACK)
   JSTACK=0
   L=1; IR=N; FX=0.
10 If (IR-L.LT.M) Then
     Do 13 J=L+1, IR
	   A=ARR(J)
	   B=Targ(J)
	   Do 11 I=J-1,1,-1
	     If (ARR(I).LE.A) Go To 12
		 ARR(I+1)=ARR(I)
		 Targ(I+1)=Targ(I)
11     Continue
       I=0
12     ARR(I+1)=A
       Targ(I+1)=B
13   Continue
     If(JSTACK.EQ.0) Return
	 IR=ISTACK(JSTACK)
	 L=ISTACK(JSTACK-1)
	 JSTACK=JSTACK-2
   else
     I=L; J=IR
	 FX=MOD(FX*FA+FC,FM) !Generate a random integer IQ between L and IR inclusive.
	 IQ=L+(IR-L+1)*(FX*FMI)
	 A=ARR(IQ)
	 B=Targ(IQ)
	 ARR(IQ)=ARR(L)
	 Targ(IQ)=Targ(L)
20   Continue
       If (J.GT.0) Then
21       If(A.LT.ARR(J)) Then
           J=J-1
		   Go To 21
         End If
	   End If
	   If (J.LE.I) Then
	     ARR(I)=A
	     Targ(I)=B
		 Go To 30
	   End If
	   ARR(I)=ARR(J)
	   Targ(I)=Targ(J)
	   I=I+1
22     If (I.LE.N) Then
         If (A.GT.ARR(I)) Then
		   I=I+1
		   Go To 22
		 End If
	   End If
	   If (J.LE.I) Then
	     ARR(J)=A
	     Targ(J)=B
		 I=J
		 Go To 30
	   End If
	   ARR(J)=ARR(I) 	  	 
	   Targ(J)=Targ(I)   
	   J=J-1	 	 
     Go To 20
30   JSTACK=JSTACK+2
     If (JSTACK.GT.NSTACK) Pause 'NSTACK too small!'
	 If (IR-I.GE.I-L) Then
	   ISTACK(JSTACK)=IR
	   ISTACK(JSTACK-1)=I+1
	   IR=I-1
	 Else
	   ISTACK(JSTACK)=I-1
	   ISTACK(JSTACK-1)=L
	   L=I+1
	 End If    	 	 	    		   
   End If
   Go To 10

end subroutine QCKSRT    
end module sorts 
    module qsort_c_module
use data_structures
implicit none
public :: QsortC
private :: Partition

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine QsortC(A, Targ)
  integer, intent(in out), dimension(:)           :: A
  type (segment), intent(in out), dimension(:) :: Targ
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, Targ, iq)
     call QsortC(A(:iq-1), Targ(:iq-1))
     call QsortC(A(iq:), Targ(iq:))
  endif
end subroutine QsortC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Partition(A, Targ, marker)
  integer, intent(in out), dimension(:) :: A
  type (segment), intent(in out), dimension(:):: Targ
  integer, intent(out) :: marker
  integer :: i, j
  real :: temp
  type(segment):: tempTarg, xTarg
  real :: x      ! pivot point
  
  x = A(1)
  !xTarg=Targ(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        tempTarg=Targ(i)
        A(i) = A(j)
        Targ(i)=Targ(j)
        A(j) = temp
        Targ(j)=tempTarg
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module             
    
    module set_up_grids
use data_structures
use torque
use breakage
use excl_volume
use sorts
use qsort_c_module 
implicit none
contains
!----------------------------------------------------------------------------- 
subroutine hinges_break_config_grids_neighbors (fibers,&
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
                                                time,&
                                                distance_factor, &
                                                simParameters)
implicit none
type(fiber), allocatable, dimension(:)  :: fibers
real(8)                                 :: E_Young, Inertia_Moment, min_curv, r_fiber, max_length, bin_length, finish2, start2
type (rod), dimension(:), allocatable   :: hinges
integer                                 :: i, j, k, m, nbr_segments, ind, numClones
integer(8)                              ::nbr_neighbors
integer, dimension(3)                   :: Nbr_bins
type (segment), dimension(:),allocatable:: ghost_segments
real(8), dimension(3)                   :: box_size,min_coor, max_coor, center
integer(8), dimension(:,:), allocatable :: neighbor_list
real(8), dimension(:,:), allocatable    :: distance_neighbors
type(cell), dimension(:), allocatable   :: cells
integer, dimension(:,:), allocatable    :: hinge_index
integer, dimension(:), allocatable      :: ghost_segment_index
integer, dimension(3)                   :: indx
real(8)                                 :: finish3, start3, gamma_dot, time, dX, distance_factor
logical                                 :: allow_breakage
type(simulationParameters)             :: simParameters
!-----------------------------------------------------------------------------
!call cpu_time(start3)

if(simParameters%IsPeriodicY ) then
    numClones =9
else
    numClones =5
end if


!call cpu_time(start2)
call bending_torque_whole(fibers, hinges, E_Young, Inertia_Moment)
if (allow_breakage)then
    call fiber_damage(fibers, hinges, min_curv, r_fiber)
end if
!print *, "set 1"
!call cpu_time(finish2)
!print *,"bemding torque and fiber damage", finish2-start2
!call cpu_time (start2)
nbr_segments=0 

! Count the number of segments
do i=1, ubound(fibers,1)       
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        nbr_segments=nbr_segments+1
    end do
end do
!call cpu_time(finish2)
!print *,"Counting segments", finish2-start2

!print *, "set 2"
!call cpu_time(start2)
if (allocated(ghost_segments)) deallocate (ghost_segments)	

allocate(ghost_segments(nbr_segments*numClones)) ! change 5 to 9 if you want periodicity in all directions

if (allocated(hinge_index)) deallocate (hinge_index)	
allocate(hinge_index(nbr_segments,2))

if (allocated(ghost_segment_index)) deallocate (ghost_segment_index)
allocate(ghost_segment_index(numClones*nbr_segments)) ! change 5 to 9 if you want periodicity in all directions
!call cpu_time(finish2)
!print *,"Allocating", finish2-start2

k=1
!print *,"set 3"
!call cpu_time(start2)

dX = mod(gamma_dot * box_size(2)/2*time,box_size(1))  ! This is how much the image boxes need to be translated for the lees edward boundaries

do i=1, ubound(fibers,1)   
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        !This is the original segment
        ghost_segments(k)%A=hinges(j)%X_i
        ghost_segments(k)%B=hinges(j+1)%X_i
        ghost_segments(k)%orig_pos(1)=j
        ghost_segments(k)%orig_pos(2)=j+1
        ghost_segments(k)%axis_loc=1
        
        ! Image segment in +x
        ghost_segments(k+1)%A(2:3)=hinges(j)%X_i(2:3)
        ghost_segments(k+1)%B(2:3)=hinges(j+1)%X_i(2:3)
        ghost_segments(k+1)%A(1)=hinges(j)%X_i(1)+box_size(1)
        ghost_segments(k+1)%B(1)=hinges(j+1)%X_i(1)+box_size(1)
        ghost_segments(k+1)%orig_pos(1)=j
        ghost_segments(k+1)%orig_pos(2)=j+1
        ghost_segments(k+1)%axis_loc=2
        
        ! Image segment in -x
        ghost_segments(k+2)%A(2:3)=hinges(j)%X_i(2:3)
        ghost_segments(k+2)%B(2:3)=hinges(j+1)%X_i(2:3)
        ghost_segments(k+2)%A(1)=hinges(j)%X_i(1)-box_size(1)
        ghost_segments(k+2)%B(1)=hinges(j+1)%X_i(1)-box_size(1)
        ghost_segments(k+2)%orig_pos(1)=j
        ghost_segments(k+2)%orig_pos(2)=j+1
        ghost_segments(k+2)%axis_loc=3
        
        ! Image segment in +z
        ghost_segments(k+3)%A(1:2)=hinges(j)%X_i(1:2)
        ghost_segments(k+3)%B(1:2)=hinges(j+1)%X_i(1:2)
        ghost_segments(k+3)%A(3)=hinges(j)%X_i(3)+box_size(3)
        ghost_segments(k+3)%B(3)=hinges(j+1)%X_i(3)+box_size(3)
        ghost_segments(k+3)%orig_pos(1)=j
        ghost_segments(k+3)%orig_pos(2)=j+1
        ghost_segments(k+3)%axis_loc=4
        
        ! Image segment in +z
        ghost_segments(k+4)%A(1:2)=hinges(j)%X_i(1:2)
        ghost_segments(k+4)%B(1:2)=hinges(j+1)%X_i(1:2)
        ghost_segments(k+4)%A(3)=hinges(j)%X_i(3)-box_size(3)
        ghost_segments(k+4)%B(3)=hinges(j+1)%X_i(3)-box_size(3)
        ghost_segments(k+4)%orig_pos(1)=j
        ghost_segments(k+4)%orig_pos(2)=j+1
        ghost_segments(k+4)%axis_loc=5
        
        if (simParameters%IsPeriodicY ) then
        ! Image segment in +Ya
        ghost_segments(k+5)%A(1)=hinges(j)%X_i(1) + dX
        ghost_segments(k+5)%B(1)=hinges(j+1)%X_i(1) + dX
        ghost_segments(k+5)%A(2)=hinges(j)%X_i(2) +box_size(2)
        ghost_segments(k+5)%B(2)=hinges(j+1)%X_i(2) +box_size(2)
        ghost_segments(k+5)%A(3)=hinges(j)%X_i(3)
        ghost_segments(k+5)%B(3)=hinges(j+1)%X_i(3)
        ghost_segments(k+5)%orig_pos(1)=j
        ghost_segments(k+5)%orig_pos(2)=j+1
        ghost_segments(k+5)%axis_loc=6
        
        ! Image segment in +Yb
        ghost_segments(k+6)%A(1)=hinges(j)%X_i(1) + dX - box_size(1)
        ghost_segments(k+6)%B(1)=hinges(j+1)%X_i(1) + dX - box_size(1)
        ghost_segments(k+6)%A(2)=hinges(j)%X_i(2) +box_size(2)
        ghost_segments(k+6)%B(2)=hinges(j+1)%X_i(2) +box_size(2)
        ghost_segments(k+6)%A(3)=hinges(j)%X_i(3)
        ghost_segments(k+6)%B(3)=hinges(j+1)%X_i(3)
        ghost_segments(k+6)%orig_pos(1)=j
        ghost_segments(k+6)%orig_pos(2)=j+1
        ghost_segments(k+6)%axis_loc=7
        
        ! Image segment in -Ya
        ghost_segments(k+7)%A(1)=hinges(j)%X_i(1) - dX
        ghost_segments(k+7)%B(1)=hinges(j+1)%X_i(1) - dX
        ghost_segments(k+7)%A(2)=hinges(j)%X_i(2) -box_size(2)
        ghost_segments(k+7)%B(2)=hinges(j+1)%X_i(2) -box_size(2)
        ghost_segments(k+7)%A(3)=hinges(j)%X_i(3)
        ghost_segments(k+7)%B(3)=hinges(j+1)%X_i(3)
        ghost_segments(k+7)%orig_pos(1)=j
        ghost_segments(k+7)%orig_pos(2)=j+1
        ghost_segments(k+7)%axis_loc=8
        
        ! Image segment in -Yb
        ghost_segments(k+8)%A(1)=hinges(j)%X_i(1) - dX + box_size(1)
        ghost_segments(k+8)%B(1)=hinges(j+1)%X_i(1) - dX + box_size(1)
        ghost_segments(k+8)%A(2)=hinges(j)%X_i(2) -box_size(2)
        ghost_segments(k+8)%B(2)=hinges(j+1)%X_i(2) -box_size(2)
        ghost_segments(k+8)%A(3)=hinges(j)%X_i(3)
        ghost_segments(k+8)%B(3)=hinges(j+1)%X_i(3)
        ghost_segments(k+8)%orig_pos(1)=j
        ghost_segments(k+8)%orig_pos(2)=j+1
        ghost_segments(k+8)%axis_loc=9
        
        end if
      
        
        
        k=k+numClones ! change 5 to 9
    end do
end do

!print *,"set 4"
min_coor=huge(0d0)
max_coor=-huge(0d0)
max_length=0
	    
do i=1, ubound(fibers,1)
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2    
        max_length=max(max_length,sqrt(dot_product((hinges(j+1)%X_i-hinges(j)%X_i),(hinges(j+1)%X_i-hinges(j)%X_i))))  
    end do
end do
!call cpu_time(finish2)
!print *,"Long chunck", finish2-start2
!print *,"set 5"
!print *,"Max_lenght",max_length
!call cpu_time(start2)
do k=1,3
    do i=1, ubound(ghost_segments,1)
        min_coor(k)=min(min_coor(k), ghost_segments(i)%A(k))
        min_coor(k)=min(min_coor(k), ghost_segments(i)%B(k))
	    max_coor(k)=max(max_coor(k), ghost_segments(i)%A(k))
	    max_coor(k)=max(max_coor(k), ghost_segments(i)%B(k))
	end do  
end do
!call cpu_time(finish2)
!print *, "Find Minumum and Maximum", finish2-start2
!print *,"Max Coor", max_coor(1), max_coor(2), max_coor(3)
!print *,"Min Coor", min_coor(1), max_coor(2), max_coor(3)
!print *,"set 6"	    


!bin_length=1.25*(max_length+2*r_fiber)!Mod 9/28/2014
bin_length=1.5*(max_length+2*r_fiber)
!min_coor= min_coor- 2*r_fiber
!max_coor= max_coor+ 2*r_fiber

min_coor= min_coor- 2*r_fiber-max_length
max_coor= max_coor+ 2*r_fiber+max_length

Nbr_bins=ceiling((max_coor-min_coor)/bin_length)

!bin_length=(max_coor(3)-min_coor(3))/Nbr_bins(3)
!min_coor=min_coor-bin_length
!max_coor=min_coor+bin_length!Mod 9/28/2014

!Nbr_bins=Nbr_bins+2

min_coor=min_coor-2*bin_length
max_coor=max_coor+2*bin_length

Nbr_bins=Nbr_bins+4


!print *,"Nbr bins", Nbr_bins
!print *,"set 7"	   
do j=1,ubound(hinges,1)
    do k=1,3
        hinges(j)%indx(k)=0
    end do
end do    	
!print *, "set 8"
!m=1


!call cpu_time(start2)
do i=1, ubound(fibers,1)
    do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
        do k=1,3
            center=(hinges(j)%X_i(k)+hinges(j+1)%X_i(k))/2
            indx(k)=floor((0.5*(hinges(j)%X_i(k)+hinges(j+1)%X_i(k))-min_coor(k))/ bin_length)+1
            if(indx(k).le.2) indx(k)=2
            if(indx(k).ge.(Nbr_bins(k)-1)) indx(k)=Nbr_bins(k)-1
	    end do
	    !print *, "indx", indx
    	!print *,"ind", indx
	    hinges(j)%ind=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
	    !print *,"hinges(j)%indx", hinges(j)%ind
	   ! m=m+1
    end do
end do
!call cpu_time(finish2)
!print *, "Clasify hinges", finish2-start2

!print *, "set 9"
m=1

do j=1, ubound(ghost_segments, 1)
    do k=1,3
        center=(ghost_segments(j)%A(k)+ghost_segments%B(k))/2
        indx(k)=floor((0.5*(ghost_segments(j)%A(k)+ ghost_segments(j)%B(k))-min_coor(k))/ bin_length)+1
        if(indx(k).le.2) indx(k)=2
        if(indx(k).ge.(Nbr_bins(k)-1)) indx(k)=Nbr_bins(k)-1
    end do
    
	ghost_segment_index(m)=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
	!print *, "ind", indx
	!print *, "ghost_segment index", ghost_segment_index(m)
	m=m+1
end do

!print *, "set 10"

!call cpu_time(start2)
call QsortC(ghost_segment_index , ghost_segments)
!call cpu_time(finish2)
!print *, "quick sort", finish2-start2

!do j=1,ubound(ghost_segment_index,1 )
!    print *,"Sorted Ghost Segment Index", ghost_segment_index(j)
!end do

!stop

!print *, "set 11"
if (allocated(cells)) deallocate(cells)
!print *, "Nbr_bins", Nbr_bins
allocate(cells(Nbr_bins(1)*Nbr_bins(2)*Nbr_bins(3)))

do i=1, ubound(cells,1)
    cells(i)%ghost_limits(1)=0
    cells(i)%ghost_limits(2)=0
end do
!print *, "set 11a"
do i=1, Nbr_bins(1)
    do j=1, Nbr_bins(2)
        do k=1, Nbr_bins(3)
            ind=(k-1)*Nbr_bins(1)*Nbr_bins(2)+(j-1)*Nbr_bins(1)+i
            cells(ind)%indx(1)=i
            cells(ind)%indx(2)=j
            cells(ind)%indx(3)=k
        end do
    end do
end do
!print *, "set 12a"
!print *, "ghost_segment_index(1)", ghost_segment_index(1)
cells(ghost_segment_index(1))%ghost_limits(1)=1
!print *, "ghost_segment_index(1)", ghost_segment_index(1)
!print *, "set 12b"
!call cpu_time(start2)
do j=2, ubound(ghost_segment_index,1)
    if (ghost_segment_index(j) .ne. ghost_segment_index(j-1)) then
        cells(ghost_segment_index(j))%ghost_limits(1)=j
        cells(ghost_segment_index(j-1))%ghost_limits(2)=j-1
    end if
end do
!call cpu_time(finish2)
!print *,"Find Limits", finish2-start2
!do j=1, ubound(cells,1)
!
!    print *,"Cells Ghost Limits",cells(j)%ghost_limits(1), cells(j)%ghost_limits(2)
!
!end do
!stop
!print *, "set 13"
!-----------------------------------------------------------------------------
!call cpu_time(start2)
call find_neighbors_new(fibers,&
                    hinges,&
                    ghost_segments,&
                    nbr_neighbors,&
                    neighbor_list,&
                    distance_neighbors,&
                    r_fiber,&
                    cells,&
                    Nbr_bins,&
                    distance_factor)
!call cpu_time(finish2)
!print *,"Internal Find Neighbors", finish2-start2
!call cpu_time(finish3)
!print *, "Time hinges_break_config_measured inside", finish3-start3
end subroutine hinges_break_config_grids_neighbors      
end module set_up_grids

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  	             	