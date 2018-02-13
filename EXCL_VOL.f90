module excl_volume
use data_structures
use torque
use omp_lib
implicit none
contains
!************************************************************************************************
!************************************************************************************************
subroutine find_neighbors_new(fibers,&
                              hinges,&
                              ghost_segments,&
                              nbr_neighbors,&
                              neighbor_list,&
                              distance_neighbors,&
                              r_fiber,&
                              cells,&
                              Nbr_bins,&
                              distance_factor)

type(fiber), dimension(:)              :: fibers
type(rod), dimension(:)                :: hinges
type(rod)                              :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:)            :: ghost_segments
integer                                :: i, j, k, m, o, ii , jj, kk, indd
integer(8)                             ::nbr_neighbors
integer(8), dimension(:,:), allocatable:: neighbor_list
real(8), dimension(:,:), allocatable   :: distance_neighbors
real(8)                                :: r_fiber, threshold, Gab_min, ex_vol_const,fric_coeff
logical                                :: flag, is_for_neigh_list
real(8)                                :: epsilon, distance_factor
type(cell), dimension(:), allocatable  :: cells
integer, dimension(3)                  :: i_cell, nbr_bins
integer, dimension(3)                  :: indx
real(8), dimension(3)                  :: rad, d


!print *, "Find 1"
epsilon=3*tiny(1d0)
!print *, "Find 2"
threshold= distance_factor*r_fiber   
ex_vol_const=1 !Not really necessary, the ex_vol_forces_moments_segs
                !is used here to find the distances without being concerned with the forces...

if (allocated(neighbor_list)) then
	deallocate(neighbor_list)
end if
!print *, "Find 3"
if (allocated(distance_neighbors)) then
	deallocate (distance_neighbors)
end if

allocate(neighbor_list(ubound(hinges,1),nbr_neighbors))
allocate(distance_neighbors(ubound(hinges,1),nbr_neighbors))
!print *, "Find 4"
neighbor_list  	  =0
distance_neighbors=1e10
is_for_neigh_list=.true.

call omp_set_num_threads(2)
$OMP PARALLEL 
$OMP DO PRIVATE(i, j, i_cell,ii,indx, jj, kk, indd, k, Gab_min, flag, m, ghost_hingeA, ghost_hingeB)
do i=1, ubound (fibers,1)
    !print *, "Find 6"
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
	   ! print *,"hinges(j)%ind", hinges(j)%ind
	    i_cell=cells(hinges(j)%ind)%indx
	   ! print *,"cells(hinges(j)%ind)%indx",  cells(hinges(j)%ind)%indx
	    !print *, "Find 7"
	    do ii=-1,1
	        indx(1)=i_cell(1)+ii
	        !print *, "Find 9"
	        do jj=-1,1
	            indx(2)=i_cell(2)+jj
	            !print *, "Find 10"
	            do kk=-1,1
	               indx(3)=i_cell(3)+kk
	              ! print *,"indx1 indx2 indx3", indx(1), indx(2), indx(3)
	               indd=(indx(3)-1)*Nbr_bins(1)*Nbr_bins(2)+(indx(2)-1)*Nbr_bins(1)+indx(1)
	               !print *,"indd", indd
	               !print *, "Find 11"
	               !print *,"indd", indd
	               do k=cells(indd)%ghost_limits(1),cells(indd)%ghost_limits(2)
	               !print *, "Find 12"
	                    if (cells(indd)%ghost_limits(1).ne.0) then 
	                        !print *, "Find 13"
	                        
	                        !d=0.5*((ghost_segments(k)%A+ghost_segments(k)%B)-(hinges(j)%X_i+hinges(j+1)%X_i))
	                        !rad=0.5*(ghost_segments(k)%A-ghost_segments(k)%B-hinges(j)%X_i+hinges(j+1)%X_i)
	                        !if (dot_product(d,d).LE. 1.5*(dot_product(rad,rad)+2*r_fiber)) then
	                        if (are_boxes_intersect(hinges(j)%X_i, hinges(j+1)%X_i,ghost_segments(k)%A, ghost_segments(k)%B, r_fiber, 6*r_fiber )) then
		                        if( (abs(dot_product(ghost_segments(k)%A-hinges(j)%X_i,ghost_segments(k)%A-hinges(j)%X_i)) .ge. epsilon) .and.&  
		                            (abs(dot_product(ghost_segments(k)%A-hinges(j+1)%X_i,ghost_segments(k)%A-hinges(j+1)%X_i)) .ge. epsilon) .and.&
		                            (abs(dot_product(ghost_segments(k)%B-hinges(j)%X_i,ghost_segments(k)%B-hinges(j)%X_i)) .ge. epsilon) .and.&
		                            (abs(dot_product(ghost_segments(k)%B-hinges(j+1)%X_i,ghost_segments(k)%B-hinges(j+1)%X_i)) .ge. epsilon)) then
			                        ghost_hingeA%X_i=ghost_segments(k)%A
			                        ghost_hingeB%X_i=ghost_segments(k)%B
					                call ex_vol_forces_moments_segs2(ghost_hingeA,&
	                                        	                    ghost_hingeB,&
	                                        	                    hinges(j),&
	                                        	                    hinges(j+1),&
		                                                            ex_vol_const,&
		                                                            threshold,&
		                                                            r_fiber,&
	                                                                Gab_min,&
				                                                    is_for_neigh_list,&
				                                                    fric_coeff)
					                flag=.false.
					                m=1
				                    do while ((m.le.nbr_neighbors) .and. (flag.eqv. .false.))
							        !print *,m
					                    if ((Gab_min .lt. Threshold) .and.&
						                        (Gab_min .lt. distance_neighbors(j,m)))then
				                            if (m.ne. nbr_neighbors) then
				  		                        distance_neighbors(j,m+1:nbr_neighbors)=&
							                    distance_neighbors(j,m:nbr_neighbors-1)
						                        neighbor_list(j,m+1:nbr_neighbors)=&
					                            neighbor_list(j,m:nbr_neighbors-1)
						                    end if
						                    flag=.true.
						                    distance_neighbors(j,m)=Gab_min
				                            neighbor_list(j,m)=k
					                    end if
				                        m=m+1
		                            end do
	                            end if
	                        end if
	                    end if
                    end do	  
	            end do
	        end do
	    end do
    end do
end do
$OMP END DO
$OMP END PARALLEL

!do i=1, ubound(hinges,1)
!	print *,"hinge", i, neighbor_list(i,1:6)
!end do


end subroutine find_neighbors_new

!************************************************************************************************
subroutine find_neighbors(fibers,&
                          hinges,&
                          hinges_with_ghosts,&
                          nbr_neighbors,&
                          neighbor_list,&
                          distance_neighbors,&
                          r_fiber)

type(fiber), dimension(:)              :: fibers
type(rod), dimension(:)                :: hinges, hinges_with_ghosts
integer(8)                             :: i, j, k, l, m, o , nbr_neighbors
integer(8), dimension(:,:), allocatable:: neighbor_list
real(8), dimension(:,:), allocatable   :: distance_neighbors
real(8)                                :: r_fiber, threshold, Gab_min, ex_vol_const,fric_coeff
logical                                :: flag, is_for_neigh_list

threshold= 500*r_fiber   
ex_vol_const=1 !Not really necessary, the ex_vol_forces_moments_segs
                !is used here to find the distances without being concerned with the forces...

if (allocated(neighbor_list)) then
	deallocate(neighbor_list)
end if

if (allocated(distance_neighbors)) then
	deallocate (distance_neighbors)
end if

allocate(neighbor_list(ubound(hinges,1),nbr_neighbors))
allocate(distance_neighbors(ubound(hinges,1),nbr_neighbors))

neighbor_list  	  =0
distance_neighbors=1e10
is_for_neigh_list=.true.



do i=1, ubound (fibers,1)
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
		do k=1, ubound (fibers,1)
			do o=0,2
				do l=fibers(k)%first_hinge+o*ubound(hinges,1), fibers(k)%first_hinge+fibers(k)%nbr_hinges+o*ubound(hinges,1)-2
					if(i.ne.k) then
					!if ((i.ne.k).or. ((i.eq.k) .and. (.not.((j.eq. l).or. (j.ne.l+1) or(j.ne.l-1))
	                
						call ex_vol_forces_moments_segs2(hinges_with_ghosts(l),&
	                                        	        hinges_with_ghosts(l+1),&
	                                        	        hinges(j),&
	                                        	        hinges(j+1),&
		                                                ex_vol_const,&
		                                                threshold,&
		                                                r_fiber,&
	                                                    Gab_min,&
				                                        is_for_neigh_list,&
				                                        fric_coeff)
						flag=.false.
						m=1
						do while ((m.le.nbr_neighbors) .and. (flag.eqv. .false.))
							!print *,m
							if ((Gab_min .lt. Threshold) .and.&
						            (Gab_min .lt. distance_neighbors(j,m)))then
								if (m.ne. nbr_neighbors) then
				  					distance_neighbors(j,m+1:nbr_neighbors)=&
									distance_neighbors(j,m:nbr_neighbors-1)
									neighbor_list(j,m+1:nbr_neighbors)=&
									neighbor_list(j,m:nbr_neighbors-1)
								end if
								flag=.true.
								distance_neighbors(j,m)=Gab_min
								neighbor_list(j,m)=l
							end if
							m=m+1
						end do
					end if
				end do
			end do
		end do	
	end do
end do


!do i=1, ubound(hinges,1)
!	print *,"hinge", i, neighbor_list(i,1:6)
!end do


end subroutine find_neighbors
!************************************************************************************************
!************************************************************************************************
subroutine ex_vol_forces_moments_total(fibers,&
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
                                        time)


real(8)                                :: L_seg, r_fiber, gamma_dot, viscosity, ex_vol_const, threshold, distanceFactor
type(fiber), dimension(:)              :: fibers
type(rod), dimension(:)                 :: hinges
integer(8)                              :: i, j,k, l, m, nbr_neighbors
integer(8), dimension(:,:), allocatable :: neighbor_list
real(8), dimension(:,:), allocatable    :: distance_neighbors
real(8), parameter                      :: pi=3.141592
real(8)                                 :: Gab_min, fric_coeff, time, dX, timeA
logical                                 :: is_for_neigh_list, flag
type(rod)                               :: ghost_hingeA, ghost_hingeB
type(segment), dimension(:), allocatable:: ghost_segments
real(8),dimension(3)                    :: box_size
threshold         = 2*r_fiber
is_for_neigh_list =.false.
dX = mod(gamma_dot * box_size(2)/2*time,box_size(1))

!$OMP PARALLEL DEFAULT(SHARED) 

!$OMP DO PRIVATE (i, j, k)
do i=1, ubound(ghost_segments,1)
!print *,"i", i   ," of ",   ubound(ghost_segments,1)
    j= ghost_segments(i)%orig_pos(1)
    k= ghost_segments(i)%orig_pos(2)
    

    select case (ghost_segments(i)%axis_loc)
        case(1) ! image +X
            ghost_segments(i)%A=hinges(j)%X_i
            ghost_segments(i)%B=hinges(k)%X_i
        case(2) ! image +X
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)+box_size(1)
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)+box_size(1)
        case(3) ! image -X
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)-box_size(1)
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)-box_size(1)
        case(4) ! image +Z
            ghost_segments(i)%A(3)=hinges(j)%X_i(3)+box_size(3)
            ghost_segments(i)%B(3)=hinges(k)%X_i(3)+box_size(3)
        case(5) ! image -Z
            ghost_segments(i)%A(3)=hinges(j)%X_i(3)-box_size(3)
            ghost_segments(i)%B(3)=hinges(k)%X_i(3)-box_size(3)
        case(6) ! image +Ya
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)+dX
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)+dX
            ghost_segments(i)%A(2)=hinges(j)%X_i(2)+box_size(2)
            ghost_segments(i)%B(2)=hinges(k)%X_i(2)+box_size(2)
        case(7) ! image +Yb
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)+dX - box_size(1)
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)+dX - box_size(1)
            ghost_segments(i)%A(2)=hinges(j)%X_i(2)+box_size(2)
            ghost_segments(i)%B(2)=hinges(k)%X_i(2)+box_size(2)
        case(8) ! image -Ya
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)-dX
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)-dX
            ghost_segments(i)%A(2)=hinges(j)%X_i(2)-box_size(2)
            ghost_segments(i)%B(2)=hinges(k)%X_i(2)-box_size(2)
        case(9) ! image -Yb
            ghost_segments(i)%A(1)=hinges(j)%X_i(1)- dX + box_size(1)
            ghost_segments(i)%B(1)=hinges(k)%X_i(1)- dX + box_size(1)
            ghost_segments(i)%A(2)=hinges(j)%X_i(2)-box_size(2)
            ghost_segments(i)%B(2)=hinges(k)%X_i(2)-box_size(2)    
            
        end select
end do
!$OMP END DO !NOWAIT



!$OMP DO PRIVATE (i)
do i=1, ubound(hinges,1)
	hinges(i)%T_Excl_Vol=0d0
	hinges(i)%F_Excl_Vol=0d0
end do
!$OMP END DO



!print *,"vecinos", nbr_neighbors
    !timeA = OMP_get_wtime()

!$OMP DO PRIVATE (i, j, k, m,l, flag, Gab_min, ghost_hingeA,ghost_hingeB) !SCHEDULE(DYNAMIC)
do i=1, ubound (fibers,1)
	do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
		m=0
		flag=.true.
		do while (flag.eqv..true.)
			m=m+1
			!print *,"test"
			l=neighbor_list(j,m)
			if (l.ne. 0) then
			    ghost_hingeA%X_i=ghost_segments(l)%A
			    ghost_hingeB%X_i=ghost_segments(l)%B
				call ex_vol_forces_moments_segs2(ghost_hingeA,&
				                ghost_hingeB,&
								hinges(j),&
								hinges(j+1),&
							    ex_vol_const,&
						        threshold,&
	                            r_fiber,&
                                Gab_min,&
			                    is_for_neigh_list,&
			                    fric_coeff)
			end if
			
			!print *,"check", i, j, m
			
			if ((m.eq.nbr_neighbors).or.(neighbor_list(j,m).eq. 0)) then
				flag=.false.
			end if
			
		end do			
	end do
end do
!$OMP END DO !NOWAIT

!!$omp single
!do i=1, ubound(hinges,1)
!	print *, hinges(i)%T_Excl_Vol
!	print *, hinges(i)%F_Excl_Vol
!end do
!!$omp end single

!$OMP END PARALLEL
    !print *, " time elpased excluded volume " , OMP_get_wtime()- timeA
end subroutine ex_vol_forces_moments_total
!************************************************************************
subroutine ex_vol_forces_moments_segs(hingesa1,&
                                      hingesa2,&
                                      hingesb1,&
                                      hingesb2,&
                                      ex_vol_const,&
                                      threshold,&
                                      R,&
                                      Gab_min,&
				                      is_for_neigh_list,&
				                      fric_coeff)

type(rod)                           :: hingesb1,&
                                       hingesb2,&
                                       hingesa1,&
                                       hingesa2
type (vec_arr)       , dimension(7) :: vec
type (simple_vec_arr), dimension(7) ::vec_a
real(8)                             :: threshold, R, Gab_min, ex_vol_const, fric_coeff, Sab
integer(8)                          :: i,j, k,l, ind1, ind2
real(8), dimension(3)               :: void, Exc_Vol_Force_Partial, v_rel
logical                             :: are_near
logical                             :: is_for_neigh_list

 call dist_segs(hingesa1%X_i,&
               hingesa2%X_i,&
               hingesb1%X_i,&
               hingesb2%X_i,&
               void,&
               vec(1)%pb,&
               vec(1)%Gab,&
               vec(1)%Gab_norm,&
	           vec(1)%Sba,&
               vec(1)%coll_course,&
               Sab)
 vec_a(1)%r=(hingesa2%X_i-hingesa1%X_i)/sqrt(dot_product(hingesa2%X_i-hingesa1%X_i, hingesa2%X_i-hingesa1%X_i))
 vec_a(1)%r=Sab*vec_a(1)%r
 
 call dist_pt_seg(hingesa1%X_i,&
                 hingesb1%X_i,&
                 hingesb2%X_i,&
                 vec(2)%pb,&
                 vec(2)%Gab,&
                 vec(2)%Gab_norm,&
                 vec(2)%Sba,&
                 vec(2)%coll_course)
 vec_a(2)%r=0

 call dist_pt_seg(hingesa2%X_i,&
                 hingesb1%X_i,&
                 hingesb2%X_i,&
                 vec(3)%pb,&
                 vec(3)%Gab,&
                 vec(3)%Gab_norm,&
		         vec(3)%Sba,&
                 vec(3)%coll_course)
vec_a(3)%r=(hingesa2%X_i-hingesa1%X_i)

do i=1,3
	vec(i)%r=vec(i)%pb*vec(i)%Sba
end do

 call dist_pt_pt(hingesa1%X_i,&
                hingesb1%X_i,&
                vec(4)%Gab,&
                vec(4)%Gab_norm)
 
vec_a(4)%r=0
 call dist_pt_pt(hingesa2%X_i,&
                hingesb1%X_i,&
                vec(5)%Gab,&
                vec(5)%Gab_norm)
  
vec_a(5)%r=(hingesa2%X_i-hingesa1%X_i)

vec(4)%r=0
vec(5)%r=0


 call dist_pt_pt(hingesa1%X_i,&
                hingesb2%X_i,&
                vec(6)%Gab,&
                vec(6)%Gab_norm)
                
vec(6)%r= hingesb2%X_i-hingesb1%X_i                
vec_a(6)%r=0

 call dist_pt_pt(hingesa2%X_i,&
                hingesb2%X_i,&
                vec(7)%Gab,&
                vec(7)%Gab_norm)
vec_a(7)%r=(hingesa2%X_i-hingesa1%X_i)
vec(7)%r= hingesb2%X_i-hingesb1%X_i


vec(4:7)%coll_course=.true. !Unlike lines, Points are always in collition course...

k=7

Gab_min=vec(7)%Gab_norm

do i=6,1,-1
	if ((vec(i)%coll_course.eqv. .true.) .and. (vec(i)%Gab_norm .lt. Gab_min)) then
		Gab_min=vec(i)%Gab_norm
		k=i
	end if 
end do

Exc_Vol_Force_Partial=0
are_near=.false.




if ((is_for_neigh_list.eqv. .false.).and.&
    (Gab_min .lt. threshold) .and.&
    (vec(k)%coll_course.eqv..true.)) then
	are_near=.true.
	!Exc_Vol_Force_Partial =-vec(k)%Gab*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))!CORRECT ONE
    
    !print *, vec(k)%Gab
	Exc_Vol_Force_Partial =-(vec(k)%Gab/(vec(k)%Gab_norm))*2*r*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))
	if (hingesb1%is_stationary==1) hingesb1%v_i=0
	if (hingesa1%is_stationary==1) hingesa1%v_i=0
	
    v_rel= -hingesb1%v_i+hingesa1%v_i -cross(hingesb1%omega,vec(k)%r)&
                                      +cross(hingesa1%omega,vec_a(k)%r)
    
    v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
    
    !print *, "v_rel", v_rel
    !print *,"fric_coeff", fric_coeff         
    if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)) ) then !Check for NaN. If NaN, vrel
        Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	end if
	
	hingesb1%F_Excl_Vol   = hingesb1%F_Excl_Vol+ Exc_Vol_Force_Partial
	hingesb1%T_Excl_Vol   = hingesb1%T_Excl_Vol +cross((vec(k)%r -hingesb1%r/2 ), Exc_Vol_Force_Partial)
    !print *, " Exc_Vol_Force_Partial " , Exc_Vol_Force_Partial
    !print *, " T " , hingesb1%T_Excl_Vol
    !print *, " r " , vec(k)%r
	!print *, "Collision case", k, "Indexes", ind1, ind2
end if



end subroutine ex_vol_forces_moments_segs
!**************************************************************************************		
subroutine dist_segs(ra,&
                     ra_end,&
                     rb,&
                     rb_end,&
                     pa,&
                     pb,&
                     Gab,&
                     Gab_norm,&
                     Sba,&
                     coll_course,&
                     Sab)

real(8), dimension(3):: ra, ra_end, rb, rb_end, pa, pb, Gab
real(8)              :: Sab, Sba, la, lb, Gab_norm
logical              :: coll_course

pa=ra_end-ra
pb=rb_end-rb

la=sqrt(dot_product(pa,pa))
lb=sqrt(dot_product(pb,pb))

pa=pa/la
pb=pb/lb

Sab=(dot_product(ra-rb, pb)*dot_product(pa,pb)-dot_product(ra-rb,pa))&
    /(1d0-dot_product(pa, pb)**2d0)

Sba=(dot_product(rb-ra, pa)*dot_product(pb,pa)-dot_product(rb-ra,pb))&
    /(1d0-dot_product(pb, pa)**2d0)

Gab=ra+Sab*pa-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))

coll_course=.false.
if (Sab.le.la .and. Sab .ge. 0.0) then
	if (Sba.le.lb .and. Sba .ge. 0.0) then
		coll_course= .true.
	end if
end if

                     end subroutine dist_segs
!************************************************************************
                     
function clamp( a, b, c)
real(8)              ::a, b, c, clamp

if ( a < b ) then
clamp = b
return
end if

if ( a > c ) then
clamp = c
return
end if

clamp = a



end function clamp

subroutine dist_segments(p1,&
                     q1,&
                     p2,&
                     q2,&
                     s,&
                     t,&
                     Gab,&
                     Gab_norm)


real(8), dimension(3):: p1, q1, p2, q2, d1, d2, Gab, r, c1, c2
real(8)              ::Gab_norm, epsilon, s, t
real(8)              ::a, e, f, c, b, denom
logical              :: coll_course

epsilon  =1D-40

d1=q1-p1
d2=q2-p2

r = p1 - p2

a = dot_product(d1,d1)
e = dot_product(d2,d2)
f = dot_product(d2,r)

if ((a<=epsilon) .and. (e<=epsilon)  ) then
    s = 0.0D0
    t = 0.0D0
    
    c1 = p1
    c2 = p2
    
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)
    return
end if
    
if ( a<= epsilon ) then
    s = 0D0
    t = f/e
    t = clamp(t, 0.0D0, 1.0D0)    
else
    c = dot_product(d1,r)
    
    if(e <= epsilon)then
        t= 0.0D0
        s = clamp(-c /a, 0.0D0, 1.0D0)
    else
        b = dot_product(d1, d2)
        denom = a*e-b*b
        
        if(denom /=	0.0D0 ) then
            s = clamp((b*f -c*e) / denom  , 0.0D0, 1.0D0) 
        else
         s =0.0D0
        end if
        
        t =(b*s +f) / e
        
        if (t < 0D0) then
            t = 0.0D0
            s =clamp(-c / a, 0.0D0 , 1.0D0)
        else if ( t> 1.0) then
            
            t= 1.0D0
            s =clamp((b-c)/a, 0.0D0, 1.0D0)
            
       end if
        
        
   end if


    
end if

    c1 = p1 + d1*s
    c2 = p2 + d2*t
    !print *, "c1 " , c1
    !print *, "c2 " , c2
     !print *, "s " , s
     !print *, "t " , t 
     !print *, "p1 " , p1 
     !print *, "p2 " , p2 
     !print *, "d1 " , d1 
     !print *, "d2 " , d2 
     !print *, " "
    Gab_norm = sqrt(dot_product(c1 - c2, c1 - c2))
    Gab =  (c1 - c2)


    end subroutine dist_segments

                     
    subroutine ex_vol_forces_moments_segs2(hingesa1,&
                                      hingesa2,&
                                      hingesb1,&
                                      hingesb2,&
                                      ex_vol_const,&
                                      threshold,&
                                      R,&
                                      Gab_min,&
				                      is_for_neigh_list,&
				                      fric_coeff)

type(rod)                           :: hingesb1,&
                                       hingesb2,&
                                       hingesa1,&
                                       hingesa2
type (vec_arr)       , dimension(7) :: vec
type (simple_vec_arr), dimension(7) ::vec_a
real(8)                             :: threshold, R, ex_vol_const, fric_coeff, Sab, fac,  s, t, Gab_min
integer(8)                          :: i,j, k,l, ind1, ind2
real(8), dimension(3)               :: void, Exc_Vol_Force_Partial, v_rel,  r1, r2, Gab
logical                             :: are_near
logical                             :: is_for_neigh_list
real(8), parameter                      :: pi=3.141592


call dist_segments(hingesa1%X_i,&
               hingesa2%X_i,&
               hingesb1%X_i,&
               hingesb2%X_i,&
                     s,&
                     t,&
                     Gab,&
                     Gab_min)

r1 = (hingesa2%X_i - hingesa1%X_i) 
r2 = (hingesb2%X_i - hingesb1%X_i) 

!print *, " gab_norm1  ", gab_min

if ((is_for_neigh_list.eqv. .false.).and.&
    (Gab_min .lt. threshold) .and. &
     Gab_min .gt. 1e-8   ) then
   

	are_near=.true.
	!Exc_Vol_Force_Partial =-vec(k)%Gab*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))!CORRECT ONE
    !fac = ex_vol_const * pi * hingesb1%ave_viscosity * hingesb1%length * r * hingesb1%gamma_dot
    fac = ex_vol_const*2*r
   Exc_Vol_Force_Partial =-(Gab/(Gab_min))*fac*exp(-2*(Gab_min/R-2))
   
   
    
    !print *, " v1  ", (hingesa2%x_i - hingesa1%x_i)
    !print *, " v2  ", (hingesb2%x_i - hingesb1%x_i)
    !print *, " r1  ", r1
    !print *, " r2  ", r2
    !print *, " s  ", s
    !print *, " t  ", t
    !print *, " (r2*t -r2/2 )  ", (r2*t -r2/2 )   
    !print *, "  "
    
	!Exc_Vol_Force_Partial =-(vec(k)%Gab/(vec(k)%Gab_norm))*2*r*ex_vol_const*exp(-2*(vec(k)%Gab_norm/R-2))
	if (hingesb1%is_stationary==1) hingesb1%v_i=0
	if (hingesa1%is_stationary==1) hingesa1%v_i=0
	
    v_rel= -hingesb1%v_i+hingesa1%v_i -cross(hingesb1%omega,r2*t)&
                                      +cross(hingesa1%omega,r1*s)
    
    if ((v_rel(1) .eq. 0.0) .and. (v_rel(2) .eq. 0.0) .and. (v_rel(3) .eq. 0.0) ) then
        v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
    else
        v_rel =0
    end if
    
    !print *, "v_rel", v_rel
    !print *,"fric_coeff", fric_coeff         
    if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)) ) then !Check for NaN. If NaN, vrel
        Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	end if
	!$OMP ATOMIC
	hingesb1%F_Excl_Vol(1)   = hingesb1%F_Excl_Vol(1)+ Exc_Vol_Force_Partial(1)
    !$OMP ATOMIC
    hingesb1%F_Excl_Vol(2)   = hingesb1%F_Excl_Vol(2)+ Exc_Vol_Force_Partial(2)
    !$OMP ATOMIC
    hingesb1%F_Excl_Vol(3)   = hingesb1%F_Excl_Vol(3)+ Exc_Vol_Force_Partial(3)
#ifdef TENSOR    
    void = cross((r2*t -r2/2 ), Exc_Vol_Force_Partial)
#else
    void = cross((r2*t  ), Exc_Vol_Force_Partial)
#endif
    !$OMP ATOMIC
	hingesb1%T_Excl_Vol(1)   = hingesb1%T_Excl_Vol(1) +void(1)
	!$OMP ATOMIC
    hingesb1%T_Excl_Vol(2)   = hingesb1%T_Excl_Vol(2) +void(2)
    !$OMP ATOMIC
    hingesb1%T_Excl_Vol(3)   = hingesb1%T_Excl_Vol(3) +void(3)
    
    !print *, "collision happened 2"
    !print *, 'fac ', fac
    !print *, 'exc forece ', exc_vol_force_partial
    !!print *," hingesb1%f_excl_vol " , hingesb1%f_excl_vol
    !!print *," hingesb1%t_excl_vol " , hingesb1%t_excl_vol
    !print *, " gab1  ", gab
    !print *, " gab_norm1  ", gab_min
    !print *, " r  ", r
    !print *, " hingesb1%length  ", hingesb1%length
    !print *, " hingesb1%ave_viscosity  ", hingesb1%ave_viscosity
    !print *, " ex_vol_const  ", ex_vol_const
    !print *, " hingesb1%gamma_dot  ", hingesb1%gamma_dot
    !print *, " v_rel ", v_rel
    end if

    
end subroutine ex_vol_forces_moments_segs2

subroutine dist_pt_seg(ra,&
                       rb,&
                       rb_end,&
		               pb,&
                       Gab,&
                       Gab_norm,&
                       Sba,&
                       coll_course)

real(8), dimension(3):: ra, rb, rb_end, pb, Gab
real(8)              :: Sba, Gab_norm, lb
logical              :: coll_course


pb=rb_end-rb

lb=sqrt(dot_product(pb,pb))

pb=pb/lb

Sba=dot_product(ra-rb, pb)

Gab=ra-rb-Sba*pb
Gab_norm=sqrt(dot_product(Gab,Gab))


coll_course=.false.

if (Sba.le.lb .and. Sba .ge. 0.0) then
	coll_course= .true.
end if

end subroutine dist_pt_seg

!************************************************************************
subroutine dist_pt_pt(ra,&
                      rb,&
                      Gab,&
                      Gab_norm)
real(8), dimension(3):: ra, rb, Gab
real(8)              :: Gab_norm


Gab=ra-rb
Gab_norm=sqrt(dot_product (Gab, Gab))

end subroutine dist_pt_pt
!************************************************************************
function are_boxes_intersect(hinge_pt1, hinge_pt2,segment_pt1, segment_pt2, radius, margin)

logical               ::are_boxes_intersect
real(8), dimension(3) ::hinge_pt1, hinge_pt2,segment_pt1, segment_pt2 
type(bound_box)       ::a, b
real(8)               :: radius, margin
are_boxes_intersect=.false.
a.X=(hinge_pt1+hinge_pt2)/2
a.W=abs(hinge_pt1-hinge_pt2)+2*radius+margin


b.X=(segment_pt1+segment_pt2)/2
b.W=abs(segment_pt1-segment_pt2)+2*radius+margin


are_boxes_intersect= (abs(a.X(1) - b.X(1)) * 2 < (a.W(1) + b.W(1))) .and.&
                     (abs(a.X(2) - b.X(2)) * 2 < (a.W(2) + b.W(2))) .and.&
                     (abs(a.X(3) - b.X(3)) * 2 < (a.W(3) + b.W(3)))  

end function are_boxes_intersect
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ex_vol_forces_moments_walls2(fibers,&
                                        hinges,&
                                        r_fiber,&
                                        ex_vol_const,&
                                        box_size,&
                                        fric_coeff,&
                                        is_fric_wall,&
				                        gamma_dot)  
                                        
integer                  :: i, j, k
type(fiber), dimension(:):: fibers
type(rod), dimension(:)  :: hinges
real(8)                  :: ex_vol_const, dist, r_fiber, vec_to_wall, threshold
real(8), dimension(3)    :: box_size, Exc_Vol_Force_Partial, r, v_rel
real(8), dimension(2)    :: wall_position
real(8)                  :: FAC, gamma_dot
logical                  :: is_fric_wall
real(8)                  :: fric_coeff

wall_position(1)= box_size(2)/2
wall_position(2)= -wall_position(1)
threshold         = 1*r_fiber

do k=1,2
    !if (k.eq.1) then
    if (k.eq.2) then
        FAC=+1;
    else
        FAC=-1;
        
    endif

    do i=1, ubound (fibers,1)
    
        j=fibers(i)%first_hinge
        Exc_Vol_Force_Partial=0;
        vec_to_wall= hinges(j)%X_i(2)-wall_position(k)
        dist=abs(vec_to_wall)
        if (dist .le. threshold) then
            r=hinges(j+1)%X_i-hinges(j)%X_i
            Exc_Vol_Force_Partial(2) =FAC*dist*ex_vol_const*exp(-2*(dist-1))
            v_rel=0
            if (hinges(j)%X_i(2).gt.0) then
                v_rel(1)= gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            else
                v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
            end if 
            v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
        
            if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)).and. is_fric_wall ) then !Check for NaN. If NaN, vrel
                 
                Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	        end if
            hinges(j)%F_Excl_Vol   = hinges(j)%F_Excl_Vol+ Exc_Vol_Force_Partial
#ifdef TENSOR            
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(-0.5*r, Exc_Vol_Force_Partial)
#else
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(r, Exc_Vol_Force_Partial)
#endif
         end if

        do j=fibers(i)%first_hinge, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
            Exc_Vol_Force_Partial=0;
            vec_to_wall= hinges(j+1)%X_i(2)-wall_position(k)
            dist=abs(vec_to_wall)
	        r=hinges(j+1)%X_i-hinges(j)%X_i
	        if (dist .le. threshold) then
	            v_rel=0
                if (hinges(j)%X_i(2).gt.0) then
                    v_rel(1)= gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                else
                    v_rel(1)= -gamma_dot*box_size(2)/2-hinges(j)%v_i(1)
                end if 
                v_rel=v_rel/(sqrt(dot_product(v_rel,v_rel)))
                Exc_Vol_Force_Partial(2) =FAC*dist*ex_vol_const*exp(-2*(dist-1))
                
               
                if ((v_rel(1) .eq. v_rel(1)).and.(v_rel(2) .eq. v_rel(2)) .and.(v_rel(3) .eq. v_rel(3)).and. is_fric_wall ) then !Check for NaN. If NaN, vrel
                    Exc_Vol_Force_Partial  = Exc_Vol_Force_Partial+ v_rel * abs(fric_coeff*Exc_Vol_Force_Partial) 
	            end if
	     
       		    hinges(j)%F_Excl_Vol   = hinges(j)%F_Excl_Vol+ Exc_Vol_Force_Partial
#ifdef TENSOR            
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(0.5*r, Exc_Vol_Force_Partial)
#else
            hinges(j)%T_Excl_Vol   = hinges(j)%T_Excl_Vol  +cross(r, Exc_Vol_Force_Partial)
#endif

	        end if  
        end do
    end do
end do

end subroutine ex_vol_forces_moments_walls2

end module excl_volume
