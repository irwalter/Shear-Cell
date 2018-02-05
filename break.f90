module breakage
use data_structures
implicit none
contains


subroutine find_curvature(p1, p2, p3, curv)
implicit none 
real(8) :: a,b,c,l1,l2,l2nb,n,rad,dotp,curv,scale1,scale2,t
real(8)   , dimension(3)  :: center, p1,p2,p3,v1,v2,v1n,v2n,v2nb
real(8)   , dimension(2)  :: p3_2d
integer :: i

center = 0
rad    = 0
v1n    = 0
v2nb   = 0
n = size(p1,1)
v1 = p2 - p1 
v2 = p3 - p1
l1 = sqrt((v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3)))
l2 = sqrt((v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)))
v1n = v1
do i=1,3 
    v1n(i) = v1n(i)/l1
end do
v2n = v2
do i=1,3 
   v2n(i) = v2n(i)/l2
end do
dotp = v2n(1)*v1n(1) + v2n(2)*v1n(2) + v2n(3)*v1n(3)
v2nb = v2n
do i=1,3 
    v2nb(i) = v2nb(i) - dotp*v1n(i)
end do
l2nb = sqrt((v2nb(1)*v2nb(1)+v2nb(2)*v2nb(2)+v2nb(3)*v2nb(3)))
do i=1,3
    v2nb(i) = v2nb(i)/l2nb
end do
do i = 1,2
p3_2d(i) = 0
end do
do i = 1,3
    p3_2d(1) = p3_2d(1) + v2(i)*v1n(i)
    p3_2d(2) = p3_2d(2) + v2(i)*v2nb(i)
end do
a = l1
b = p3_2d(1)
c = p3_2d(2)
t = 0.5*(a-b)/c
scale1 = b/2 + c*t
scale2 = c/2 - b*t
do i = 1,3
center(i) = 0
end do 
do i=1,3
    center(i) = p1(i) + scale1*v1n(i) + scale2*v2nb(i)
end do
rad = sqrt((center(1)-p1(1))**2+(center(2)-p1(2))**2+(center(3)-p1(3))**2)
curv = rad



end subroutine find_curvature

subroutine fiber_damage(fibers, hinges, min_curv, r_fiber)
implicit none
type (fiber), dimension(:), allocatable :: fibers, fibers_temp
type (rod)  , dimension(:), allocatable :: hinges, hinges_temp
real(8)                    :: max_alpha, r_fiber, fac, alpha,min_curv,curv
integer                    :: k, i, j, m, n


fac=0.07*r_fiber

do i=1, ubound(hinges,1)
	hinges(i)%is_broken    =.false.
	hinges(i)%is_separated =.false.
end do

k=0

!print *, "number of fibers" ,  ubound(fibers,1)
do i=1, ubound(fibers,1)
	if(fibers(i)%nbr_hinges.ge.3) then 
		do j=fibers(i)%first_hinge+1, fibers(i)%first_hinge+fibers(i)%nbr_hinges-2
            
            call find_curvature(hinges(j-1)%X_i, hinges(j)%X_i,hinges(j+1)%X_i, curv)
            if (.not. isnan(curv) .and. curv.le.min_curv) then 
                !print *, "curv", curv
 				hinges(j)%is_broken = .true.		 				
				k=k+1
			end if
		end do
	end if
end do
!stop
do i=1, ubound(fibers,1)
	hinges(fibers(i)%first_hinge)%is_separated=.true.
end do

!PRINT *, "Monitor Hinges"
!do i=1, ubound(hinges,1)
!	print *, hinges(i)%X_i, hinges(i)%is_separated, hinges(i)%is_broken
!end do

!PRINT *, "Monitor Fibers"
!do i=1, ubound(fibers,1)
!	print*, fibers(i)%first_hinge, fibers(i)%nbr_hinges
!end do

allocate(hinges_temp(ubound(hinges,1)+k))
allocate(fibers_temp(ubound(fibers,1)+k))
!print *, "Neo Ubound *****", (ubound(fibers,1)+k), (ubound(hinges,1)+k)

k=0

m=1
do i=1, ubound(hinges,1)
	if(hinges(i)%is_broken .eqv. .true.) then
		k=k+1
	end if
	if((hinges(i)%is_broken .eqv. .true.).or.(hinges(i)%is_separated .eqv. .true.)) then
		fibers_temp(m)%first_hinge=k+i
		m=m+1
	end if

end do

do i=1, ubound(fibers_temp,1)-1
	fibers_temp(i)%nbr_hinges=fibers_temp(i+1)%first_hinge-fibers_temp(i)%first_hinge
end do
fibers_temp(ubound(fibers_temp,1))%nbr_hinges=ubound(hinges_temp,1)-fibers_temp(ubound(fibers_temp,1))%first_hinge+1

!do i=1, ubound(fibers_temp,1)
!	print *," NEON FIBERS", fibers_temp(i)%first_hinge
!end do



k=1
do i=1, ubound(hinges,1)
	if (hinges(i)%is_broken .eqv. .true.)then
		!print *, "BREAKAGE"
        
        
        !!!!!!!!!!!!!!!!!!!
        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        !hinges_temp(k) = hinges(i)
        !!!!!!!!!!!!!!!!!!
        
        
        ! updating position to ensure not overlapping
		hinges_temp(k)%X_i= hinges(i)%X_i - fac*(hinges(i)%X_i-hinges(i-1)%X_i)&
                                             /sqrt(dot_product(hinges(i)%X_i-hinges(i-1)%X_i,hinges(i)%X_i-hinges(i-1)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary
		k = k + 1
        
        
        
        !!!!!!!!!!!!!!!!!!
        !adding to alocate all properties of hinges(i) to hinges_temp(i)
        !hinges_temp(k) = hinges(i)
        !!!!!!!!!!!!!!!!!!
        
        
        
        
        ! updating position to ensure not overlapping
		hinges_temp(k)%X_i= hinges(i)%X_i + fac*(hinges(i+1)%X_i-hinges(i)%X_i)&
                                             /sqrt(dot_product( hinges(i+1)%X_i-hinges(i)%X_i, hinges(i+1)%X_i-hinges(i)%X_i))
        hinges_temp(k)%is_stationary = hinges(i)%is_stationary
		k = k + 1
	else
		hinges_temp(k)=hinges(i)
		k = k + 1
	end if

end do


deallocate(hinges)
deallocate(fibers)

allocate(hinges(ubound(hinges_temp,1)))
allocate(fibers(ubound(fibers_temp,1)))

hinges=hinges_temp
fibers=fibers_temp



deallocate(hinges_temp)
deallocate(fibers_temp)
!PRINT *, "Monitor Hinges"
!do i=1, ubound(hinges,1)
!	print *, hinges(i)%X_i
!end do

!PRINT *, "Monitor Fibers"
!do i=1, ubound(fibers,1)
!	print *, fibers(i)%first_hinge, fibers(i)%nbr_hinges
!end do

do i=1, ubound(hinges,1)
	hinges(i)%alpha    =0
end do

end subroutine fiber_damage
end module breakage

