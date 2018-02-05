module output
use data_structures
implicit none
contains
subroutine data_out(fibers, hinges, frame, printVelocities)
type (rod), allocatable, dimension(:)  :: hinges
type (fiber), allocatable, dimension(:):: fibers
integer(8)                             :: i, j, k,frame
logical                                :: printVelocities
  open(4,file='OUTPUT/nbr_frames.txt')
  	!write (4,*), ubound(hinges,1)
	write (4,*), frame
  close(4)

  !write (3,*), frame
  k=1
  write (3,*), ubound(fibers,1)
  do i=1, ubound(fibers,1)
  	write (3,*), fibers(i)%nbr_hinges
  	do j=k, k+fibers(i)%nbr_hinges-1
  		write (3,*), hinges(k)%X_i(1), hinges(k)%X_i(2), hinges(k)%X_i(3)
        if (printVelocities) then
            write (5,*), hinges(k)%v_i(1), hinges(k)%v_i(2), hinges(k)%v_i(3)
        end if
        
		k=k+1
  	end do
  end do	

end subroutine data_out






end module output
