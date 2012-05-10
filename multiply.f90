program multiply
  use cl_collop
  implicit none
  
  integer :: n, k, i
  integer :: sout = 4 ! length of output array
  integer, dimension(:), allocatable :: ns
  integer, dimension(:,:), allocatable :: times
  real, dimension(4) :: outp_times
  
  ! sanity check
  call runClCompare(2,outp_times)
  
!  allocate(ns(2))
!  allocate(times(sout,2))
!  ns = (/ 2, 5 /)
!  n = 2
!  
!  ! this includes the 'cl' library
!  do k = 1,size(ns)
!    call runClCompare(ns(k), outp_times)
!    do i = 1,sout
!      times(k,i) = outp_times(i);
!    end do
!  end do
!  ! measured times output
!  do k = 1,sout
!    print *, "measured time: ", times(k,2)
!  end do
! 
!  deallocate(ns)
!  deallocate(times)
  
end program multiply