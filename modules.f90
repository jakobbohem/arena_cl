!     
! File:   modules.f90
! Author: jakob
!
! Created on 25 April 2012, 14:27
!

!> module for storing elliptic integrals, calculated by ellipf
module comek
  implicit none
  !> number of points at which elliptic integrals are stored
  integer, parameter :: max_num_of_elliptic_integral_points = 502
  real :: ef(max_num_of_elliptic_integral_points)
  real :: kf(max_num_of_elliptic_integral_points)
end module comek