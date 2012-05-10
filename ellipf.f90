subroutine ellipf(pm, pe, pk)
!> Elliptic integrals of first and second kind
  implicit none
  integer :: ierr
  real :: pm, pe, pk, rd, rf,  zrf, zrd, zarg
  zarg = 1.0 - pm
  if (zarg .lt. 1.0e-25) zarg = 1.0e-25
  zrf = rf(0.0, zarg, 1.0, ierr)
  zrd = rd(0.0, zarg, 1.0, ierr)

  pe = zrf - 0.33333333 * pm * zrd
  pk = zrf

  return
end