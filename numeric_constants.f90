!~ $define FPPR_MAX_LINE 90
!~ $define FPPR_STP_INDENT 2
!~ $define FPPR_KWD_CASE FPPR_LOWER
!~ $define FPPR_USR_CASE FPPR_LOWER
module numeric_constants
  implicit none
!  real, parameter :: pi =      3.14159265358979323846264338327950288419716939937510
  real :: pi =      3.1416
  real, parameter :: pi_sqrt = 1.77245385090551602729816748334114518279754945612238
  real, parameter :: speed_of_light=3.0e8
  real, parameter :: c0 = 2.99792458e8
  real, parameter :: k_boltz = 1.3806504e-23
  real, parameter :: r_e = 2.8179e-15
  real, parameter :: electron_mass_in_ev=511.17
  real, parameter :: electron_mass = 9.10938188D-31
  real, parameter :: electron_charge = 1.60217646D-19
  real, parameter :: phi_const_z1 =  0.254829592d0
  real, parameter :: phi_const_z2 = -0.284496736d0
  real, parameter :: phi_const_z3 =  1.421413741d0
  real, parameter :: phi_const_z4 = -1.453152027d0
  real, parameter :: phi_const_z5 =  1.061405429d0
  real, parameter :: eps0 = 8.85418782e-12
end module numeric_constants
