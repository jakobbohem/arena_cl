!~ $define FPPR_MAX_LINE 90
!~ $define FPPR_STP_INDENT 2
!~ $define FPPR_KWD_CASE FPPR_LOWER
!~ $define FPPR_USR_CASE FPPR_LOWER
module collision_operator
    use numeric_constants
!    use globals
!    use comdis
    implicit none
    
    contains
!DEC$ ATTRIBUTES FORCEINLINE :: Phi
    real function phi (x, exp_x2)
      implicit none
      real :: x, zt, exp_x2

      zt = 1.0d0 / (1.0d0+0.3275911d0*x)
      phi = 1.0d0 - ((((phi_const_z5*zt+phi_const_z4)*zt+phi_const_z3)*zt+&
     & phi_const_z2)*zt+phi_const_z1) * zt * exp_x2
      return
    end function phi
!
!DEC$ ATTRIBUTES FORCEINLINE :: Psi
    real function psi (x, x2, exp_x2)
      real :: x, x2, exp_x2
      psi = (phi (x, exp_x2)-x*1.128379167D0*exp_x2) / (2.0*x2)
      return
    end function psi
!
!DEC$ ATTRIBUTES FORCEINLINE :: Psi_prime
    real function psi_prime (x, x2, exp_x2)
      real :: x, x2, x3, phip, phipp, exp_x2
      x3 = x * x2
      phip = 2.0 / pi_sqrt * exp_x2
      phipp = - 4.0 * x / pi_sqrt * exp_x2
      psi_prime = 0.5 * (-phipp*x2-2.0*(phi(x, exp_x2)-x*phip)) / x3
      return
    end function psi_prime
!
    subroutine mccoll_v1 (peps, pp, plam, pi1, pi2, ppcdot, pacp, plcdot, paclam, temperature, zeff, lnLam)
!      use globals
!      use comdis
      implicit none
      real, intent (in) :: peps, pp, plam, pi1, pi2
      real, intent (in) :: temperature, lnLam, zeff
      real, intent (out) :: ppcdot, pacp, plcdot, paclam
      real :: pacp1, pacp2
      real :: ppcdot1, ppcdot2
      real :: zbeta, zdnu, znu
      real :: zp, zp2, zp3, zp4
      real :: zpb, zpb2, zpb3, zpb4
      real :: zpc, zpc4, zpcpa
      real :: zpth, zpth2
      real :: zzpth, zzpth2
      real :: zt1, zt2, zt3, zt4
      real :: ztheta
      real :: zz, tempe
      real :: beta, nu_d
      real :: zx, zx2, zg, zve, znuve, zb1, zb2, zalpha, exp_x2
      
      ! logisv and logder were logicals in the merljin and lge versions since their
      ! first letter is 'L'. However they are integers here, because Gfortran 
      ! has problems of converting real-to-logical, but can do real-to-integer conversion.
      ! This is sufficient here since their value can be in the interval [0;1],
      ! which upon conversion will be 0 or 1, which are equialent to logical
      ! .false. or .true., so in this special case using logical or integer gives 
      ! the same results
      integer :: logisv, logder

!      mccoll_calls = mccoll_calls + 1 
!
!      temperature=tempe(peps)
      !MJ: normalised thermal velocity
      zpth = sqrt (temperature/electron_mass_in_ev)
      !MJ: epsilon (norm. therm. vel. squared)
      zpth2 = temperature / electron_mass_in_ev
      !more modern definition of therm. vel.
      zzpth = zpth * sqrt (2.0)
!
      zp = pp + 1.0E-6
      zp2 = zp * zp
      zp3 = zp2 * zp
      zp4 = zp3 * zp
!
      zpc = 0.4 * zpth
      zpc4 = zpc * zpc * zpc * zpc
      zpcpa = 2.0 * zpth
!
      if (zp .lt. zpcpa) then
        zpb = zpcpa
        zpb2 = zpcpa * zpcpa
        zpb3 = zpb2 * zpcpa
        zpb4 = zpb3 * zpcpa
      else
        zpb2 = zp * zp
        zpb3 = zp2 * zp
        zpb4 = zp3 * zp
      end if
!
      ztheta = 1.9563E-3 * temperature
!
      zt1 = (1.0+zp2)
      zt2 = (1.0+log(zp+1.0)/lnLam)
      zt3 = zpc4 + zp4
      zt4 = sqrt (zt1)
!
      zz = ztheta * (1.0+2.0*zpb2) / (zpb2*(1.0+zpb2))
      zz = min (zz, 1.0)
!
      zbeta = 2.0 * (1.0+zeff-zz) * sqrt (1.0+zpb2) / zpb3 * (1.0+log(zp+1.0)/lnLam)

      znu = zt1 * zt2 / zp3
      znu = znu * zp4 / zt3
      zdnu = (1.0+3.0*zp2-4.0*(1.0+zp2)*zp4/zt3) / zt3
!
      ppcdot1 = - zp * znu + ztheta * &
     & (((2.0*znu/zp+zdnu)*zt4+znu*zp/zt4)*zt2+znu*zt4/((zp+1.0)*lnLam))
!
      pacp1 = sqrt (2.0*znu*ztheta*zt4)
      

! OUTPUT LAMBDA
      plcdot = - zbeta * (0.5*plam-pi2/pi1)
      paclam = sqrt (2.0*zbeta*plam*pi2/pi1)
      if (pi1 == 0) then
        print *, "WARNING: dividing by 0 will yeild an infinite quantity (NaN)"
        print *, "DBG: paclam = sqrt(2*zbeta*plam*pi2/pi1)"
        print *, "DBG: paclam = sqrt(2*",zbeta,"*",plam,"*",pi2,"/",pi1,") =",paclam
      end if
!
! MJ  LOW ENERGY COLL OP.
      zx = zp / zzpth
      zx2 = zx * zx
      zg = gstix (zx)
      zzpth2 = zzpth * zzpth
      zve = speed_of_light * zzpth / sqrt (1.0+zzpth2)
      znuve = speed_of_light * speed_of_light * speed_of_light / (zve*zve*zve)
      exp_x2 = exp (-zx2)
! 
      pacp2 = sqrt (2.0*znuve*zzpth2*zg/zx)
      ppcdot2 = znuve * zzpth * (-zg*(2.0+1.0/zx2)+1.1248*exp_x2/zx)
!
! MJ  THE DERRIVATIVE PART OF L
      zalpha = zp / sqrt (2*zpth2*(1+zp2))
      zb1 = zpth2 * zt1 * zt4 / zp3
      zb2 = zpth2 * zt1 * zt4 * gstix (zalpha) / zp
!
      !logistics function evaluation to speed things up
      logisv = logistic (zp, temperature)
      !derrivative of logistics function
      logder = logisv * (1-logisv)
! OUTPUT P
      pacp = pacp1 * logisv + pacp2 * (1-logisv)
      ppcdot = ppcdot1 * logisv + ppcdot2 * (1-logisv) + zb1 * logder - zb2 * logder
!      print *, "pacp = pacp1 * logisv + pacp2 * (1-logisv)"
!      print *, "pacp = ",pacp1," * ",logisv," + ",pacp2," * (1-",logisv,") =",pacp
      
      
      return
    end subroutine mccoll_v1

!>  Function G
!>  Purpose: calculate gstix(x)=(erf(x) - x * erf'(x) / (2 * x**2)
    real function gstix (x)
      implicit none
      real :: x,x2,exp_x2
      x2=x*x
      if (x .lt. 0.01) then
        gstix = 1.1284 * (0.33333D0*x+0.2D0*x2*x)
      else
        if (x2 .lt. 20.0D0) then
          gstix=psi (x, x2, exp(-x2))
        else
          gstix = 1.0D0 / (2.0D0*x2)
        end if
      end if
      return
    end function gstix

!>  Function LOGISTIC
!  PURPOSE: CALCULATE logistic(x) = 1/(1 + exp(-(x-x0)/deltax))
!    x0     = centering parameter, to be set by the user
!    deltax = width parameter, to be set by the user
    real function logistic (px, temperature)
      implicit none
      real :: px, px0, pdx, temperature

      !MJ: treshold for transition between low and high energy region
      px0 = 3 * sqrt (temperature/electron_mass_in_ev)
      pdx = 0.1 * px0

      logistic = 1 / (1+exp(-(px-px0)/pdx))
      return
    end function logistic
end module collision_operator
