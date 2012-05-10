!~ $define FPPR_MAX_LINE 90
!~ $define FPPR_STP_INDENT 2
!~ $define FPPR_KWD_CASE FPPR_LOWER
!~ $define FPPR_USR_CASE FPPR_LOWER
module globals
  implicit none
  
  type :: simulation_data_t
    ! from command line
    character(255) :: input_file_name
    ! from input
    real :: timestep_len
    real :: subtimestep_length
    real :: p_norm_max
    real :: p_norm_min
    real :: frac_init_runaways
    real :: weight_boundary_in_pth
    real :: runaway_boundary_in_pth
    integer :: number_of_particles
    integer :: number_of_timesteps
    integer :: number_of_flux_surfaces
    integer :: sw_selfcons_electric_field
    integer :: sw_source_term
    integer :: sw_sync
    integer :: averaged_timesteps_in_SC
    integer :: do_maxwellian_test
    integer :: plot_interval
    integer :: quicktest = 0
    logical :: read_from_CPOs = .false.
    ! physical values, which cannot be obtained from CPOs
!    real :: Coulomb_logarithm ! no CPO
    real :: T_e_boundary ! no CPO
    real :: T_e_decay_time ! no CPO
    real :: T_e_min ! no CPO
    real :: dB_per_B ! no cpo data for it, do i have to calculate it in a shot?
    real :: efdbck ! we don't know what is this
    ! not from input
    double precision :: elapsed_time=0.d0
    integer :: calc_status=0
    integer(kind=8) :: mccoll_calls=0
    integer(kind=8) :: mcop_calls=0
  end type simulation_data_t
  
  type :: physical_data_t
    ! from input and the possible CPO
    real :: epsilon_min ! calculated from equilibrium(1)%profiles_1d%rho_tor
    real :: epsilon_max ! calculated from equilibrium(1)%profiles_1d%rho_tor
    real :: epsilon_iron !long term: limiter/limiter_unit, short term: epsilon_max*1.05
    real :: R_major !coreprof(1)%toroid_field%r0
    real :: E_ext_norm2cr !coreprof/profiles1d/vloop <- this is E, not E/Ec!
!    real :: Z_eff !coreprof%profiles_1d%zeff
    real :: B0 !equilibrium%global_param/mag_axis/bphi
    real :: n_e_bulk !coreprof/ne
    real :: n_e_bulk_prof ! calculated from coreprof/ne
    real :: n_e_bulk_prof_exp ! calculated from coreprof/ne
    real :: T_e !coreprof/te
    real :: T_e_bulk_prof ! calculated from coreprof/te
    real :: T_e_bulk_prof_exp ! calculated from coreprof/te

    ! optional input
    real, allocatable :: T_e_numerical(:)
    real, allocatable :: n_e_numerical(:)
    ! not from input
    real :: B2ECR
    real :: RRAT
    real :: TAU
    real :: E_Dreicer
    real :: EBIN
    real :: E_critical
    real :: p_th
    real :: number_of_electrons
    ! not from input, optional
    real, allocatable :: epsilon_numerical(:)
  end type physical_data_t
    
  type(simulation_data_t) :: codeparams
  type(physical_data_t) :: physparams

  ! variables to be used in the codeparams reading routine
  real :: timestep_length
  real :: p_norm_max
  real :: p_norm_min
  real :: frac_init_runaways
  real :: fwb
  real :: rcpn
  integer :: number_of_particles
  integer :: number_of_timesteps
  integer :: number_of_flux_surfaces
  integer :: sw_selfcons_electric_field
  integer :: sw_source_term
  integer :: sw_sync
  integer :: averaged_timesteps_in_SC
  integer :: do_maxwellian_test
  integer :: plot_interval
  integer :: quicktest = 0
  real :: Coulomb_logarithm
  real :: T_e_boundary
  real :: T_e_decay_time
  real :: T_e_min
  real :: dB_per_B
  real :: efdbck

  contains
  subroutine print_global_parameters
    print *, '***DEBUG PRINT OF GLOBALS PARAMS***'
    print *, 'codeparams%input_file_name=', trim(codeparams%input_file_name)

    print *, 'codeparams%timestep_len=', codeparams%timestep_len
    print *, 'codeparams%p_norm_max=', codeparams%p_norm_max
    print *, 'codeparams%p_norm_min=', codeparams%p_norm_min
    print *, 'codeparams%frac_init_runaways=', codeparams%frac_init_runaways
    print *, 'codeparams%weight_boundary_in_pth=', codeparams%weight_boundary_in_pth
    print *, 'codeparams%runaway_boundary_in_pth=', codeparams%runaway_boundary_in_pth
    print *, 'codeparams%number_of_particles=', codeparams%number_of_particles
    print *, 'codeparams%number_of_timesteps=', codeparams%number_of_timesteps
    print *, 'codeparams%number_of_flux_surfaces=', codeparams%number_of_flux_surfaces
    print *, 'codeparams%sw_selfcons_electric_field=', codeparams%sw_selfcons_electric_field
    print *, 'codeparams%sw_source_term=', codeparams%sw_source_term
    print *, 'codeparams%sw_sync=', codeparams%sw_sync
    print *, 'codeparams%averaged_timesteps_in_SC=', codeparams%averaged_timesteps_in_SC
    print *, 'codeparams%do_maxwellian_test=', codeparams%do_maxwellian_test
    print *, 'codeparams%plot_interval=', codeparams%plot_interval
    print *, 'codeparams%quicktest=', codeparams%quicktest
    print *, 'codeparams%read_from_CPOs=', codeparams%read_from_CPOs

!    print *, 'codeparams%Coulomb_logarithm=', codeparams%Coulomb_logarithm
    print *, 'codeparams%T_e_boundary=', codeparams%T_e_boundary
    print *, 'codeparams%T_e_decay_time=', codeparams%T_e_decay_time
    print *, 'codeparams%T_e_min=', codeparams%T_e_min
    print *, 'codeparams%dB_per_B=', codeparams%dB_per_B
    print *, 'codeparams%efdbck=', codeparams%efdbck

    print *, 'physparams%epsilon_min=', physparams%epsilon_min
    print *, 'physparams%epsilon_max=', physparams%epsilon_max
    print *, 'physparams%epsilon_iron=', physparams%epsilon_iron
    print *, 'physparams%R_major=', physparams%R_major
    print *, 'physparams%E_ext_norm2cr=', physparams%E_ext_norm2cr
!    print *, 'physparams%Z_eff=', physparams%Z_eff
    print *, 'physparams%B0=', physparams%B0
    print *, 'physparams%n_e_bulk=', physparams%n_e_bulk
    print *, 'physparams%n_e_bulk_prof=', physparams%n_e_bulk_prof
    print *, 'physparams%n_e_bulk_prof_exp=', physparams%n_e_bulk_prof_exp
    print *, 'physparams%T_e=', physparams%T_e
    print *, 'physparams%T_e_bulk_prof=', physparams%T_e_bulk_prof
    print *, 'physparams%T_e_bulk_prof_exp=', physparams%T_e_bulk_prof_exp
    print *, '***DEBUG END***'

  end subroutine print_global_parameters
end module globals
