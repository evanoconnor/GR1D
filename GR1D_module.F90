! -*-f90-*-
module GR1D_module

  character*80 :: input_file = "parameters"

  ! parameters to be read
  character*128 jobname
  character*128 outdir
  character*128 profile_name
  character*256 yeprofile_name
  character*256 opacity_table     
  character*256 eos_table_name
  character*256 restart_file_name
  logical force_restart_output
  logical wipe_outdir
  logical do_restart

  integer profile_type   !1 for .short
  logical WHW02

  !grid variables
  integer n1
  integer radial_zones
  real*8 grid_custom_rad1
  real*8 grid_custom_dx1
  real*8 grid_custom_inner
  integer grid_custom_number
  real*8 grid_rmax
  real*8 smallest_dx
  logical do_profile_rmax
  real*8 rho_cut

  !eos desired precision
  real*8,parameter :: eos_rf_prec = 1.0d-9

  !number of ghost zones
  integer ghosts1

  real*8 time,timep,timepp,time_c
  real*8 dt,dtp,dtpp
  real*8 :: dt_max = 1.0d0
  real*8 :: dt_min = 1.0d-10
  real*8 :: dt_reduction_factor = 1.0d0

  !logicals for output handling
  logical :: OutputFlag = .false.
  logical :: CheckFlag = .false.
  logical :: OutputFlagScalar  = .false.
  logical :: OutputFlagRestart = .false.

  integer nt,ntmax
  !output interval settings
  integer ntout,ntout_scalar,ntinfo,ntout_restart
  real*8 dtout,dtout_scalar,dtout_restart
  !original settings from parameter file
  integer ntout0,ntout_scalar0,ntinfo0,ntout_restart0
  real*8 dtout0,dtout_scalar0,dtout_restart0
  logical small_output

  real*8 tdump,tdump_scalar,tdump_restart

  ! status variable that tells us in what RK step we are in
  integer :: rkstep = 0

  real*8 tend

  real*8 cffac
  real*8 cffac_original
  real*8 tiny

  real*8 t_bounce

  real*8 totalmass
  real*8 binding_energy_total
  real*8 binding_energy_envelope

  integer accretion_radii(11)
  real*8 accretion_rates(11)
  real*8 accreted_mass(11)

  real*8 hybridgamma_th
  real*8 hybridgamma1
  real*8 hybridgamma2

  real*8 atmo_rho_abs_min
  real*8 atmo_rho_rel_min
  real*8 atmo_fac

  integer n_cons

  ! flags
  integer eoskey
  integer iorder_hydro
  character*80 hydro_formulation
  character*80 reconstruction_method
  character*80 M1_reconstruction_method

  character*80 flux_type

  character*80 tvd_limiter

  character*80 initial_data

  character*80 gridtype

  logical gravity_active

  logical dynamic_output_control
  logical vs_mass

  integer geometry
  integer ppm_origin_TVD

  logical :: do_leak_ros = .false.
  logical :: do_heating = .false.
  logical :: do_NNBrem = .false.
  logical :: do_hydro = .true.
  logical :: fake_neutrinos = .false.
  logical :: do_nupress = .false.
  logical :: do_ye_of_rho = .false.
  logical :: do_yeofrhofit = .false.
  logical :: accretion_analysis = .true.
  logical :: bounce = .false.
  logical :: switch1 = .false.
  logical :: switch2 = .false.
  logical :: switch3 = .false.

  integer ishock(1)
  real*8 shock_radius
  real*8 gain_radius
  integer igain(1)
  real*8 mgravX,mgrav12
  real*8 mbaryX,mbary12
  real*8 mass_inner_core
  real*8 rXmax,r12max,r11max
  real*8 heat_fac

  !rotation
  logical :: do_rotation = .false.
  logical :: set_omega = .false.
  real*8 omega_c
  real*8 omega_A
  real*8 angular_momentum
  real*8 rotational_energy
  real*8 internal_energy
  real*8,allocatable,save :: ToverW(:)

  !testcases variables
  integer :: shocktube_problem

  !M1 scheme stuff
  logical :: do_M1 = .false.
  integer :: M1_prev_phase = 1

  integer :: v_order
  real*8 :: M1_maxradii
  real*8 :: M1_extractradii
  integer :: M1_imaxradii
  integer :: M1_iextractradii
  integer :: number_species
  integer :: number_species_to_evolve
  integer :: number_groups
  integer :: number_eas
  logical :: include_epannihil_kernels
  logical :: include_bremsstrahlung_kernels
  logical :: include_nes_kernels
  logical :: include_Ielectron_exp
  logical :: include_Ielectron_imp
  logical :: include_energycoupling_exp
  logical :: include_energycoupling_imp

  character*80 :: M1closure
  integer :: M1_testcase_number
  real*8 :: total_energy_radiated
  real*8 :: total_energy_absorped
  real*8 :: total_net_heating
  real*8 :: total_net_deintdt
  real*8 :: total_mass_gain

  !M1 controls
  integer :: M1_do_backwardfix
  logical :: M1_control_flag = .false.
  real*8 :: M1_phase1phase2_density
  real*8 :: M1_phase2phase3_pbtime
  character*80 M1_phase1_reconstruction
  character*80 M1_phase2_reconstruction
  character*80 M1_phase3_reconstruction
  real*8 :: M1_phase1_cffac
  real*8 :: M1_phase2_cffac
  real*8 :: M1_phase3_cffac
  integer :: M1_phase1_ns
  integer :: M1_phase2_ns
  integer :: M1_phase3_ns
  integer :: M1_phase1_ies_way
  integer :: M1_phase2_ies_way
  integer :: M1_phase3_ies_way
  integer :: M1_phase1_encpl_way 
  integer :: M1_phase2_encpl_way
  integer :: M1_phase3_encpl_way


  ! all evolution variables: 3 timelevels; p suffix: previous timelevels
  ! pp suffix: previous previous timelevel

  ! #######################################################
  ! GR VARIABLES
  logical :: GR = .false.
  logical :: do_effectivepotential = .false.
  real*8,allocatable,save :: phi(:),phii(:),alp(:),X(:),W(:), &
       alpp(:),alpm(:),Xp(:),Xm(:),Wp(:),Wm(:)
  real*8,allocatable,save :: dphidr(:),mgrav(:),mgravi(:)

  ! romero's velocity
  real*8,allocatable,save :: v(:),vm(:),vp(:)
  real*8,allocatable,save :: v_prev(:)
  real*8,allocatable,save :: vphi(:),vphim(:),vphip(:)

  ! #######################################################

  ! radial coordinate 
  real*8,allocatable,save :: x1(:)
  ! inner interface coordinate
  real*8,allocatable,save :: x1i(:)
  ! density, cell center and at interfaces
  real*8,allocatable,save :: rho(:), rhop(:),rhom(:)
  ! specific internal energy, cell center and at interfaces
  real*8,allocatable,save :: eps(:),epsp(:),epsm(:)
  real*8,allocatable,save :: energy_nu(:)
  real*8,allocatable,save :: mom_nu(:)
  ! mass interior cell center, mass of cell, volume of cell
  real*8,allocatable,save :: mass(:), mass1(:), volume(:)
  ! pressure 
  real*8,allocatable,save :: press(:),pressm(:),pressp(:)
  real*8,allocatable,save :: press_nu(:),dnupdr(:)
  ! speed of sound squared
  real*8,allocatable,save :: cs2(:), cs2m(:), cs2p(:)
  ! temperature 
  real*8,allocatable,save :: temp(:)
  ! entropy 
  real*8,allocatable,save :: entropy(:)
  ! velocity & rotation
  real*8,allocatable,save :: v1(:),v1m(:),v1p(:)
  real*8,allocatable,save :: vphi1(:),vphi1m(:),vphi1p(:)
  real*8,allocatable,save :: omega(:)
  ! ye 
  real*8,allocatable,save :: ye(:),yem(:),yep(:),ye_prev(:)
  real*8,allocatable,save :: dyedt_hydro(:),dyedt_neutrino(:)
  real*8,allocatable,save :: depsdt(:),dyedt(:)
  real*8,allocatable,save :: ynu(:)
  !chemical potential
  real*8,allocatable,save :: nuchem(:),elechem(:)
  !mass fractions
  real*8,allocatable,save :: massfrac_p(:)
  real*8,allocatable,save :: massfrac_n(:)
  real*8,allocatable,save :: massfrac_a(:)
  real*8,allocatable,save :: massfrac_h(:)
  real*8,allocatable,save :: massfrac_abar(:)
  real*8,allocatable,save :: massfrac_zbar(:)

  ! conserved variables
  real*8,allocatable,save :: q(:,:),qm(:,:),qp(:,:)
  real*8,allocatable,save :: qold(:,:)
  real*8,allocatable,save :: q_hat_old(:,:)
  real*8,allocatable,save :: q_hat(:,:)

  !M1 conserved variables and such
  real*8,allocatable,save :: q_M1(:,:,:,:) !energy density, flux,
  !and pressure tensor
  !(for each species and
  !energy)
  real*8,allocatable,save :: q_M1_prev(:,:,:,:) !energy density, flux,
  !and pressure tensor
  !(for each species and
  !energy)
  real*8,allocatable,save :: q_M1p(:,:,:,:,:) !reconstructed, plus
  !interfaces
  real*8,allocatable,save :: q_M1m(:,:,:,:,:) !reconstructed, minus
  !interfaces

  real*8,allocatable,save :: q_M1_extra(:,:,:,:) !extra moments
  !from the closure
  !(P_tt, Q_rrr,
  !Q_ffr, for
  !example)
  real*8,allocatable,save :: q_M1_extrap(:,:,:,:,:) !extra moments
  !reconstructed,
  !plus
  real*8,allocatable,save :: q_M1_extram(:,:,:,:,:) !extra moments
  !reconstructed,
  !minus
  real*8,allocatable,save :: q_M1_old(:,:,:,:) !energy density, flux,
  !and pressure tensor
  !(for each species and
  !energy)
  real*8,allocatable,save :: B_M1(:,:,:,:) !variables used to store
  !explicit parts of IMEX
  !update
  real*8,allocatable,save :: C_M1(:,:,:,:) !variables used to store
  !explicit parts of IMEX
  !update
  real*8,allocatable,save :: D_M1(:,:,:,:) !variables used to store
  !explicit parts of IMEX
  !update
  real*8,allocatable,save :: flux_M1(:,:,:,:) !used to store
  !calculated flux
  !temporarly
  real*8,allocatable,save :: flux_M1_energy(:,:,:,:) !used to store
  !calculated flux
  !temporarly
  real*8,allocatable,save :: flux_M1_scatter(:,:,:,:) !used to store
  !calculated flux
  !temporarly
  real*8,allocatable,save :: q_M1_fluid(:,:,:,:) !energy density, flux,
  !and pressure tensor
  !(for each species and
  !energy)
  real*8,allocatable,save :: eas(:,:,:,:) !emissivity,absorbive
  !crosssection, scattering
  !crosssection (for each
  !species and energy)
  real*8,allocatable,save :: ies(:,:,:,:,:) !inelastic kernals for
  !electron scattering
  !zeroth and first moment
  real*8,allocatable,save :: ies_sourceterm(:,:,:,:) !for keeping
  !track of
  !source term
  !from ies for
  !matter
  real*8,allocatable,save :: epannihil(:,:,:,:,:) !kernals for
  !ep-annihilation zeroth
  !and first moment
  real*8,allocatable,save :: epannihil_sourceterm(:,:,:,:) !for keeping
  !track of
  !source term
  !from epannihil for
  !matter
  real*8,allocatable,save :: bremsstrahlung(:,:,:,:,:) !kernals for
  !brem zeroth
  !and first moment
  real*8,allocatable,save :: bremsstrahlung_sourceterm(:,:,:,:) !for keeping
  !track of
  !source term
  !from brem for
  !matter
  real*8,allocatable,save :: M1_matter_source(:,:) !source terms
  !for matter
  !interaction
  real*8,allocatable,save :: M1_moment_to_distro(:) !conversion factor
  real*8,allocatable,save :: M1_moment_to_distro_inverse(:) !inverse conversion factor

  ! metric stuff
  real*8,allocatable,save :: sqrt_gamma(:)

  ! flux differences
  real*8,allocatable,save :: flux_diff(:,:)
  ! gravity sources
  real*8,allocatable,save :: gravsource(:,:)
  ! pressure source terms
  real*8,allocatable,save :: presssource(:,:)
  real*8,allocatable,save :: pressth(:)
  ! cooling source terms
  real*8,allocatable,save :: coolingsource(:,:)
  real*8,allocatable,save :: denergyloss(:)
  !atmosphere
  integer,allocatable,save :: atmo(:)
  ! analysis
  real*8,allocatable,save :: eps_kin(:)
  real*8,allocatable,save :: binding_energy(:)

  ! constants
  real*8,parameter :: pi = 3.14159265358979d0
  real*8,parameter :: pi4 = 12.5663706143592d0
  real*8,parameter :: ggrav = 6.673d-8
  real*8,parameter :: clite = 2.99792458d10
  real*8,parameter :: clite_g = 2.99792458d10

  real*8,parameter :: rho_gf = 1.61930347d-18
  real*8,parameter :: press_gf = 1.80171810d-39
  real*8,parameter :: eps_gf = 1.11265006d-21
  real*8,parameter :: time_gf = 2.03001708d+05
  real*8,parameter :: mass_gf = 5.02765209d-34
  real*8,parameter :: length_gf = 6.77140812d-06
  real*8,parameter :: energy_gf = 5.59424238d-55
  real*8,parameter :: lum_gf = 2.7556091d-60

  real*8,parameter :: mev_to_erg = 1.60217733d-6
  real*8,parameter :: erg_to_mev = 6.24150636d5
  real*8,parameter :: amu_cgs = 1.66053873d-24
  real*8,parameter :: massn_cgs = 1.674927211d-24
  real*8,parameter :: amu_mev = 931.49432d0
  real*8,parameter :: kb_erg = 1.380658d-16
  real*8,parameter :: kb_mev = 8.61738568d-11
  real*8,parameter :: temp_mev_to_kelvin = 1.1604447522806d10
  real*8,parameter :: planck = 6.626176d-27
  real*8,parameter :: avo = 6.0221367d23
  real*8,parameter :: hbarc_mevcm = 1.97326966d-11

  real*8,parameter :: nulib_emissivity_gf = 5.59424238d-55/(6.77140812d-06**3*2.03001708d+05)
  real*8,parameter :: nulib_opacity_gf = 1.0d0/6.77140812d-06
  real*8,parameter :: nulib_energy_gf = 1.60217733d-6*5.59424238d-55
  real*8,parameter :: nulib_kernel_gf = 6.77140812d-06**3/2.03001708d+05

  real*8,parameter :: twothirds = 2.0d0/3.0d0
  real*8,parameter :: onethird = 1.0d0/3.0d0

end module GR1D_module
