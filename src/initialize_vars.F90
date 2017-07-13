! -*-f90-*-
subroutine initialize_vars
  
  use GR1D_module
      
  implicit none

  !don't worry, we call input parser after this to set variable to desired values.

  ntout = 100
  ntout_scalar = 100
  small_output = .false.
  force_restart_output = .false.

  profile_type = 1
  WHW02 = .false.

  dt = 0.0d0
  
  tiny = 1.d-20

  hydro_formulation = "conservative_p_source"
     
  ppm_origin_TVD = 0
      
  reconstruction_method = " "
  M1_reconstruction_method = " "
  tvd_limiter = " "
  flux_type = " "
  fake_neutrinos = .false.
  do_profile_rmax = .false.
  do_nupress = .false.
  do_ye_of_rho = .false.
  do_hydro = .true.
  do_rotation = .false.
  do_leak_ros = .false.
  do_heating = .false.
  vs_mass = .false.

  do_M1 = .false.
  M1_do_backwardfix = 0
  v_order = -1
  M1_maxradii = 0.0d0
  M1_imaxradii = 0
  number_species = 0
  number_species = -1
  number_species_to_evolve = -1
  number_groups = 0
  number_eas = 0
  include_epannihil_kernels = .false.
  include_nes_kernels = .false.
  include_Ielectron_exp = .false.
  include_Ielectron_imp = .false.
  include_energycoupling_exp = .false.
  include_energycoupling_imp = .false.
  opacity_table = ""
  M1closure = 'ME'
  M1_testcase_number = 0
  total_energy_radiated = 0.0d0
  total_energy_absorped = 0.0d0
  total_net_heating = 0.0d0

  shock_radius = 0.0d0      
  ishock(1) = ghosts1+1
  mgravX = 0.0d0
  mgrav12 = 0.0d0
  mbaryX = 0.0d0
  mbary12 = 0.0d0
  rXmax = 0.0d0
  r12max = 0.0d0
  r11max = 0.0d0
  mass_inner_core = 0.0d0
  heat_fac = 0.0d0

  angular_momentum = 0.0d0
  rotational_energy = 0.0d0
  internal_energy = 0.0d0
  omega_c = 0.0d0
  omega_A = 1.0d0
  ToverW(:) = 0.0d0

  eoskey = 0
  
  initial_data = " "

  geometry = 0
  
  gridtype = " "
  
  iorder_hydro = 2
  
  bounce = .false.
  t_bounce = 0.0d0
  
  binding_energy_total = 0.0d0
  binding_energy_envelope = 0.0d0

  dynamic_output_control = .true.
  gravity_active = .false.
  
  grid_rmax = 0.0d0
  rho_cut = 0.0d0
  grid_custom_inner = 0.0d0
  grid_custom_dx1 = 0.0d0
  grid_custom_rad1 = 0.0d0
  grid_custom_number = 0
  
  call initialize_arrays

end subroutine initialize_vars

subroutine initialize_arrays

  use GR1D_module
  
  implicit none

  mass(:) = 0.0d0
  mass1(:) = 0.0d0
  volume(:) = 0.0d0
  
  x1(:) = 0.0d0
  x1i(:) = 0.0d0
  
  cs2(:) = 0.0d0
  cs2p(:) = 0.0d0
  cs2m(:) = 0.0d0

  rho(:) = 0.0d0
  rhop(:) = 0.0d0
  rhom(:) = 0.0d0
  
  eps(:) = 0.0d0
  epsp(:) = 0.0d0
  epsm(:) = 0.0d0

  energy_nu(:) = 0.0d0
  mom_nu(:) = 0.0d0
  
  eps_kin(:) = 0.0d0
  binding_energy(:) = 0.0d0

  accretion_radii(:) = ghosts1+1
  accretion_rates(:) = 0.0d0
  accreted_mass(:) = 0.0d0

  press(:) = 0.0d0
  pressth(:) = 0.0d0
  pressp(:) = 0.0d0
  pressm(:) = 0.0d0
  press_nu(:) = 0.0d0
  dnupdr(:) = 0.0d0
  
  ye(:) = 0.0d0
  ye_prev(:) = 0.0d0
  dyedt_hydro(:) = 0.0d0
  dyedt_neutrino(:) = 0.0d0
  yep(:) = 0.0d0
  yem(:) = 0.0d0
  
  ynu(:) = 0.0d0

  v1(:) = 0.0d0
  v1p(:) = 0.0d0
  v1m(:) = 0.00d0

  if(do_rotation) then
     vphi1(:) = 0.0d0
     vphi1p(:) = 0.0d0
     vphi1m(:) = 0.0d0
     vphi(:) = 0.0d0
     vphip(:) = 0.0d0
     vphim(:) = 0.0d0
     omega(:) = 0.0d0
  endif

  v(:) = 0.0d0
  vp(:) = 0.0d0
  vm(:) = 0.0d0
  v_prev(:) = 0.0d0
  
  q(:,:) = 0.0d0
  qold(:,:) = 0.0d0
  qm(:,:) = 0.0d0
  qp(:,:) = 0.0d0
  
  q_hat_old(:,:) = 0.0d0
  q_hat(:,:) = 0.0d0

  eas(:,:,:,:) = 0.0d0
  ies(:,:,:,:,:) = 0.0d0
  ies_sourceterm(:,:,:,:) = 0.0d0
  epannihil(:,:,:,:,:) = 0.0d0
  epannihil_sourceterm(:,:,:,:) = 0.0d0
  q_M1(:,:,:,1) = 1.0d-100
  q_M1(:,:,:,2) = 1.0d-110
  q_M1(:,:,:,3) = 0.5
  q_M1p(:,:,:,:,:) = 0.0d0
  q_M1m(:,:,:,:,:) = 0.0d0
  q_M1_extra(:,:,:,:) = 1.0d-90
  q_M1_extrap(:,:,:,:,:) = 0.0d0
  q_M1_extram(:,:,:,:,:) = 0.0d0
  q_M1_fluid(:,:,:,1) = 1.0d-90
  q_M1_fluid(:,:,:,2) = 1.0d-110
  M1_matter_source(:,:) = 0.0d0
  M1_moment_to_distro(:) = 0.0d0
  M1_moment_to_distro_inverse(:) = 0.0d0

  nuchem(:) = 0.0d0
  elechem(:) = 0.0d0
  massfrac_n(:) = 0.0d0
  massfrac_p(:) = 0.0d0
  massfrac_a(:) = 0.0d0
  massfrac_h(:) = 0.0d0
  massfrac_abar(:) = 0.0d0
  massfrac_zbar(:) = 0.0d0

  eas(:,:,:,:) = 0.0d0
  ies(:,:,:,:,:) = 0.0d0
  ies_sourceterm(:,:,:,:) = 0.0d0
  epannihil(:,:,:,:,:) = 0.0d0
  epannihil_sourceterm(:,:,:,:) = 0.0d0
  q_M1(:,:,:,1) = 1.0d-100
  q_M1(:,:,:,2) = 1.0d-110
  q_M1(:,:,:,3) = 0.5
  q_M1p(:,:,:,:,:) = 0.0d0
  q_M1m(:,:,:,:,:) = 0.0d0
  q_M1_extra(:,:,:,:) = 1.0d-90
  q_M1_extrap(:,:,:,:,:) = 0.0d0
  q_M1_extram(:,:,:,:,:) = 0.0d0
  q_M1_fluid(:,:,:,1) = 1.0d-90
  q_M1_fluid(:,:,:,2) = 1.0d-110
  M1_matter_source(:,:) = 0.0d0
  M1_moment_to_distro(:) = 0.0d0
  M1_moment_to_distro_inverse(:) = 0.0d0

  sqrt_gamma(:) = 0.0d0
  
  flux_diff(:,:) = 0.0d0
  gravsource(:,:) = 0.0d0
  presssource(:,:) = 0.0d0
  coolingsource(:,:) = 0.0d0
  denergyloss(:) = 0.0d0

  atmo(:) = 0
!##########################################
! GR VARIABLES
  phi(:) = 0.0d0
  phii(:) = 0.0d0
  alp(:) = exp(phi(:))
  alpp(:) = exp(phi(:))
  alpm(:) = exp(phi(:))
  X(:) = 1.0d0
  Xp(:) = 1.0d0
  Xm(:) = 1.0d0
  W(:) = 1.0d0
  Wp(:) = 1.0d0
  Wm(:) = 1.0d0
  mgrav(:) = 0.0d0
  mgravi(:) = 0.0d0

end subroutine initialize_arrays
  

