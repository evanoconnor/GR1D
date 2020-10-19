! -*-f90-*-
subroutine start

  use GR1D_module
  use ppm
#if HAVE_LEAK_ROS
  use leakage_rosswog,only: initialize_leakage_rosswog
#endif
  implicit none

  character(len=128) cpstring
  character(len=128) rmstring
  logical :: outdirthere

  write(*,*) "Parsing input file..."
  !this time to get details of arrays sizes etc
  call input_parser

  if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     ! Ott EOS routines: read table
     call readtable(eos_table_name)
#else
     stop "start.F90: NUC_EOS code not present but want eoskey=3"
#endif
  endif

  !total zones
  n1 = radial_zones+ghosts1*2

  if(do_rotation) then
     n_cons = 5
  else
     n_cons = 4
  endif

  !allocate & initialize variables
  call allocate_vars
  call initialize_vars

  !this time to set all variables requested values
  call input_parser
  write(*,*) "Your job name: ", trim(adjustl(jobname))
  write(*,*) "Using ", radial_zones, " radial zones and ", ghosts1, " ghosts zones"
  write(*,*) "Using ",trim(adjustl(reconstruction_method))," reconstruction"
  write(*,*) "Using ",trim(adjustl(flux_type)), " for hydro"

  ! check if output directory exists
#if __INTEL_COMPILER
  inquire(directory=trim(adjustl(outdir)),exist=outdirthere)
#else
  inquire(file=trim(adjustl(outdir)),exist=outdirthere)
#endif
  if(.not.outdirthere) then
     write(6,*) "*** Output directory does not exist! Please mkdir the output directory!"
     write(6,*) trim(adjustl(outdir))
     stop "Aborting!"
  endif

  ! wipe output directory
!~   rmstring="rm -rf "//trim(adjustl(outdir))//"/*"
!~   call system(rmstring)
  ! copy parameter file
  cpstring="cp parameters "//trim(adjustl(outdir))
  call system(cpstring)

  !setting up initial data
  call problem

  !Collapse specific setups
  if(initial_data.eq."Collapse") then
     write(*,*) "Finding envelope binding energy"
     call map_envelope_binding_energy(profile_name)
     write(*,*) "Finding accretion radii"
     call findaccretionradii
  endif

  if(reconstruction_method.eq.'ppm'.or.M1_phase1_reconstruction.eq.'ppm' &
       .or.M1_phase2_reconstruction.eq.'ppm'.or.M1_phase3_reconstruction.eq.'ppm') then
     write(*,*) "Setting up PPM coefficients"
     call ppm_coefficients
  endif

  !setup dynamic output control
  if(dynamic_output_control) then
     call output_control
  endif

  if(do_M1) then
     write(*,*) "Setting up M1 transport"
     if (M1_control_flag) then
        call M1_control
     endif
     call M1_init
  endif

  !if leaking, initialize variables etc.
#if HAVE_LEAK_ROS
  call initialize_leakage_rosswog
#endif

  !setup dumping variables
  tdump_scalar = dtout_scalar
  tdump = dtout
  tdump_restart = 0.0d0 !make a restart file at the beginning
  
  time = 0.0d0
  time_c = 0.0d0
  nt = 0

  OutputFlag = .false.
  OutputFlagScalar = .false.

  !if doing restart here is where we call to read in file.
  if (do_restart) then
     write(*,*) "Doing restart: ", restart_file_name
     call restart_init_h5
     !warning, this will adjust out put times
     if (force_restart_output) then
        tdump_restart = time+dtout_restart
        tdump = time
        tdump_scalar = time
     endif
     if (M1_control_flag) then
        call M1_control
     endif
     call restart_output_h5
  endif


end subroutine start
