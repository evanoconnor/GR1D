! -*-f90-*-
subroutine input_parser

  use GR1D_module
  use ye_of_rho
  implicit none

  integer :: tempint
  
  call get_string_parameter('jobname',jobname)
  call get_string_parameter('outdir',outdir)
  call get_logical_parameter('GR',GR)
  if (.not.GR) then
     call get_logical_parameter('do_effpot',do_effectivepotential)
  endif
  
  call get_string_parameter('initial_data',initial_data)
  if (initial_data.eq.'Shocktube') then
     call get_integer_parameter('shocktube_problem',shocktube_problem)
  endif
  if (initial_data.eq.'Collapse') then
     call get_string_parameter('profile_name',profile_name)
     call get_logical_parameter('WHW02profile',WHW02)
     call get_integer_parameter('profile_type',profile_type)
     if (profile_type.ne.1) stop "Unknown profile"
  endif
  call get_integer_parameter('geometry',geometry)
  if(geometry.eq.1) then
     if (initial_data.ne.'Shocktube') stop "geometry not matching initial data"
  else
     if (initial_data.eq.'Shocktube') stop "geometry not matching initial data"
  endif

  call get_string_parameter('gridtype',gridtype)
  if(gridtype.eq.'log') then
     call get_double_parameter('grid_custom_dx1',grid_custom_dx1)
  endif
  if(gridtype.eq.'custom') then
     call get_double_parameter('grid_custom_rad1',grid_custom_rad1)
     call get_double_parameter('grid_custom_dx1',grid_custom_dx1)
  endif
  if(gridtype.eq.'custom2') then
     call get_double_parameter('grid_custom_rad1',grid_custom_rad1)
     call get_double_parameter('grid_custom_dx1',grid_custom_dx1)     
     call get_double_parameter('grid_custom_inner',grid_custom_inner)
     call get_integer_parameter('grid_custom_number',grid_custom_number)
  endif
  if (initial_data.eq.'Collapse') then
     call get_logical_parameter('rmax_from_profile',do_profile_rmax)
     if (do_profile_rmax) then
        call get_double_parameter('rho_cut',rho_cut)
     else
        call get_double_parameter('grid_rmax',grid_rmax)
     endif
  endif

  call get_integer_parameter('radial_zones',radial_zones)
  call get_integer_parameter('ghosts',ghosts1)

  call get_logical_parameter('do_hydro',do_hydro)
  call get_logical_parameter('do_M1',do_M1)
  call get_logical_parameter('fake_neutrinos',fake_neutrinos)

  call get_integer_parameter('ntmax',ntmax)
  call get_double_parameter('tend',tend)
  if (initial_data.eq.'Collapse') then
     call get_logical_parameter('vs_mass',vs_mass)
     call get_logical_parameter('small_output',small_output)
     call get_logical_parameter('dynamic_output_control',&
          dynamic_output_control)
     if (.not.dynamic_output_control) then
        stop "Collapse should have dynamic control"
     endif
  endif
  call get_double_parameter('dtout',dtout0)
  dtout = dtout0
  call get_double_parameter('dtout_scalar',dtout_scalar0)
  dtout_scalar = dtout_scalar0
  call get_double_parameter('dtout_restart',dtout_restart0)
  dtout_restart = dtout_restart0
  call get_integer_parameter('ntout',ntout0)
  ntout = ntout0
  call get_integer_parameter('ntout_scalar',ntout_scalar0)
  ntout_scalar = ntout_scalar0
  call get_integer_parameter('ntout_restart',ntout_restart0)
  ntout_restart = ntout_restart0
  call get_integer_parameter('ntinfo',ntinfo0)
  ntinfo = ntinfo0
  call get_double_parameter('cffac',cffac)
  call get_integer_parameter('iorder_hydro',iorder_hydro)
  call get_logical_parameter('gravity_active',gravity_active)

  call get_string_parameter('reconstruction_method',&
       reconstruction_method)
  if(reconstruction_method.eq.'ppm') then
     call get_integer_parameter('ppm_origin_TVD',ppm_origin_TVD)
  endif
  call get_string_parameter('tvd_limiter',tvd_limiter)
  call get_string_parameter('flux_type',flux_type)
  if (flux_type.ne.'HLLE') stop "Flux type not recognized"
  call get_integer_parameter('eoskey',eoskey)
  if(eoskey.eq.3) then
     call get_string_parameter('eos_table_name',eos_table_name)
  endif
  if(eoskey.eq.1) then
     call get_double_parameter('hybridgamma_th',hybridgamma_th)
     call get_double_parameter('hybridgamma1',hybridgamma1)
     call get_double_parameter('hybridgamma2',hybridgamma2)
  endif
  if (initial_data.eq.'Collapse') then
     call get_logical_parameter('fake_neutrinos',fake_neutrinos)
     if(fake_neutrinos) then
        if (eoskey.ne.3) stop "Neutrinos only work with nuc_eos"
        call get_logical_parameter('do_leak_ros',do_leak_ros)
        call get_logical_parameter('do_heating',do_heating)
        call get_double_parameter('heat_fac',heat_fac)
        call get_logical_parameter('do_NNBrem',do_NNBrem)
        call get_logical_parameter('neutrino_pressure',do_nupress)
        call get_logical_parameter('ye_of_rho', do_ye_of_rho)
        if (do_ye_of_rho) then
           call get_logical_parameter('do_yeofrhofit', do_yeofrhofit)
           if (do_yeofrhofit) then
              call get_double_parameter('yeofrho_rho1',yeofrho_logrho1)
              yeofrho_logrho1 = log10(yeofrho_logrho1)
              call get_double_parameter('yeofrho_rho2',yeofrho_logrho2)
              yeofrho_logrho2 = log10(yeofrho_logrho2)
              call get_double_parameter('yeofrho_ye1',yeofrho_ye1)
              call get_double_parameter('yeofrho_ye2',yeofrho_ye2)
              call get_double_parameter('yeofrho_yec',yeofrho_yec)
              call get_logical_parameter('do_highcorrection',do_highcorrection)
              if (do_highcorrection) then
                 call get_double_parameter('yeofrho_rhoH',yeofrho_logrhoH)
                 yeofrho_logrhoH = log10(yeofrho_logrhoH)
                 call get_double_parameter('yeofrho_yeH',yeofrho_yeH)
              end if
           else
              call get_string_parameter('ye_profile_name',yeprofile_name)
           endif
        else
           if (do_nupress) stop "Neutrino pressure without ye_of_rho is not good"
        endif
     endif
  endif
  if(do_M1) then
     call get_integer_parameter('v_order',v_order)
     if (fake_neutrinos) stop "M1 neutrino's are not fake, how dare you (fake_neutrino needs to be 0 for M1)" 
     call get_double_parameter('evolution_radii',M1_maxradii)
     call get_double_parameter('extraction_radii',M1_extractradii)
     call get_integer_parameter('number_species',number_species)
     call get_integer_parameter('number_groups',number_groups)
     call get_integer_parameter('number_eas',number_eas)
     call get_logical_parameter('include_epannihil_kernels',include_epannihil_kernels)
     call get_logical_parameter('include_nes_kernels',include_nes_kernels)
     call get_integer_parameter('nes_evolution_type',tempint)
     if (tempint.eq.0) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .false.
     else if (tempint.eq.1) then
        include_Ielectron_exp = .true.
        include_Ielectron_imp = .false.
        if (include_nes_kernels.eqv..false.) then
           stop "please enable nes kernel read"
        endif
     else if (tempint.eq.2) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .true.
        if (include_nes_kernels.eqv..false.) then
           stop "please enable nes kernel read"
        endif
     else
        stop "unknown value for nes_evolution_exp"
     endif
     call get_integer_parameter('energycoupling_evolution_type',tempint)
     if (tempint.eq.0) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .false.
     else if (tempint.eq.1) then
        include_energycoupling_exp = .true.
        include_energycoupling_imp = .false.
     else if (tempint.eq.2) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .true.
     else
        stop "unknown value for energycoupling_evolution_exp"
     endif

     call get_string_parameter('M1closure',M1closure)     
     call get_string_parameter('opacity_table',opacity_table)
     call get_integer_parameter('testcase',M1_testcase_number)
     if (M1_testcase_number.gt.0.and.initial_data.ne."M1test") stop "set initial_data to `M1test'"
     if (M1_testcase_number.gt.9) then
        stop "add this in to input_parser"
     endif

     call get_logical_parameter('M1_control',M1_control_flag)
     if (M1_control_flag) then
        !read in phase settings
        call get_double_parameter('M1_phase1phase2_density',M1_phase1phase2_density)
        call get_double_parameter('M1_phase2phase3_pbtime',M1_phase2phase3_pbtime)
        !phase 1
        call get_string_parameter('M1_phase1_reconstruction',M1_phase1_reconstruction)
        call get_double_parameter('M1_phase1_cffac',M1_phase1_cffac)
        call get_integer_parameter('M1_phase1_ns',M1_phase1_ns)
        call get_integer_parameter('M1_phase1_ies_way',M1_phase1_ies_way)
        call get_integer_parameter('M1_phase1_encpl_way',M1_phase1_encpl_way)

        if ((M1_phase1_ies_way.eq.0).and.(include_Ielectron_exp.or.include_Ielectron_imp)) then
           write(*,*) "Initial M1 inelastic electron scattering method not the same as specified in main parameters"
           stop
        endif
        if (M1_phase1_ies_way.eq.1.and.(include_Ielectron_imp.or..not.include_Ielectron_exp)) then
           write(*,*) "Initial M1 inelastic electron scattering method not the same as specified in main parameters"
           stop
        endif
        if (M1_phase1_ies_way.eq.2.and.(include_Ielectron_exp.or..not.include_Ielectron_exp)) then
           write(*,*) "Initial M1 inelastic electron scattering method not the same as specified in main parameters"
           stop
        endif
        if (M1_phase1_encpl_way.eq.0.and.(include_energycoupling_exp.or.include_energycoupling_imp)) then
           write(*,*) "Initial M1 energy coupling method not the same as specified in main parameters"
           stop
        endif
        if (M1_phase1_encpl_way.eq.1.and.(include_energycoupling_imp.or..not.include_energycoupling_exp)) then
           write(*,*) "Initial M1 energy coupling method not the same as specified in main parameters"
           stop
        endif
        if (M1_phase1_encpl_way.eq.2.and.(include_energycoupling_exp.or..not.include_energycoupling_imp)) then
           write(*,*) "Initial M1 energy coupling method not the same as specified in main parameters"
           stop
        endif
        !phase 2
        call get_string_parameter('M1_phase2_reconstruction',M1_phase2_reconstruction)
        call get_double_parameter('M1_phase2_cffac',M1_phase2_cffac)
        call get_integer_parameter('M1_phase2_ns',M1_phase2_ns)
        call get_integer_parameter('M1_phase2_ies_way',M1_phase2_ies_way)
        call get_integer_parameter('M1_phase2_encpl_way',M1_phase2_encpl_way)
        !phase 3
        call get_string_parameter('M1_phase3_reconstruction',M1_phase3_reconstruction)
        call get_double_parameter('M1_phase3_cffac',M1_phase3_cffac)
        call get_integer_parameter('M1_phase3_ns',M1_phase3_ns)
        call get_integer_parameter('M1_phase3_ies_way',M1_phase3_ies_way)
        call get_integer_parameter('M1_phase3_encpl_way',M1_phase3_encpl_way)
     endif
  endif

  call get_double_parameter('atmo_rho_abs_min',atmo_rho_abs_min)
  call get_double_parameter('atmo_rho_rel_min',atmo_rho_rel_min)
  call get_double_parameter('atmo_fac',atmo_fac)
  call get_logical_parameter('do_restart',do_restart)
 
  if (initial_data.eq.'Collapse') then
     call get_logical_parameter('do_rotation',do_rotation)
     if (do_rotation) then
        if (geometry.ne.2) stop "Rotation in 1D?"
        call get_logical_parameter('set_omega',set_omega)
        call get_double_parameter('omega_c',omega_c)
        call get_double_parameter('omega_A',omega_A)
        !convert to code units, input should be in cgs
        omega_A = omega_A*length_gf
        omega_c = omega_c/time_gf
     endif
  endif
     
  if(do_restart) then
     call get_string_parameter('restart_file_name',restart_file_name)
     call get_logical_parameter('force_restart_dump',force_restart_output)
  endif

contains

 subroutine get_string_parameter(parname,par)

   implicit none
   character*(*) parname
   character(*) par
   character*(200) line_string
   integer i,j,l,ll
   integer isokay
   character*(200) temp_string

   open(unit=27,file='parameters',status='unknown')

10 continue
   read(27,'(a)',end=19) line_string

   ! separator is an equal sign '=', # is comment
   i = index(line_string,'=')
   j = index(line_string,'#')

   if (i.eq.0.or.j.eq.1) goto 10
   !   if the whole line is a comment or there is no
   !   equal sign, then go on to the next line    

   if(j.gt.0.and.j.lt.i) goto 10
   !   if there is an equal sign, but it is in a comment
   !   then go on to the next line


   ! is this the right parameter? If not, cycle
   temp_string=trim(adjustl(line_string(1:i-1)))
   l=len(parname)
   if(parname.ne.temp_string(1:l)) goto 10

   !  If there is a comment in the line, exclude it!
   l = len(line_string)
   if (j.gt.0) l = j - 1

   par = line_string(i+1:l)
   ! now remove potential crap!
   do ll=1,len(par)
      if(par(ll:ll).eq.'\t') par(ll:ll) = ' '
      if(par(ll:ll).eq.'"') par(ll:ll) = ' '
      if(par(ll:ll).eq."'") par(ll:ll) = ' '
   enddo
   do ll=1,len(par) 
      if(isokay(par(ll:ll)).ne.1) then
         par(ll:ll) = ' '
      endif
   enddo

   ! adjust left...
   ll = len(par)
   do while (par(1:1).eq.' '.or.par(1:1).eq.'\t') 
      par(1:ll-1) = par(2:ll)
   end do
   ! get rid of trailing blanks
   j = index(par," ")
   par = par(1:j-1)

   ! now look for " or ' and remove them
   j=index(par,'"')
   if(j.ne.0) stop "No quotes in my strings, please!"

   j=index(par,"'")
   if(j.ne.0) stop "No quotes in my strings, please!"

   close(27)
   return


19 continue
   write(6,*) "Fatal problem in input parser:"
   write(6,*) "Parameter ",parname
   write(6,*) "could not be read!"
   write(6,*) 
   call flush(6)
   stop

 end subroutine get_string_parameter

 subroutine get_double_parameter(parname,par)

   implicit none

   character(*) parname
   character*256 line_string
   real*8 par

   call get_string_parameter(parname,line_string)

   if(index(line_string,'.').eq.0) then
      write(6,*) "Uh. Bad double parameter ",trim(parname)
      write(6,*) "Please check input file!"
      call flush(6)
      stop
   endif

   read(line_string,*) par

 end subroutine get_double_parameter

 subroutine get_integer_parameter(parname,par)

   implicit none

   character(*) parname
   character*256 line_string
   integer par

   call get_string_parameter(parname,line_string)

   read(line_string,*) par

 end subroutine get_integer_parameter

 subroutine get_logical_parameter(parname,par)

   implicit none

   character*(*) parname
   character*(50) value_string
   integer temp_par
   logical par


   call get_string_parameter(parname,value_string)

   value_string = trim(adjustl(value_string))
   read(value_string,*) temp_par

   if(temp_par.ne.0) then
      par = .true.
   else
      par = .false.
   endif

 end subroutine get_logical_parameter

 

end subroutine input_parser

 function isokay(checkme)
   implicit none
   integer isokay
   character(1) checkme
   character(len=67) okaychars
   integer i
   isokay = 0
   okaychars = "1234567890.-+/_abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
   do i=1,len(okaychars)
      if(checkme(1:1) .eq. okaychars(i:i)) then
         isokay = 1
         exit
      endif
   enddo

 end function isokay
