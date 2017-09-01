!-*-f90-*-
!A primitive control system for handling different phases of
!core-collapse (collapse, bounce, postbounce)
subroutine M1_control 

  use GR1D_module
  implicit none

  integer :: phase

  !this routine handle controling the M1 parameters as a function of
  !time.  Crucial to maintain efficiency.

  if (initial_data.ne."Collapse") stop "add in controls for this initial_data"

  !determine the phase you are in:
  if (rho(ghosts1+1)/rho_gf.lt.M1_phase1phase2_density.and.M1_prev_phase.eq.1) then
     !in the collapse phase
     phase = 1
  else
     if (.not.bounce.or.(bounce.and.((time-t_bounce).lt.M1_phase2phase3_pbtime))) then
        !bounce phase
        phase = 2
     else
        !postbounce phase
        phase = 3
     endif
  endif

  if (phase.eq.1) then
     M1_reconstruction_method = M1_phase1_reconstruction
     cffac = M1_phase1_cffac
     number_species_to_evolve = M1_phase1_ns
     if (M1_phase1_ies_way.eq.0) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .false.
     else if (M1_phase1_ies_way.eq.1) then
        include_Ielectron_exp = .true.
        include_Ielectron_imp = .false.
     else if (M1_phase1_ies_way.eq.2) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .true.
     else
        stop "ies_way not supported"
     endif
     if (M1_phase1_encpl_way.eq.0) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .false.
     else if (M1_phase1_encpl_way.eq.1) then
        include_energycoupling_exp = .true.
        include_energycoupling_imp = .false.
     else if (M1_phase1_encpl_way.eq.2) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .true.
     else
        stop "encpl_way not supported"
     endif
  else if (phase.eq.2) then
     M1_reconstruction_method = M1_phase2_reconstruction
     cffac = M1_phase2_cffac 
     number_species_to_evolve = M1_phase2_ns
     if (M1_phase2_ies_way.eq.0) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .false.
     else if (M1_phase2_ies_way.eq.1) then
        include_Ielectron_exp = .true.
        include_Ielectron_imp = .false.
     else if (M1_phase2_ies_way.eq.2) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .true.
     else
        stop "ies_way not supported"
     endif
     if (M1_phase2_encpl_way.eq.0) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .false.
     else if (M1_phase2_encpl_way.eq.1) then
        include_energycoupling_exp = .true.
        include_energycoupling_imp = .false.
     else if (M1_phase2_encpl_way.eq.2) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .true.
     else
        stop "encpl_way not supported"
     endif
  else
     M1_reconstruction_method = M1_phase3_reconstruction  
     cffac = M1_phase3_cffac
     number_species_to_evolve = M1_phase3_ns
     if (M1_phase3_ies_way.eq.0) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .false.
     else if (M1_phase3_ies_way.eq.1) then
        include_Ielectron_exp = .true.
        include_Ielectron_imp = .false.
     else if (M1_phase3_ies_way.eq.2) then
        include_Ielectron_exp = .false.
        include_Ielectron_imp = .true.
     else
        stop "ies_way not supported"
     endif
     if (M1_phase3_encpl_way.eq.0) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .false.
     else if (M1_phase3_encpl_way.eq.1) then
        include_energycoupling_exp = .true.
        include_energycoupling_imp = .false.
     else if (M1_phase3_encpl_way.eq.2) then
        include_energycoupling_exp = .false.
        include_energycoupling_imp = .true.
     else
        stop "encpl_way not supported"
     endif
  endif

  M1_prev_phase = phase

end subroutine M1_control
