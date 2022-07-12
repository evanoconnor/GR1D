! -*-f90-*-
subroutine SetTimeStep
  
  use GR1D_module
  implicit none

  real*8 sound,dtnew
  integer keytemp,eosflag,keyerr
  integer i
  logical nan,inf
  
  if (GR) then
     call findnaninf(v,n1,nan,inf)
  else
     call findnaninf(v1,n1,nan,inf)
  endif
  
  if(nan) stop "NaNs found :-/"
    
  if (initial_data.eq."M1test".and.M1_testcase_number.gt.9) then
     stop "make sure dt works with grid"
     return
  endif

  ! sets time step
  ! Note: so far only cfl-stability conditions is implemented
  ! Other conditions may be added
  ! get the speed of sound from the eos
  do i=1,n1
     keytemp = 0 ! not coming in with temperature (that would reset the energy), needs to be zero for the hybrid/poly/ideal EOS
     eosflag = 6 ! we want cs2 to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in eos: cs2 at SetTimeStep"
     endif
  enddo
  dtp = dt

  dtnew = max(1.0d0*time_gf,dt_max*time_gf)
  do i=ghosts1+1, n1-ghosts1
     if (do_M1) then
        sound = 1.0d0
     else
        sound = sqrt(cs2(i))
     endif
     if (GR) then 
        if(.not.do_rotation) then
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(abs(v(i) + sound),abs(v(i) - sound)))
        else 
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(max(abs(v(i) + sound),abs(v(i) - sound)),&
                max(abs(vphi(i)+sound),abs(vphi(i)-sound))))
        endif
     else
        if(.not.do_rotation) then
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(abs(v1(i) + sound),abs(v1(i) - sound)))
        else
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(max(abs(v1(i) + sound),abs(v1(i) - sound)),&
                max(abs(vphi1(i)+sound),abs(vphi1(i)-sound))))
        endif
        
     endif
  enddo

  ! if not adaptive, set time step dt
  dt = min(dt_reduction_factor*cffac*dtnew,1.05d0*dtp)
  
  if (dtnew.eq.dt_max*time_gf) then
     dt_max = 1.01*dt_max
  endif

end subroutine SetTimeStep

subroutine handle_output
  
  use timers
  use GR1D_module
  implicit none
  
  if(dynamic_output_control) then
     call output_control
  endif
  if(mod(nt,ntinfo).eq. 0) then
     call outinfo
  endif
  
  if (nt.ge.ntmax) then
     write(*,*) "Done! :-) ntmax reached"
     call output_all(1)
     call output_all(2)
     call output_timers
     call restart_output_h5
     stop
  endif
  
  if (time.ge.tend) then
     write(*,*) "Done! :-) tend reached"
     call output_all(1)
     call output_all(2)
     call output_timers
     call restart_output_h5
     open(unit=666,file=trim(adjustl(outdir))//"/done",status="unknown")
     write(666,*) 1
     close(666)
     stop
  endif

  !!   Output/Checking
  if ( mod(nt,ntout) .eq. 0 .and. ntout.ne.-1) OutputFlag = .true.
  
  if ( mod(nt,ntout_restart) .eq. 0 .and. &
       ntout_restart.ne.-1) OutputFlagRestart = .true.
  
  if ( mod(nt,ntout_scalar).eq.0 .and. ntout_scalar.ne.-1) &
       OutputFlagScalar = .true.
  
  if ( time.ge.tdump) then
     tdump=tdump+dtout
     OutputFlag = .true.
  endif
  
  if ( time.ge.tdump_restart) then
     OutputFlagRestart = .true.
  endif
  
  if ( time.ge.tdump_scalar) then
     tdump_scalar=tdump_scalar+dtout_scalar
     OutputFlagScalar = .true.
  endif
  
  if(nt.eq.0) then
     OutputFlag = .true.
     OutputFlagScalar = .true.
  endif
  
  if (OutputFlag) then
     call output_all(1)
     call output_timers
     OutputFlag = .false.
  endif
  
  if (OutputFlagRestart) then
     call restart_output_h5
     tdump_restart=tdump_restart+dtout_restart
     OutputFlagRestart = .false.
  endif
  
  if (OutputFlagScalar) then
     call output_all(2)
     OutputFlagScalar = .false.
  endif
  
  if ((time+dt/time_gf).gt.tend) then
    dt = (tend-time)*time_gf
  endif
    
end subroutine handle_output

subroutine postStep_analysis
  
  use GR1D_module
  implicit none

  !for local EOS calls
  real*8 lrho,ltemp,lye,eosdummy(14)
  integer i,keytemp,keyerr
  
  if (initial_data.eq."Collapse".or.initial_data.eq."Collapse_inflow") then
     call get_shock_radius
     call mass_analysis    
     call get_binding_energy
     if (do_rotation) then
        call rotation_analysis
     endif
     if (M1_control_flag) then
        call M1_control
     endif
  endif
  
  nt=nt+1
  time = time+dt/time_gf
  time_c = time_c + dt/time_gf*alp(ghosts1+1)
  if (GR) then
     v_prev = v
  else
     v_prev = v1
  endif
  ye_prev = ye

  !fill mass fractions
#if HAVE_NUC_EOS
  keytemp = 1 !keep eps constant
  keyerr = 0
  do i=ghosts1+1,n1-ghosts1
     lrho = rho(i)/rho_gf
     ltemp = temp(i)
     lye = ye(i)
     call nuc_eos_full(lrho,ltemp,lye,eosdummy(1),eosdummy(2),eosdummy(3), &
          eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),massfrac_a(i),massfrac_h(i), &
          massfrac_n(i),massfrac_p(i),massfrac_abar(i),massfrac_zbar(i),eosdummy(8), &
          eosdummy(9),eosdummy(10),eosdummy(11),keytemp,keyerr,eos_rf_prec)
     if(keyerr.ne.0) then
        write(6,*) "############################################"
        write(6,*) "EOS PROBLEM in poststep analysis:"
        write(6,*) "timestep number: ",nt
        write(6,"(i4,1P10E15.6)") i,x1(i)/length_gf,lrho,ltemp,lye
        stop "Shouldn't fair here...."
     endif
  enddo
#endif



end subroutine postStep_analysis
