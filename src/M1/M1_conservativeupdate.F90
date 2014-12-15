!-*-f90-*- 
!this routine handles the updating of the conservative (or primative
!for Newtonian) variables after the radiation step is done.  it also
!computes the gain properties and monitors changes in ye for
!instability (near black hole fomration)

subroutine M1_conservativeupdate(dts)

  use GR1D_module

  implicit none

  real*8 :: dDye,dtau,dSr,dts
  integer k
  real*8 :: oneX,maxye
  integer gaintracker,maxyeloc
  logical nogain 
  integer keyerr,keytemp
  real*8 eosdummy(14)

  nogain = .true.
  total_net_heating = 0.0d0
  igain(1) = ghosts1+1
  gain_radius = 0.0d0
  maxye = 0.0d0
  maxyeloc = 0
  do k=ghosts1+1,M1_imaxradii

     if (GR) then
        oneX = X(k)
     else
        oneX = 1.0d0
     endif

     !convert this to number from energy, note since we evolve E, not
     !\sqrt{\gamma}E when determining the change in DYe we need to 
     dDye =  4.0d0*pi*dts*M1_matter_source(k,4)*oneX*(amu_cgs*mass_gf) ! \alp S^\alp u_\alp * X * 4*pi*dts
     dtau =  4.0d0*pi*dts*M1_matter_source(k,3) !-\alpha^2S^t * 4*pi*dts
     dSr = 4.0d0*pi*dts*M1_matter_source(k,2) ! -\alpha X S^r * 4*pi*dts

     !careful when determining gain radius, we require 10 zones of heating
     if (dtau.gt.0.0d0.and.nogain) then
        gaintracker = gaintracker+1
        if (gaintracker.eq.10) then
           igain(1) = k-9
           gain_radius = x1(igain(1))
           nogain = .false.
        endif
     else if (nogain) then 
        gaintracker = 0
        total_net_heating = 0.0d0
     else
        !have gain region, don't reset gaintracker
     endif

     if (gaintracker.gt.0.or..not.nogain) then
        if (dtau.gt.0.0d0) then
           total_net_heating = total_net_heating + &
                dtau*volume(k)/(energy_gf*dts/time_gf)
        endif
     endif

     total_energy_absorped = total_energy_absorped + &
         dtau*volume(k)/energy_gf/(dts/time_gf)

     q(k,2) = qold(k,2) + dSr
     q(k,3) = qold(k,3) + dtau
     q(k,4) = qold(k,4) + dDye

     if (abs(dDye/q(k,1)).gt.abs(maxye)) then
        maxye = dDye/q(k,1)
        maxyeloc = k
     endif
  enddo

  if(abs(maxye).gt.0.02d0) then
     !sometimes, near black hole formation, ye evolution can go
     !unstable right at nuclear matter density, likely something with
     !the opacities dramatically changing and maybe some heavies
     !kicking around, if this is the case, reduce time step a little,
     !throw a warning.  I think this is a sign that implicit coupling
     !of the matter and neutrinos is needed, but since it is such a
     !special case (black hole formation) this is fine for now

     dt_reduction_factor = dt_reduction_factor*0.9d0
     write(*,*) "Warning, ye seems unstable, reducing time step to compensate",dts,dt_reduction_factor
  endif

  !reset source terms for next step
  M1_matter_source = 0.0d0

  !update GR variables, resolve for primatives, update EOS
  !resolve for alpha, boundaries, update mass, 
  if (GR.and.gravity_active) then
     !find mgrav & X
     call con2GR
  endif

  !reconstruct primatives
  call con2prim

  ! eos update, eps fixed, find temp,entropy,cs2 etc.
  do k=ghosts1+1,n1-ghosts1
     keyerr = 0
     keytemp = 0
     call eos_full(k,rho(k),temp(k),ye(k),eps(k),press(k),pressth(k), & 
          entropy(k),cs2(k),eosdummy(2),&
          eosdummy(3),eosdummy(4),massfrac_a(k),massfrac_h(k), &
          massfrac_n(k),massfrac_p(k),massfrac_abar(k),massfrac_zbar(k), &
          elechem(k),eosdummy(12),eosdummy(13),eosdummy(14), &
          keytemp,keyerr,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        ! -> Issues with the EOS, this can happen around bounce
        !    and is due to very large temperature gradients
        !    in the bouncing inner core in adiabatic collapse
        !    for which the EOS was not really designed. The
        !    problems seen here should not show up for leakage/ye_of_rho
        !    runs.
        write(6,*) "############################################"
        write(6,*) "EOS PROBLEM in M1_conservativeupdate:"
        write(6,*) "timestep number: ",nt
        write(6,"(i4,1P10E15.6)") k,x1(k),rho(k)/rho_gf,temp(k),eps(k)/eps_gf,ye(k)
     endif
  enddo

  !GR gravity updates that rely on primitive variables
  if (GR.and.gravity_active) then
     call GR_alp
     call GR_boundaries
  endif
  call mass_interior
  call boundaries(0,0)   

end subroutine M1_conservativeupdate
