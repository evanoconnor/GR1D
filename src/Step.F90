! -*-f90-*-
subroutine Step(dts)

  use GR1D_module
  use ye_of_rho
#if HAVE_LEAK_ROS
  use leakage_rosswog
#endif
  implicit none
  
  real*8 dts !time step that is passed in

  real*8 beta_rk,alpha_rk,cooling_rk_1,cooling_rk_2

  integer rkindex,i,m,gi,itc,j,k
  integer(kind=4) :: eosflag,keyerr,keytemp
  real*8 eosdummy(15)
  real*8 tempeps1(n1), tempeps2(n1)
  real*8 epsin0
  
  logical nan,inf

  !M1 stuff
  real*8 implicit_factor

  !set up conserved variables
  call prim2con

  !GR, do not need sqrt_gamma
  if (GR) then
     q_hat(:,:) = q(:,:)
     q_hat_old(:,:) = q(:,:)
  else
     qold(:,:) = q(:,:)
     do m=1,n_cons
        q_hat(:,m) = sqrt_gamma(:) * q(:,m)
     enddo
     q_hat_old(:,:) = q_hat(:,:)
  endif

  if(iorder_hydro.eq.2) then
     alpha_rk = 2.0d0
     beta_rk = 1.0d0
     cooling_rk_1 = 0.5d0
     cooling_rk_2 = 0.5d0
  else if(iorder_hydro.eq.3) then
     alpha_rk = 4.0d0
     beta_rk = 3.0d0
     cooling_rk_1 = 0.25d0
     cooling_rk_2 = 1.0d0/12.0d0
  else
     alpha_rk = 1.0d0
     beta_rk = 1.0d0
     cooling_rk_1 = 1.0d0
     cooling_rk_2 = 0.0d0 !not used
  endif

  denergyloss(:) = 0.0d0

  if(.not.do_hydro) goto 123

  do rkindex=1,iorder_hydro

     rkstep = rkindex

     gravsource(:,:) = 0.0d0
     flux_diff(:,:) = 0.0d0
     coolingsource(:,:) = 0.0d0
     presssource(:,:) = 0.0d0

     call reconstruct 
     call boundaries(0,0)

     if(flux_type .eq. "HLLE") then
        call flux_differences_hlle
     else
        stop "Sorry. Don't have the Riemann solver you want..."
     endif

     if(gravity_active.and.(geometry.eq.2)) then
       call gravity
     endif

     if(hydro_formulation.eq."conservative_p_source".and.(geometry.eq.2)) then
        call press_sources 
     endif
     
     if (GR.and.gravity_active) then
        if (do_nupress) then
           call nu_press_sources
        endif
     endif

#if HAVE_LEAK_ROS
     if (do_leak_ros.and.bounce) then
        call leak_rosswog
     endif
#endif

     if (rkindex .eq. 1 ) then
        do i=ghosts1,n1-1
           
           ! rho,D
           q_hat(i,1) = q_hat_old(i,1) + dts * ( - flux_diff(i,1) )
           
           ! rho*v, S
           q_hat(i,2) = q_hat_old(i,2) + dts * ( - flux_diff(i,2) &
                + gravsource(i,2) + presssource(i,2) + &
                coolingsource(i,2))
           
           ! energy, tau
           q_hat(i,3) = q_hat_old(i,3) + dts * ( - flux_diff(i,3) &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3))
           denergyloss(i) = cooling_rk_1*coolingsource(i,3)*dts
           
           ! ye
           q_hat(i,4) = q_hat_old(i,4) + dts * ( - flux_diff(i,4) &
                + coolingsource(i,4))
           
        enddo
        
        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = q_hat_old(i,5) + dts * ( - flux_diff(i,5) &
                   + presssource(i,5) )
           enddo
        endif
        
     elseif (rkindex .eq. 2 ) then
        do i=ghosts1,n1-1
           q_hat(i,1) = ( beta_rk * q_hat_old(i,1) + q_hat(i,1)  &
                + dts * ( - flux_diff(i,1) ) ) / alpha_rk
           
           q_hat(i,2) = ( beta_rk * q_hat_old(i,2) + q_hat(i,2)  &
                + dts * ( - flux_diff(i,2)    &
                + gravsource(i,2) & 
                + presssource(i,2) + coolingsource(i,2)) ) / alpha_rk
           
           q_hat(i,3) = ( beta_rk * q_hat_old(i,3) + q_hat(i,3)  &
                + dts * ( - flux_diff(i,3)    &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3) ) ) / alpha_rk
           denergyloss(i) = denergyloss(i) + cooling_rk_2 * &
                coolingsource(i,3)*dts
           
           q_hat(i,4) = ( beta_rk * q_hat_old(i,4)           &
                + q_hat(i,4)                            &
                + dts * ( - flux_diff(i,4)    &
                + coolingsource(i,4) ) ) / alpha_rk
           
        enddo

        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = ( beta_rk * q_hat_old(i,5)           &
                   + q_hat(i,5)                            &
                   + dts * ( - flux_diff(i,5)    &
                   + presssource(i,5) ) ) / alpha_rk
           enddo
        endif
        
     elseif (rkindex .eq. 3) then
        do i=ghosts1,n1-1
           q_hat(i,1) = ( q_hat_old(i,1) + 2.0d0*q_hat(i,1)  &
                + 2.0d0*dts * ( - flux_diff(i,1)) ) / 3.0d0
           
           q_hat(i,2) = ( q_hat_old(i,2) + 2.0d0*q_hat(i,2)  &
                + 2.0d0*dts * ( - flux_diff(i,2)    &
                + gravsource(i,2) & 
                + presssource(i,2) + coolingsource(i,2)) ) / 3.0d0
           
           q_hat(i,3) = ( q_hat_old(i,3) + 2.0d0*q_hat(i,3)  &
                + 2.0d0*dts * ( - flux_diff(i,3)    &
                + gravsource(i,3) + presssource(i,3) + &
                coolingsource(i,3) ) ) / 3.0d0
           denergyloss(i) = denergyloss(i) + 2.0d0/3.0d0 * &
                coolingsource(i,3)*dts
           
           q_hat(i,4) = ( q_hat_old(i,4) + 2.0d0*q_hat(i,4)  &
                + 2.0d0*dts * ( - flux_diff(i,4)    &
                + coolingsource(i,4) ) ) / 3.0d0
        enddo
        
        if(do_rotation) then
           do i=ghosts1,n1-1
              q_hat(i,5) = ( q_hat_old(i,5) + 2.0d0*q_hat(i,5)  &
                   + 2.0d0*dts * ( - flux_diff(i,5)    &
                   + presssource(i,5) ) ) / 3.0d0
           enddo
        endif
     else 
        stop 'Only iorder_hydro = 1, 2, and 3 implemented!'
     endif
     
     do m=1,n_cons
        if (GR) then
           q(:,m) = q_hat(:,m)
	else 
           q(:,m) = q_hat(:,m) / sqrt_gamma(:)
        endif
     enddo

     if (GR.and.gravity_active) then
        !find mgrav & X
        call con2GR
     endif
     
     !reconstruct primatives
     call con2prim

     ! eos update, eps fixed, find temp,entropy,cs2 etc.
     do i=ghosts1+1,n1-ghosts1
        keyerr = 0
        keytemp = 0
        tempeps1(i) = eps(i)
        call eos_full(i,rho(i),temp(i),ye(i),eps(i),press(i),pressth(i), & 
             entropy(i), &
             cs2(i), & 
             eosdummy(2),&
             eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
             eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
             eosdummy(11),eosdummy(12),eosdummy(13),nuchem(i), &
             keytemp,keyerr,eoskey,eos_rf_prec)
        tempeps2(i) = eps(i)
        if(keyerr.ne.0) then
           ! -> Issues with the EOS, this can happen around bounce
           !    and is due to very large temperature gradients
           !    in the bouncing inner core in adiabatic collapse
           !    for which the EOS was not really designed. The
           !    problems seen here should not show up for leakage/ye_of_rho
           !    runs.
           write(6,*) "############################################"
           write(6,*) "EOS PROBLEM in Step:"
           write(6,*) "timestep number: ",nt
           write(6,"(i4,1P10E15.6)") i,x1(i),rho(i)/rho_gf,temp(i),eps(i)/eps_gf,ye(i)
           write(6,*) "keyerr: ",keyerr
           call flush(6)
           if(.not.fake_neutrinos.and..not.do_leak_ros) then
              ! let's pump in a tiny bit of energy < 0.01*epsin_orig
              itc = 0
              epsin0 = eps(i)
              do while(keyerr.ne.0.and.itc.lt.10) 
                 itc = itc + 1
                 eps(i) = eps(i) + epsin0 * 1.0001d0
                 call eos_full(i,rho(i),temp(i),ye(i),eps(i),press(i),pressth(i), & 
                      entropy(i), &
                      cs2(i), & 
                      eosdummy(2),&
                      eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                      eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                      eosdummy(11),eosdummy(12),eosdummy(13),nuchem(i), &
                      keytemp,keyerr,eoskey,eos_rf_prec)
              enddo
              write(6,*) itc,keyerr
              write(6,*) "############################################"
              if(keyerr.ne.0) then
                 stop "problem in reconstruct: Step, could not be fixed."
              endif
           else
              stop "problem in reconstruct: Step"
           endif
        endif
     enddo
     
     !GR gravity updates that rely on primitive variables
     if (GR.and.gravity_active) then
        call GR_alp
        call GR_boundaries
     elseif (GR) then
        if (geometry.eq.2) then
           !GR Sedov, reflective boundaries on inside
           call GR_boundaries
        else 
           !planer, shocktube
           gi = 0
           do i=ghosts1,1,-1
              gi=gi+1
              v(i) = 0.0d0
              vp(i) = 0.0d0
              vm(i) = 0.0d0   
              W(i) = 1.0d0
           enddo
           do i=n1-ghosts1,n1
              gi=n1-ghosts1-1
              v(i) = v(gi)
              W(i) = W(gi)
           enddo
        endif
     endif

     call mass_interior

     if (do_nupress) then
        call neutrino_pressure
     endif

    call boundaries(0,0)
    
 enddo

123 continue

 !do operator split here
 !M1
 if (do_M1) then

    !we need to find the new plus/minus states, GR (alp,X) boundaries are done
    call reconstruct
    call boundaries(0,0)
    !If Newtonian, need to set v to v1 for velocity terms, shouldn't need to be
    if (.not.GR) then
       v = v1
       vp = v1p
       vm = v1m
    endif

    dyedt_hydro(:) = (ye(:) - ye_prev(:))/dts
    ye_prev(:) = ye(:)

    implicit_factor = 1.0d0

    q_M1_old = q_M1
    qold = q
       
    !reset source term
    M1_matter_source = 0.0d0

    !update interaction rates
    call M1_updateeas

    !reconstruct energy and flux in space and energy
    call M1_reconstruct

    !0. update closure variables
    call M1_closure
    
    !1. get explicit fluxes at time t^(n)
    call M1_explicitterms(dts,implicit_factor)

    if (M1_do_backwardfix.eq.1) then
       do j=1,number_groups
          do i=1,number_species
             do k=ghosts1+1,M1_imaxradii
                if ((eas(k,i,j,2)+eas(k,i,j,3))*(x1i(k+1)-x1i(k)).lt.0.01d0) then
                   B_M1(k,i,j,1:2) = -flux_M1(k,i,j,1:2)*0.5d0 !on the RHS now   
                else
                   B_M1(k,i,j,1:2) = -flux_M1(k,i,j,1:2) !on the RHS now       
                endif
             enddo
          enddo
       enddo
    else
       B_M1 = -flux_M1 !on the RHS now      
    endif
    C_M1 = -flux_M1_energy !explicit momentum flux, on the RHS now 
    D_M1 = flux_M1_scatter !explicit scattering, stays in the RHS

    !2. do implicit step for RHS source terms and calculate matter source terms
    call M1_implicitstep(dts,implicit_factor)

    !3. update matter conservatively
    if (do_hydro.or.(M1_testcase_number.eq.1.and.time.gt.0.0012d0)) then
       call M1_conservativeupdate(dts)
    endif

    !code for backward euler explicit flux fix
    if (M1_do_backwardfix.eq.1) then
       !reconstruct energy and flux in space and energy
       call M1_reconstruct
       !0. update closure variables
       call M1_closure
       !1. get explicit fluxes at time t^(n)
       call M1_explicitterms(dts,implicit_factor)

       do j=1,number_groups
          do i=1,number_species
             do k=ghosts1+1,M1_imaxradii
                if ((eas(k,i,j,2)+eas(k,i,j,3))*(x1i(k+1)-x1i(k)).lt.0.01d0) then
                   q_M1(k,i,j,1:2) = q_M1(k,i,j,1:2) - B_M1(k,i,j,1:2) - &
                        flux_M1(k,i,j,1:2)
                   if (q_M1(k,i,j,1).lt.0.0d0) then
                      write(*,*) k,i,j,q_M1(k,i,j,1)
                      stop "negative en after flux correct"
                   endif

                   if (abs(q_M1(k,i,j,2)/X(k)).gt.q_M1(k,i,j,1)) then
                      !fix it
                      q_M1(k,i,j,2) = X(k)*q_M1(k,i,j,2) / abs((1.0d0+1.0d-10)*q_M1(k,i,j,2)/q_M1(k,i,j,1))
                   endif
                endif
             enddo
          enddo
       enddo
    endif

    dyedt_neutrino(:) = (ye(:) - ye_prev(:))/dts

 endif

 !ye of rho prescription
 if(bounce) then
    if (do_ye_of_rho.and.(time.lt.t_bounce+0.005d0)) then
       call adjust_ye
       !things slightly change so redo all the variables
       call prim2con
       call boundaries(0,0)
       if(GR) then
          call con2GR
          call GR_alp
          call GR_boundaries
       endif
    endif
 else
    if (do_ye_of_rho) then
       call adjust_ye
       !things slightly change so redo all the variables
       call prim2con
       call boundaries(0,0)
       if(GR) then
          call con2GR
          call GR_alp
          call GR_boundaries
       endif
    endif
 endif

end subroutine Step

