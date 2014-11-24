!-*-f90-*-
!this routine calculates all of the explicit terms (spatial
!flux,energy flux, and scattering terms) It returns the terms in the
!flux_M1 type variables flux_M1 (for spatial), flux_M1_energy (for
!energy), flux_M1_scatter (for scattering)
subroutine M1_explicitterms(dts,implicit_factor)
  
  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies,nulibtable_ewidths, &
       nulibtable_energies,nulibtable_etop,nulibtable_ebottom, &
       nulibtable_logenergies,nulibtable_logetop
  implicit none

  !inputs:
  real*8 :: dts !time step
  real*8 :: implicit_factor !time step

  !local, spatial
  real*8 :: M1en_space(n1),M1flux_space(n1)
  real*8 :: M1chi_space(n1),M1eddy_space(n1)
  real*8 :: M1en_space_plus(n1),M1en_space_minus(n1)
  real*8 :: M1flux_space_plus(n1),M1flux_space_minus(n1)
  real*8 :: M1chi_space_plus(n1),M1chi_space_minus(n1)
  real*8 :: M1eddy_space_plus(n1),M1eddy_space_minus(n1)
  real*8 :: l_min_thin(2),l_max_thin(2)
  real*8 :: l_min_thick(2),l_max_thick(2) 
  real*8 :: l_min(3),l_max(3)
  real*8 :: p,discrim,sqrtdiscrim
  real*8 :: M1flux_interface(n1,2),M1flux_diff(n1,2)
  real*8 :: rm,rp,dx

  !local, energy
  real*8 :: div_v(n1),dvdt(n1)
  real*8 :: h,dmdr,dmdt,dXdr,dXdt,Kdownrr,dWdt,dWdr,dWvuprdt,dWvuprdr
  real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6),heattermff_NL(6)
  real*8 :: velocity_coeffs(6,2)
  real*8 :: M1en_energy(number_groups),M1flux_energy(number_groups)
  real*8 :: M1eddy_energy(number_groups),M1pff_energy(number_groups)
  real*8 :: M1en_energy_fluid(number_groups),M1chi_energy(number_groups)
  real*8 :: M1qrrr_energy(number_groups),M1qffr_energy(number_groups)
  real*8 :: littlefactors_int(number_groups,6)
  real*8 :: littlefactors(number_groups,6)
  real*8 :: velocity(number_groups,2)
  real*8 :: velocity_top(number_groups,2)
  real*8 :: log_distro(number_groups)
  real*8 :: log_energy_interface_distroj(number_groups)
  real*8 :: energy_interface_distroj(number_groups)
  real*8 :: energy_interface_M1en(number_groups)
  real*8 :: M1flux_energy_interface(number_groups,2)
  real*8 :: em,ep,de
  real*8 :: M1flux_diff_energy(number_groups,2)
  real*8 :: M1_moment_to_distro_nobinwidth(number_groups) 
  real*8 :: M1_moment_to_distro_top_nobinwidth(number_groups) 
  real*8 :: nulibtable_etop_logged(number_groups)
  real*8 :: nulibtable_energies_logged(number_groups)

  real*8 :: a_asym,kappa_inter,ipeclet_mean
  real*8 :: Jkplus1,Jk,diffusive_flux,advected_energy

  real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in,ies_temp,species_factor
  real*8 :: ies_sourceterms(2*number_groups)
  real*8 :: local_M(2,2),local_J(number_groups),local_H(number_groups,2)
  real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
  real*8 :: local_Ltilde(number_groups,2,2),JoverE(number_groups)
  real*8 :: JoverF(number_groups),HoverE(number_groups,2),HoverF(number_groups,2)
  real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)
  real*8 :: invalp,invalp2,invX,invX2,X2,alp2,W2,v2,invr,invr2,onev,oneW,onealp,oneX,oneWm,oneWp
  real*8 :: local_u(2),local_littleh(2,2),local_littlehupup(2,2),local_uup(2)

  !mueller
  real*8 :: logdistro(number_groups)
  real*8 :: loginterface_distroj(number_groups)
  real*8 :: interface_distroj(number_groups)
  real*8 :: xi(number_groups)
  real*8 :: FL(number_groups,2),FR(number_groups,2)
  real*8 :: temp_term
  
  real*8 :: fluxtemp_1,fluxtemp_2
  real*8 :: limitingflux

  !counters
  integer i,j,k,j_prime,jj,ii

  !reset flux to zero
  flux_M1 = 0.0d0
  flux_M1_energy = 0.0d0
  flux_M1_scatter = 0.0d0

!#################################################################
!#################################################################
!########################Spatial##################################
!#################################################################
!#################################################################

  !first spatial. This is essentially a reimann solve + corrections
  !over all the radial zones. Variables are already interpolated

  !$OMP PARALLEL DO PRIVATE(i,M1en_space,M1flux_space,M1chi_space,M1eddy_space, &
  !$OMP M1en_space_plus,M1flux_space_plus,M1chi_space_plus,M1eddy_space_plus, &
  !$OMP M1en_space_minus,M1flux_space_minus,M1chi_space_minus,M1eddy_space_minus, &
  !$OMP k,oneWm,oneWp,kappa_inter,ipeclet_mean,a_asym,l_min_thin,l_max_thin, &
  !$OMP p,discrim,sqrtdiscrim,l_min_thick,l_max_thick,l_min,l_max,Jkplus1,Jk,diffusive_flux, &
  !$OMP advected_energy,M1flux_interface,limitingflux,rm,rp,dx,M1flux_diff)
  do j=1,number_groups
     do i=1,number_species_to_evolve

        M1en_space = q_M1(:,i,j,1)
        M1flux_space = q_M1(:,i,j,2)
        M1chi_space = q_M1_extra(:,i,j,4)
        M1eddy_space = q_M1(:,i,j,3)

        M1en_space_plus = q_M1p(:,i,j,1,1)
        M1en_space_minus = q_M1m(:,i,j,1,1)
        M1flux_space_plus = q_M1p(:,i,j,2,1)
        M1flux_space_minus = q_M1m(:,i,j,2,1)

        M1eddy_space_plus = q_M1p(:,i,j,3,1)
        M1eddy_space_minus = q_M1m(:,i,j,3,1)
        M1chi_space_plus = q_M1_extrap(:,i,j,1,1)
        M1chi_space_minus = q_M1_extram(:,i,j,1,1)

        !now find speeds at each interface
        do k=ghosts1,M1_imaxradii

           !only used in when GR.eq.0 and v_order.eq.-1
           oneWm = 1.0d0/sqrt(1.0d0-v1m(k+1)**2)
           oneWp = 1.0d0/sqrt(1.0d0-v1p(k)**2)

           !determine the regime for the flux calculation
           kappa_inter = sqrt((eas(k,i,j,2)+eas(k,i,j,3))*(eas(k+1,i,j,2)+eas(k+1,i,j,3)))
           if (v_order.eq.-1) then
              if (GR) then
                 ipeclet_mean = 1.0d0/(Wm(k+1)**3*(1.0d0+vm(k+1))*Xm(k+1)**2* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              else
                 ipeclet_mean = 1.0d0/(oneWm**3*(1.0d0+v1m(k+1))* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              endif
                 
           else if (v_order.eq.0) then
              if (GR) then
                 ipeclet_mean = 1.0d0/(Xm(k+1)**2* &
                      (kappa_inter)*(x1(k+1)-x1(k))) 
              else
                 ipeclet_mean = 1.0d0/(kappa_inter*(x1(k+1)-x1(k))) 
              endif
           else 
              stop "add in v order, peclet number"
           endif

           !luke's term
           a_asym = tanh(ipeclet_mean)
           if (a_asym.gt.1.0d0) then
              a_asym = 1.0d0
           endif

           !speeds
           !minus interface (k+1 zone)
           !thin limit:
           if (GR) then
              l_min_thin(1) = -Xm(k+1)*alpm(k+1)
              l_max_thin(1) = Xm(k+1)*alpm(k+1)
           else
              l_min_thin(1) = -1.0d0
              l_max_thin(1) = 1.0d0
           endif
           
           !thick limit:
           if (v_order.eq.-1) then
              if (GR) then
                 p = alpm(k+1)*vm(k+1)*Xm(k+1)
                 discrim = alpm(k+1)**2*Xm(k+1)**2*(2.0d0*Wm(k+1)**2+1.0d0)-2.0d0*Wm(k+1)**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(1) = min((2.0d0*Wm(k+1)**2*p - sqrtdiscrim)/(2.0d0*Wm(k+1)**2+1.0d0),p)
                 l_max_thick(1) = max((2.0d0*Wm(k+1)**2*p + sqrtdiscrim)/(2.0d0*Wm(k+1)**2+1.0d0),p)
              else
                 p = v1m(k+1)
                 discrim = (2.0d0*oneWm**2+1.0d0)-2.0d0*oneWm**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(1) = min((2.0d0*oneWm**2*p - sqrtdiscrim)/(2.0d0*oneWm**2+1.0d0),p)
                 l_max_thick(1) = max((2.0d0*oneWm**2*p + sqrtdiscrim)/(2.0d0*oneWm**2+1.0d0),p)
              endif


           else if (v_order.eq.0) then
              p = 0.0d0
              if (GR) then
                 discrim = alpm(k+1)**2*Xm(k+1)**2*3.0d0
              else
                 discrim = 3.0d0
              endif
              sqrtdiscrim = sqrt(discrim)
              l_min_thick(1) = -sqrtdiscrim/3.0d0
              l_max_thick(1) = sqrtdiscrim/3.0d0
           else
              stop "add me in"
           endif

           !actual speed
           l_min(1) = (3.0d0*M1chi_space_minus(k+1)-1.0d0)*0.5d0*l_min_thin(1) + &
                (1.0d0 - M1chi_space_minus(k+1))*1.5d0*l_min_thick(1)

           l_max(1) = (3.0d0*M1chi_space_minus(k+1)-1.0d0)*0.5d0*l_max_thin(1) + &
                (1.0d0 - M1chi_space_minus(k+1))*1.5d0*l_max_thick(1)

           !plus interface (k zone)
           !thin limit:
           if (GR) then
              l_min_thin(2) = -Xp(k)*alpp(k)
              l_max_thin(2) = Xp(k)*alpp(k)
           else
              l_min_thin(2) = -1.0d0
              l_max_thin(2) = 1.0d0
           endif
              
           
           !thick limit:
           if (v_order.eq.-1) then
              if (GR) then
                 p = alpp(k)*vp(k)*Xp(k)
                 discrim = alpp(k)**2*Xp(k)**2*(2.0d0*Wp(k)**2+1.0d0)-2.0d0*Wp(k)**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(2) = min((2.0d0*Wp(k)**2*p - sqrtdiscrim)/(2.0d0*Wp(k)**2+1.0d0),p)
                 l_max_thick(2) = max((2.0d0*Wp(k)**2*p + sqrtdiscrim)/(2.0d0*Wp(k)**2+1.0d0),p)
              else
                 p = v1p(k)
                 discrim = (2.0d0*Wp(k)**2+1.0d0)-2.0d0*Wp(k)**2*p*p
                 sqrtdiscrim = sqrt(discrim)
                 l_min_thick(2) = min((2.0d0*oneWp**2*p - sqrtdiscrim)/(2.0d0*oneWp**2+1.0d0),p)
                 l_max_thick(2) = max((2.0d0*oneWp**2*p + sqrtdiscrim)/(2.0d0*oneWp**2+1.0d0),p)
              endif

           else if (v_order.eq.0) then
              p = 0.0d0
              if (GR) then
                 discrim = alpp(k)**2*Xp(k)**2*3.0d0
              else
                 discrim = 3.0d0
              endif

              sqrtdiscrim = sqrt(discrim)
              l_min_thick(2) = -sqrtdiscrim/3.0d0
              l_max_thick(2) = sqrtdiscrim/3.0d0
           else
              stop "add me in"
           endif
           
           !actualy speed
           l_min(2) = (3.0d0*M1chi_space_plus(k)-1.0d0)*0.5d0*l_min_thin(2) + &
                (1.0d0 - M1chi_space_plus(k))*1.5d0*l_min_thick(2)
           l_max(2) = (3.0d0*M1chi_space_plus(k)-1.0d0)*0.5d0*l_max_thin(2) + &
                (1.0d0 - M1chi_space_plus(k))*1.5d0*l_max_thick(2)


           !check for NaNs
           if (l_min(1).ne.l_min(1)) then
              write(*,*) "NaNs in speeds 1"
              stop
           else if(l_min(2).ne.l_min(2)) then
              write(*,*) "NaNs in speeds 2"
              stop
           else if(l_max(1).ne.l_max(1)) then
              write(*,*) "NaNs in speeds 3"
              stop
           else if(l_max(2).ne.l_max(2)) then
              write(*,*) "NaNs in speeds 4"
              stop
           endif

           l_max(3) = max(l_max(1),l_max(2))
           l_min(3) = max(l_min(1),l_min(2))

           !Lukes idea/suggestion
           if (v_order.eq.-1) then
              if (GR) then
!                 Jkplus1 = W(k+1)**2*(M1en_space(k+1)*(1.0d0+v(k+1)**2* &
!                      M1eddy_space(k+1)/X(k+1)**2)-2.0d0*M1flux_space(k+1)*v(k+1)/X(k+1))
!                 Jk = W(k)**2*(M1en_space(k)*(1.0d0+v(k)**2* &
!                      M1eddy_space(k)/X(k)**2)-2.0d0*M1flux_space(k)*v(k)/X(k))

                 !shibata's approximation to above, probably should use above.
                 Jk = 3.0d0/(2.0d0*W(k)**2+1.0d0)*((2.0d0*W(k)**2-1.0d0)*M1en_space(k) - &
                      2.0d0*W(k)**2*v(k)/X(k)*M1flux_space(k))
                 Jkplus1 = 3.0d0/(2.0d0*W(k+1)**2+1.0d0)*((2.0d0*W(k+1)**2-1.0d0)*M1en_space(k+1) - &
                      2.0d0*W(k+1)**2*v(k+1)/X(k+1)*M1flux_space(k+1))

                 if (a_asym.eq.1.0d0) then
                    diffusive_flux = 0.0d0
                 else
                    diffusive_flux = -W(k)/(3.0d0*kappa_inter*X(k)**2)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif
                 
                 if (vp(k)*vm(k+1).gt.0.0d0) then
                    if (vp(k).lt.0.0d0) then
                       advected_energy = 4.0d0*Wm(k+1)**2*vm(k+1)*Xm(k+1)*Jkplus1*onethrd
                    else
                       advected_energy = 4.0d0*Wp(k)**2*vp(k)*Xp(k)*Jk*onethird
                    endif
                 else
                    advected_energy = 0.0d0
                 endif
                 
              else
                 Jkplus1 = (1.0d0/(1.0d0-v1(k+1)**2))*(M1en_space(k+1)*(1.0d0+v1(k+1)**2* &
                      M1eddy_space(k+1))-2.0d0*M1flux_space(k+1)*v1(k+1))
                 Jk = (1.0d0/(1.0d0-v1(k)**2))*(M1en_space(k)*(1.0d0+v1(k)**2* &
                      M1eddy_space(k))-2.0d0*M1flux_space(k)*v1(k))
              
                 if (a_asym.eq.1.0d0) then
                    diffusive_flux = 0.0d0
                 else
                    diffusive_flux = -1.0d0/(sqrt(1.0d0-v1(k)**2)*3.0d0*kappa_inter)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif

                 if (v1p(k)*v1m(k+1).gt.0.0d0) then
                    if (v1p(k).lt.0.0d0) then
                       advected_energy = 4.0d0*oneWm**2*v1m(k+1)*Jkplus1*onethird
                    else
                       advected_energy = 4.0d0*oneWp**2*v1p(k)*Jk*onethird
                    endif
                 else
                    advected_energy = 0.0d0
                 endif

              endif
                 
           else if (v_order.eq.0) then
              Jkplus1 = M1en_space(k+1)
              Jk = M1en_space(k)
              
              if (a_asym.eq.1.0d0) then
                 diffusive_flux = 0.0d0
              else
                 if (GR) then
                    diffusive_flux = -1.0d0/(3.0d0*kappa_inter*X(k)**2)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 else
                    diffusive_flux = -1.0d0/(3.0d0*kappa_inter)*(Jkplus1-Jk)/(x1(k+1)-x1(k))
                 endif
              endif
                 
              advected_energy = 0.0d0

           else
              stop "add v order"
           endif

           !flux at interface has two componants, asympotic part, free streaming part.
           M1flux_interface(k,1) = a_asym*( &
                ((l_max(3)*M1flux_space_plus(k)-&
                l_min(3)*M1flux_space_minus(k+1)) + l_max(3)*l_min(3)* &
                (M1en_space_minus(k+1)-M1en_space_plus(k)))/(l_max(3)-l_min(3)) &
                ) + &
                (1.0d0-a_asym)*( &
                diffusive_flux + advected_energy &
                )

           !shift interface flux to geometric mean is the difference
           !is too large.  This helps the deleptonziation burst
           !stablely evolve
           if (M1en_space_plus(k)/M1en_space_minus(k+1).gt.10.0d0) then
              limitingflux = sqrt(M1eddy_space_plus(k)*M1en_space_plus(k)* &
                   M1eddy_space_minus(k+1)*M1en_space_minus(k+1))
           else
              limitingflux = (M1eddy_space_plus(k)*M1en_space_plus(k)+ &
                   M1eddy_space_minus(k+1)*M1en_space_minus(k+1))/2.0d0
           endif

           M1flux_interface(k,2) = ( &
                a_asym*((l_max(3)*M1eddy_space_plus(k)*M1en_space_plus(k)- &
                l_min(3)*M1eddy_space_minus(k+1)*M1en_space_minus(k+1)) + &
                l_max(3)*l_min(3)*(M1flux_space_minus(k+1)- &
                M1flux_space_plus(k)))/(l_max(3)-l_min(3)) &
                ) + &
                (1.0d0-a_asym)*( &
                limitingflux &
                )
           
        enddo

        do k=ghosts1+1,M1_imaxradii
           rm = x1i(k)
           rp = x1i(k+1)
           dx = (rp-rm)

           if (GR) then
              M1flux_diff(k,1) = (alpp(k)/Xp(k)**2*x1i(k+1)**2*M1flux_interface(k,1)- &
                   alpm(k)/Xm(k)**2*x1i(k)**2*M1flux_interface(k-1,1))/(dx*x1(k)**2)
              M1flux_diff(k,2) = (alpp(k)/Xp(k)**2*x1i(k+1)**2*M1flux_interface(k,2)- &
                   alpm(k)/Xm(k)**2*x1i(k)**2*M1flux_interface(k-1,2))/(dx*x1(k)**2)
           else
              M1flux_diff(k,1) = (x1i(k+1)**2*M1flux_interface(k,1)- &
                   x1i(k)**2*M1flux_interface(k-1,1))/(dx*x1(k)**2)
              M1flux_diff(k,2) = (x1i(k+1)**2*M1flux_interface(k,2)- &
                   x1i(k)**2*M1flux_interface(k-1,2))/(dx*x1(k)**2)
           endif

           flux_M1(k,i,j,1) = dts*implicit_factor*M1flux_diff(k,1)
           flux_M1(k,i,j,2) = dts*implicit_factor*M1flux_diff(k,2)

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO! end do
  
  !let things get away from IC for a few time steps
  if (time.lt.0.000001d0) return

!################################################################
!################################################################
!######################Explicit Energy###########################
!################################################################
!################################################################

  if (include_energycoupling_exp) then
     !$OMP PARALLEL DO PRIVATE(W2,v2,oneW,onev,dmdr,dmdt,dXdr,dXdt,Kdownrr, &
     !$OMP dWdt,dWdr,dWvuprdt,dWvuprdr,Z,Yupr,Xuprr,Xupff,heatterm_NL,heattermff_NL, &
     !$OMP velocity_coeffs,i,M1en_energy,M1en_energy_fluid,M1flux_energy,M1eddy_energy, &
     !$OMP M1pff_energy,M1qrrr_energy,M1qffr_energy,j,littlefactors,velocity,M1flux_energy_interface, &
     !$OMP logdistro,loginterface_distroj,interface_distroj,xi,temp_term,FL,FR,ep,M1flux_diff_energy, &
     !$OMP em,de)
     do k=ghosts1+1,M1_imaxradii
        if (GR) then
           if (k.eq.ghosts1+1) then
              div_v(k) = v(k)/x1(k)
           else
              div_v(k) = (v(k+1)-v(k-1))/(x1(k+1)-x1(k-1))
           endif
        else
           if (k.eq.ghosts1+1) then
              div_v(k) = v1(k)/x1(k)
           else
              div_v(k) = (v1(k+1)-v1(k-1))/(x1(k+1)-x1(k-1))
           endif
        endif
        if (nt.eq.0) then
           dvdt(k) = 0.0d0
        else
           if (GR) then
              dvdt(k) = (v(k)-v_prev(k))/dts
           else
              dvdt(k) = (v1(k)-v_prev(k))/dts
           endif

        endif

        if (v_order.eq.-1) then
           if (GR) then
              W2 = W(k)**2
              v2 = v(k)**2
              oneW = W(k)
              onev = v(k)
           else
              W2 = 1.0d0/(1.0d0-v1(k)**2)
              v2 = v1(k)**2
              oneW = sqrt(W2)
              onev = v1(k)
           endif

        else if (v_order.eq.0) then
           W2 = 1.0d0
           oneW = 1.0d0
           v2 = 0.0d0
           onev = 0.0d0
           div_v(k) = 0.0d0
           dvdt(k) = 0.0d0
        endif

        if (GR) then
           !this W^2 has the kinetic energy contriubiton to the mass,
           !nothing to do with the transport, so we'll leave it in
           !regardless of v_order.
           dmdr = 4.0d0*pi*x1(k)**2*(rho(k)*h*W(k)**2-press(k)) 
           dmdt = -4.0d0*pi*x1(k)**2*alp(k)*rho(k)*h*W(k)**2*v(k)/X(k)
           dXdr = X(k)**3*(dmdr/x1(k)-mgrav(k)/x1(k)**2)
           dXdt = X(k)**3*dmdt/x1(k)
           Kdownrr = -X(k)*dXdt/alp(k)
           dWdt = W2*oneW*onev*dvdt(k)
           dWdr = W2*oneW*onev*div_v(k)
           dWvuprdt = onev/X(k)*dWdt - oneW*onev/X(k)**2*dXdt + oneW/X(k)*dvdt(k)
           dWvuprdr = onev/X(k)*dWdr - oneW*onev/X(k)**2*dXdr + oneW/X(k)*div_v(k)
        else
           dmdr = 0.0d0
           dmdt = 0.0d0
           dXdr = 0.0d0
           dXdt = 0.0d0
           Kdownrr = 0.0d0
           dWdt = W2*oneW*onev*dvdt(k)
           dWdr = W2*oneW*onev*div_v(k)
           dWvuprdt = onev*dWdt + oneW*dvdt(k)
           dWvuprdr = onev*dWdr + oneW*div_v(k)
        endif

        Z = 0.0d0
        Yupr = 0.0d0
        Xuprr = 0.0d0
        Xupff = 0.0d0
        heatterm_NL = 0.0d0
        heattermff_NL = 0.0d0

        !g_rr v^r v^r = v^2 -> v = Xv^r -> v = v_r/X -> v_r = vx -> v^r = v/x
        if (GR) then
           Z(1) = 1.0d0/oneW !E
           Z(2) = onev/(X(k)*oneW) !F
           Z(3) = v2/(X(k)**2*oneW) !prr
           Z(5) = v2*onev/X(k)**3 !qrrr
           Yupr(2) = 1.0d0/(X(k)**2*oneW) !F
           Yupr(3) = onev/(X(k)**3*oneW)
           Yupr(5) = v2/X(k)**4
           Xuprr(3) = 1.0d0/(X(k)**4*oneW)
           Xuprr(5) = onev/X(k)**5
           Xupff(4) = 1.0d0/(x1(k)**2*oneW)
           Xupff(6) = onev/(X(k)*x1(k)**2)
           heatterm_NL(5) = 1.0d0 
           heattermff_NL(6) = 1.0d0
        else
           Z(1) = 1.0d0/oneW !E
           Z(2) = onev/oneW !F
           Z(3) = v2/oneW !prr
           Z(5) = v2*onev !qrrr
           Yupr(2) = 1.0d0/oneW !F
           Yupr(3) = onev/oneW
           Yupr(5) = v2
           Xuprr(3) = 1.0d0/oneW
           Xuprr(5) = onev
           Xupff(4) = 1.0d0/(x1(k)**2*oneW)
           Xupff(6) = onev/x1(k)**2
           heatterm_NL(5) = 1.0d0 
           heattermff_NL(6) = 1.0d0
        endif

        !determine velocity coefficients
        !energy
        if (GR) then
           velocity_coeffs(:,1) = alp(k)*(oneW*((Z(:)*onev/X(k)-Yupr(:))*dphidr(k) - &
                Xuprr(:)*onev*dXdr - 2.0d0*Xupff(:)*onev/X(k)*x1(k) + &
                Xuprr(:)*Kdownrr) + Z(:)/alp(k)*dWdt + Yupr(:)*dWdr - &
                X(k)**2*Yupr(:)/alp(k)*dWvuprdt - X(k)**2*Xuprr(:)*dWvuprdr)

           !flux
           velocity_coeffs(:,2) = alp(k)*(oneW*((Yupr(:)*onev*X(k)-Xuprr(:)*X(k)**2)*dphidr(k) - &
                heatterm_NL(:)/X(k)**4*onev*dXdr - 2.0d0*heattermff_NL(:)/x1(k)**2*onev/X(k)*x1(k) + &
                heatterm_NL(:)/X(k)**4*Kdownrr) + Yupr(:)*X(k)**2/alp(k)*dWdt +  &
                Xuprr(:)*X(k)**2*dWdr - Xuprr(:)*X(k)**4/alp(k)*dWvuprdt - &
                heatterm_NL(:)/X(k)**2*dWvuprdr)
        else
           velocity_coeffs(:,1) = (oneW*(-2.0d0*Xupff(:)*onev*x1(k)) + &
                Z(:)/alp(k)*dWdt + Yupr(:)*dWdr - &
                Yupr(:)*dWvuprdt - Xuprr(:)*dWvuprdr)

           !flux
           velocity_coeffs(:,2) = (oneW*(-2.0d0*heattermff_NL(:)/x1(k)**2*onev*x1(k)) + &
                Yupr(:)*dWdt + Xuprr(:)*dWdr - Xuprr(:)*dWvuprdt - &
                heatterm_NL(:)*dWvuprdr)
        endif

        do i=1,number_species_to_evolve
           !get E and F for this species and spatial point.  If this
           !term was ignored, the \Delta \epsilon term in the
           !equations completely factors out, so we include it right
           !from the beginning (i.e. we change \eta when we read in
           !the NuLib table so \eta has units of erg/cm^3/s/srad
           !rather than erg/cm^3/s/srad/MeV).  But the energy
           !derivative prevents that from being true. Here, we replace
           !E and F with E/\Delta \epsilon and F/\Delta_\epsilon.  Out
           !side of the energy derivative term we multiple through by
           !\Delta_\epsilon to get it back in the right units (since
           !there is a / \Delta_\epsilon from the derivative, these
           !terms actually just cancel)
           M1en_energy = q_M1(k,i,:,1)/(nulibtable_etop(:)-nulibtable_ebottom(:))
           M1en_energy_fluid = q_M1_fluid(k,i,:,1)
           M1flux_energy = q_M1(k,i,:,2)/(nulibtable_etop(:)-nulibtable_ebottom(:))
           M1eddy_energy = q_M1(k,i,:,3)
           M1pff_energy = q_M1_extra(k,i,:,1)
           M1qrrr_energy = q_M1_extra(k,i,:,2)/(nulibtable_etop(:)-nulibtable_ebottom(:))
           M1qffr_energy = q_M1_extra(k,i,:,3)/(nulibtable_etop(:)-nulibtable_ebottom(:))

           do j=1,number_groups

              littlefactors(j,1) = M1en_energy(j)
              littlefactors(j,2) = M1flux_energy(j)
              littlefactors(j,3) = M1eddy_energy(j)*M1en_energy(j)
              littlefactors(j,4) = M1pff_energy(j)*M1en_energy(j)
              littlefactors(j,5) = M1qrrr_energy(j)
              littlefactors(j,6) = M1qffr_energy(j)

              velocity(j,1) = sum(littlefactors(j,:)*velocity_coeffs(:,1))
              velocity(j,2) = sum(littlefactors(j,:)*velocity_coeffs(:,2))

           enddo

           M1flux_energy_interface = 0.0d0

           !get geometric interpolated distributions functions for Mueller
           logdistro = log(M1en_energy_fluid(1:number_groups)*M1_moment_to_distro(:))
           do j=1,number_groups-1
              loginterface_distroj(j) = logdistro(j) + &
                   (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
                   (logdistro(j+1) - logdistro(j))/ &
                   (nulibtable_logenergies(j+1)-nulibtable_logenergies(j))
              interface_distroj(j) = exp(loginterface_distroj(j))
           enddo
           j=number_groups
           loginterface_distroj(j) = logdistro(j) + &
                (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
                (logdistro(j) - logdistro(j-1))/ &
                (nulibtable_logenergies(j)-nulibtable_logenergies(j-1))
           interface_distroj(j) = exp(loginterface_distroj(j))

           !find xi_i's
           xi(1) = 1.0d0
           do j=2,number_groups
              xi(j) = interface_distroj(j)/(interface_distroj(j)+interface_distroj(j-1))
           enddo

           temp_term = (nulibtable_etop(1)-nulibtable_ebottom(1))/ &
                (1.0d0-nulibtable_energies(1)/nulibtable_energies(1+1))*xi(1)
           FL(1,1) = velocity(1,1)*temp_term
           FL(1,2) = velocity(1,2)*temp_term
           
           M1flux_energy_interface(1,1) = M1flux_energy_interface(1,1) + FL(1,1)/nulibtable_etop(1)
           M1flux_energy_interface(1,2) = M1flux_energy_interface(1,2) + FL(1,2)/nulibtable_etop(1)
           do j=2,number_groups
              temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
                   (1.0d0-nulibtable_energies(j)/nulibtable_energies(j+1))*xi(j)
              FL(j,1) = velocity(j,1)*temp_term
              FL(j,2) = velocity(j,2)*temp_term

              temp_term = (nulibtable_etop(j)-nulibtable_ebottom(j))/ &
                   (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0)*(1.0d0-xi(j))
              FR(j,1) = velocity(j,1)*temp_term
              FR(j,2) = velocity(j,2)*temp_term

              M1flux_energy_interface(j-1,1) = M1flux_energy_interface(j-1,1) + FR(j,1)/nulibtable_ebottom(j)
              M1flux_energy_interface(j-1,2) = M1flux_energy_interface(j-1,2) + FR(j,2)/nulibtable_ebottom(j)
              M1flux_energy_interface(j,1) = M1flux_energy_interface(j,1) + FL(j,1)/nulibtable_etop(j)
              M1flux_energy_interface(j,2) = M1flux_energy_interface(j,2) + FL(j,2)/nulibtable_etop(j)
           enddo

           !find flux differences to get actual flux term.
           j=1
           ep = nulibtable_etop(j)
           M1flux_diff_energy(j,1) = ep*M1flux_energy_interface(j,1)!/ep but not b/c of comment above
           M1flux_diff_energy(j,2) = ep*M1flux_energy_interface(j,2)!/ep but not b/c of comment above

           do j=2,number_groups
              em = nulibtable_ebottom(j)
              ep = nulibtable_etop(j)
              de = (ep-em)           

              M1flux_diff_energy(j,1) = (ep*M1flux_energy_interface(j,1)- &
                   em*M1flux_energy_interface(j-1,1))!/de but not b/c of comment above
              M1flux_diff_energy(j,2) = (ep*M1flux_energy_interface(j,2)- &
                   em*M1flux_energy_interface(j-1,2))!/de but not b/c of comment above

              if (M1flux_diff_energy(j,1).ne.M1flux_diff_energy(j,1)) then
                 write(*,*) M1flux_diff_energy(j,1),j,i,k
                 stop "flux NaNing...1"
              else if (M1flux_diff_energy(j,2).ne.M1flux_diff_energy(j,2)) then
                 write(*,*) M1flux_diff_energy(j,2),j,i,k
                 stop "flux NaNing... 2"
              endif
           enddo

           !set cell fluxes, add to existing value
           flux_M1_energy(k,i,:,1) = dts*implicit_factor*M1flux_diff_energy(:,1)
           flux_M1_energy(k,i,:,2) = dts*implicit_factor*M1flux_diff_energy(:,2)
        enddo
     enddo
     !$OMP END PARALLEL DO! end do
  endif

!################################################################
!################################################################
!###################Explicit Scattering##########################
!################################################################
!################################################################


  if (include_Ielectron_exp) then
     !$OMP PARALLEL DO PRIVATE(i,M1en_energy,M1flux_energy,M1eddy_energy, &
     !$OMP M1pff_energy,M1qrrr_energy,M1qffr_energy,M1chi_energy,local_M, &
     !$OMP local_J,local_H,local_L,ies_sourceterms,h,invalp,invalp2,alp2,onealp, &
     !$OMP invX,invX2,X2,oneX,W2,oneW,v2,onev,invr,invr2,local_u,local_uup, &
     !$OMP local_littleh,local_littlehupup,local_Hdown,local_Ltilde,JoverE, &
     !$OMP JoverF,HoverE,HoverF,LoverE,LoverF,j,j_prime,nucubed,nucubedprime, &
     !$OMP R0out,R1out,R0in,R1in,ies_temp,species_factor)
     do k=ghosts1+1,M1_imaxradii
        do i=1,number_species_to_evolve

           if (i.eq.3.and.number_species.eq.3) then
              species_factor = 4.0d0
           else if (i.eq.3) then
              stop "add in correct species terms"
           else
              species_factor = 1.0d0
           endif
           
           M1en_energy = q_M1(k,i,:,1)/species_factor
           M1flux_energy = q_M1(k,i,:,2)/species_factor
           M1eddy_energy = q_M1(k,i,:,3)
           M1pff_energy = q_M1_extra(k,i,:,1)
           M1qrrr_energy = q_M1_extra(k,i,:,2)/species_factor
           M1qffr_energy = q_M1_extra(k,i,:,3)/species_factor
           M1chi_energy = q_M1_extra(k,i,:,4)

           local_M = 0.0d0
           local_J = 0.0d0
           local_H = 0.0d0
           local_L = 0.0d0
           ies_sourceterms = 0.0d0

           h = (1.0d0+eps(k)+press(k)/rho(k))

           if (GR) then
              invalp = 1.0d0/alp(k)
              invalp2 = invalp**2
              alp2 = alp(k)**2
              onealp = alp(k)
              invX = 1.0d0/X(k)
              invX2 = invX**2
              X2 = X(k)**2
              oneX = X(k)
           else
              invalp = 1.0d0
              invalp2 = 1.0d0
              alp2 = 1.0d0
              onealp = 1.0d0
              invX = 1.0d0
              invX2 = 1.0d0
              X2 = 1.0d0
              oneX = 1.0d0
           endif

           if (v_order.eq.-1) then
              if (GR) then
                 W2 = W(k)**2
                 oneW = W(k)
                 v2 = v(k)**2
                 onev = v(k)
              else
                 W2 = 1.0d0/(1.0d0-v1(k)**2)
                 oneW = sqrt(W2)
                 v2 = v1(k)**2
                 onev = v1(k)

              endif
           else if (v_order.eq.0) then
              W2 = 0.0d0
              oneW = 1.0d0
              v2 = 0.0d0
              onev = 0.0d0
           else
              stop "add in v order"
           endif
           invr = 1.0d0/x1(k)
           invr2 = invr*invr

           local_u(1) = -oneW*onealp
           local_u(2) = oneW*onev*oneX
           local_uup(1) = -local_u(1)*invalp2
           local_uup(2) = local_u(2)*invX2
           
           local_littleh(1,1) = -v2*W2 !1.0d0-local_u(1)*local_u(1)/alp(k)**2
           local_littleh(2,1) = local_u(2)*invX2*local_u(1)
           local_littleh(1,2) = -local_u(1)*invalp2*local_u(2)
           local_littleh(2,2) = W2 !1.0d0+local_u(2)*local_u(2)/X(k)**2

           local_littlehupup(1,1) = -invalp2 + local_uup(1)*local_uup(1)
           local_littlehupup(1,2) = local_uup(1)*local_uup(2)
           local_littlehupup(2,1) = local_uup(2)*local_uup(1)
           local_littlehupup(2,2) = invX2 + local_uup(2)*local_uup(2)

           !linearize the block terms
           do j=1,number_groups
              !M^{ab}
              local_M(1,1) = M1en_energy(j)*invalp2
              local_M(1,2) = M1flux_energy(j)*invX2*invalp
              local_M(2,1) = local_M(1,2)
              local_M(2,2) = M1eddy_energy(j)*M1en_energy(j)*invX2**2

              do ii=1,2
                 do jj=1,2
                    local_J(j) = local_J(j) + local_M(ii,jj)*local_u(ii)*local_u(jj)

                    !H^{a}
                    local_H(j,1) = local_H(j,1) - local_M(ii,jj)*local_u(ii)*local_littleh(1,jj)
                    local_H(j,2) = local_H(j,2) - local_M(ii,jj)*local_u(ii)*local_littleh(2,jj)
                 enddo
              enddo
                   
              do ii=1,2
                 do jj=1,2
                    !L^{ab}
                    local_L(j,1,1) = local_L(j,1,1) + &
                         (local_M(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj))* &
                         (1.5d0*M1chi_energy(j)-0.5d0) + &
                         (local_littlehupup(1,1)*local_J(j)/3.0d0)* &
                         (1.5d0-1.5d0*M1chi_energy(j))
                    local_L(j,1,2) = local_L(j,1,2) + &
                         (local_M(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj))* &
                         (1.5d0*M1chi_energy(j)-0.5d0) + &
                         (local_littlehupup(1,2)*local_J(j)/3.0d0)* &
                         (1.5d0-1.5d0*M1chi_energy(j))
                    local_L(j,2,1) = local_L(j,2,1) + &
                         (local_M(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj))* &
                         (1.5d0*M1chi_energy(j)-0.5d0) + &
                         (local_littlehupup(2,1)*local_J(j)/3.0d0)* &
                         (1.5d0-1.5d0*M1chi_energy(j))
                    local_L(j,2,2) = local_L(j,2,2) + &
                         (local_M(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj))* &
                         (1.5d0*M1chi_energy(j)-0.5d0) + &
                         (local_littlehupup(2,2)*local_J(j)/3.0d0)* &
                         (1.5d0-1.5d0*M1chi_energy(j))
                 enddo
              enddo

              local_Hdown(j,1) = local_H(j,1)*local_littleh(1,1)*(-alp2) + &
                   local_H(j,2)*local_littleh(1,2)*(-alp2)
              local_Hdown(j,2) = local_H(j,1)*local_littleh(2,1)*X2 + &
                   local_H(j,2)*local_littleh(2,2)*X2
              
              !L^a_b tilde as in Shibata et al. 2011
              local_Ltilde(j,1,1) = -local_L(j,1,1)*alp2 - local_J(j)*local_littleh(1,1)*onethird
              local_Ltilde(j,1,2) = local_L(j,1,2)*X2 - local_J(j)*local_littleh(1,2)*onethird
              local_Ltilde(j,2,1) = -local_L(j,2,1)*alp2 - local_J(j)*local_littleh(2,1)*onethird
              local_Ltilde(j,2,2) = local_L(j,2,2)*X2 - local_J(j)*local_littleh(2,2)*onethird

              !get parts of J,H,L that go like E,F for implicit solve
              JoverE(j) = W2*(1.0d0+M1eddy_energy(j)*v2*invX2)
              JoverF(j) = -2.0d0*W2*onev*invX

              HoverE(j,1) = oneW*local_littleh(1,1)*invalp - &
                   local_littleh(1,2)*local_u(2)*M1eddy_energy(j)*invX2**2
              HoverE(j,2) = oneW*local_littleh(2,1)*invalp - &
                   local_littleh(2,2)*local_u(2)*M1eddy_energy(j)*invX2**2
              HoverF(j,1) = oneW*local_littleh(1,2)*invX2 - &
                   local_littleh(1,1)*local_u(2)*invX2*invalp
              HoverF(j,2) = oneW*local_littleh(2,2)*invX2 - &
                   local_littleh(2,1)*local_u(2)*invX2*invalp
              
              LoverE(j,1,1) = local_littleh(1,1)**2*invalp2+ &
                   M1eddy_energy(j)*invX2**2*local_littleh(1,2)**2
              LoverE(j,1,2) = invalp2*local_littleh(1,1)*local_littleh(2,1) + &
                   M1eddy_energy(j)*invX2**2*local_littleh(1,2)*local_littleh(2,2)
              LoverE(j,2,1) = LoverE(j,1,2)
              LoverE(j,2,2) = invalp2*local_littleh(2,1)**2 + &
                   M1eddy_energy(j)*invX2**2*local_littleh(2,2)**2
              LoverF(j,1,1) = 2.0d0*invX2*invalp*local_littleh(1,2)*local_littleh(1,1)
              LoverF(j,1,2) = invX2*invalp*(local_littleh(1,1)*local_littleh(2,2) + &
                   local_littleh(1,2)*local_littleh(2,1))
              LoverF(j,2,1) = LoverF(j,1,2)
              LoverF(j,2,2) = 2.0d0*invX2*invalp*local_littleh(2,1)*local_littleh(2,2)
              
           enddo

           do j=1,number_groups !\omega
              do j_prime=1,number_groups !\omega^\prime, integrate over this
                 
                 !\nu^3 in units of E,F.  If f(fprime)=1,
                 !nucubed(nucubedprime) = J(Jprime). I do not include
                 !the 4\pi seen in Shibata et al. (2011) as my
                 !variables are /srad.  nucubed has the \delta E for
                 !the bin width in it already, so do the E and J,
                 !therefore we don't have to worry about any bin
                 !widths, just add up the individual contributions
                 
                 !J should never be bigger than nu3, if it is we are
                 !over populated, this happens in high velocity,
                 !degenerate conditions.  It never gets too big.
                 
                 !multiply by 4pi because we just integrated over
                 !incoming neutrino. We only integrated over energy, this
                 !is integraing over solid angle (which is done kinda,
                 !but still off by 4\pi)

                 nucubed = M1_moment_to_distro_inverse(j) 
                 nucubedprime = M1_moment_to_distro_inverse(j_prime)

                 R0out = 0.5d0*ies(k,i,j,j_prime,1)
                 R1out = 1.5d0*ies(k,i,j,j_prime,2)
                 if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
                 R0in = 0.5d0*ies(k,i,j_prime,j,1)
                 R1in = 1.5d0*ies(k,i,j_prime,j,2)
                 
                 !Note, I quickly tested the values of these cutoffs.
                 !Turns out, 1e13 works better than 5e12.  It does not
                 !lead to instability in the core but does give better
                 !agreement when the neutrinosphere start climbing in
                 !density.  5e13 gives instability in the core right
                 !when rho_c gets to ~few*10^13
                 if (rho(k).gt.5.0d12*rho_gf) then
                    R0out = 0.5d0*ies(k,i,j,j_prime,1)*(5.0d12*rho_gf/rho(k))**1.5d0
                    R1out = 1.5d0*ies(k,i,j,j_prime,2)*(5.0d12*rho_gf/rho(k))**1.5d0
                    if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
                    R0in = 0.5d0*ies(k,i,j_prime,j,1)*(5.0d12*rho_gf/rho(k))**1.5d0
                    R1in = 1.5d0*ies(k,i,j_prime,j,2)*(5.0d12*rho_gf/rho(k))**1.5d0
                 endif

                 !implicit terms for j equation for energy density,
                 !proto j_prime energy density
                 ies_temp = species_factor*implicit_factor*dts*alp2*4.0d0*pi*( & 
                      ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
                      local_H(j,1))*R0in*JoverE(j_prime) + &
                      R1in*(HoverE(j_prime,1)*(  &
                      (nucubed-local_J(j))*local_littleh(1,1)*onethird - &
                      local_u(1)*(-invalp2)*local_Hdown(j,1) - local_Ltilde(j,1,1)) + &
                      HoverE(j_prime,2)*(  &
                      (nucubed-local_J(j))*local_littleh(1,2)*onethird - &
                      local_u(1)*(-invalp2)*local_Hdown(j,2) - local_Ltilde(j,1,2))) &
                      )*nulibtable_inv_energies(j_prime)

                 ies_sourceterms(j) = ies_sourceterms(j) + ies_temp*M1en_energy(j_prime)

                 !implicit terms for j equation for energy density,
                 !proto j_prime flux density
                 ies_temp = species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                      ((nucubed-local_J(j))*(-local_u(1)*invalp2) - &
                      local_H(j,1))*R0in*JoverF(j_prime) + &
                      R1in*(HoverF(j_prime,1)*(  &
                      (nucubed-local_J(j))*local_littleh(1,1)*onethird - &
                      local_u(1)*(-invalp2)*local_Hdown(j,1) - local_Ltilde(j,1,1)) + &
                      HoverF(j_prime,2)*(  &
                      (nucubed-local_J(j))*local_littleh(1,2)*onethird - &
                      local_u(1)*(-invalp2)*local_Hdown(j,2) - local_Ltilde(j,1,2))) &
                      )*nulibtable_inv_energies(j_prime)

                 ies_sourceterms(j) = ies_sourceterms(j) + &
                      ies_temp*M1flux_energy(j_prime)

                 !implicit terms for j equation for energy density,
                 !proto to j energy density
                 ies_temp = species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                      -R0out*(nucubedprime-local_J(j_prime))* &
                      (JoverE(j)*local_u(1)*(-invalp2)+HoverE(j,1)) + &
                      R1out*(local_Hdown(j_prime,1)*(HoverE(j,1)*local_u(1)*(-invalp2) + &
                      LoverE(j,1,1)) + local_Hdown(j_prime,2)* &
                      (HoverE(j,2)*local_u(1)*(-invalp2) + LoverE(j,1,2))) &
                      )*nulibtable_inv_energies(j_prime)

                 ies_sourceterms(j) = &
                      ies_sourceterms(j) + ies_temp*M1en_energy(j)

                 !implicit terms for j equation for energy density,
                 !proto to j flux density
                 ies_temp = species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                      -R0out*(nucubedprime-local_J(j_prime))* &
                      (JoverF(j)*local_u(1)*(-invalp2)+HoverF(j,1)) + &
                      R1out*(local_Hdown(j_prime,1)*(HoverF(j,1)*local_u(1)*(-invalp2) + &
                      LoverF(j,1,1)) + local_Hdown(j_prime,2)* &
                      (HoverF(j,2)*local_u(1)*(-invalp2) + LoverF(j,1,2))) &
                      )*nulibtable_inv_energies(j_prime)
                 
                 ies_sourceterms(j) = ies_sourceterms(j) + ies_temp*M1flux_energy(j)

                 !now these equations are for the evolution of F
                 !implicit terms for j equation for flux density,
                 !proto j_prime energy density
                 ies_temp = species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                      ((nucubed-local_J(j))*(local_u(2)*invX2) - &
                      local_H(j,2))*R0in*JoverE(j_prime) + &
                      R1in*(HoverE(j_prime,1)*(  &
                      (nucubed-local_J(j))*local_littleh(2,1)*onethird - &
                      local_u(2)*X2*local_Hdown(j,1) - local_Ltilde(j,2,1)) + &
                      HoverE(j_prime,2)*(  &
                      (nucubed-local_J(j))*local_littleh(2,2)*onethird - &
                      local_u(2)*X2*local_Hdown(j,2) - local_Ltilde(j,2,2))) &
                      )*nulibtable_inv_energies(j_prime)
                 
                 ies_sourceterms(j+number_groups) = &
                      ies_sourceterms(j+number_groups) + ies_temp*M1en_energy(j_prime)

                 !implicit terms for j equation for flux density,
                 !proto j_prime flux density
                 ies_temp = species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                      ((nucubed-local_J(j))*(local_u(2)*invX2) - &
                      local_H(j,2))*R0in*JoverF(j_prime) + &
                      R1in*(HoverF(j_prime,1)*(  &
                      (nucubed-local_J(j))*local_littleh(2,1)*onethird - &
                      local_u(2)*X2*local_Hdown(j,1) - local_Ltilde(j,2,1)) + &
                      HoverF(j_prime,2)*(  &
                      (nucubed-local_J(j))*local_littleh(2,2)*onethird - &
                      local_u(2)*X2*local_Hdown(j,2) - local_Ltilde(j,2,2))) &
                      )*nulibtable_inv_energies(j_prime)
                 
                 ies_sourceterms(j+number_groups) = &
                      ies_sourceterms(j+number_groups) + ies_temp*M1flux_energy(j_prime)                  

                 !implicit terms for j equation for flux density,
                 !proto to j energy density
                 ies_temp = species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                      -R0out*(nucubedprime-local_J(j_prime))* &
                      (JoverE(j)*local_u(2)*(invX2)+HoverE(j,2)) + &
                      R1out*(local_Hdown(j_prime,1)*(HoverE(j,1)*local_u(2)*invX2 + &
                      LoverE(j,2,1)) + local_Hdown(j_prime,2)* &
                      (HoverE(j,2)*local_u(2)*invX2 + LoverE(j,2,2))) &
                      )*nulibtable_inv_energies(j_prime)
                 
                 ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
                      ies_temp*M1en_energy(j)

                 !implicit terms for j equation for flux density,
                 !proto to j flux density
                 ies_temp = species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                      -R0out*(nucubedprime-local_J(j_prime))* &
                      (JoverF(j)*local_u(2)*invX2+HoverF(j,2)) + &
                      R1out*(local_Hdown(j_prime,1)*(HoverF(j,1)*local_u(2)*invX2 + &
                      LoverF(j,2,1)) + local_Hdown(j_prime,2)* &
                      (HoverF(j,2)*local_u(2)*invX2 + LoverF(j,2,2))) &
                      )*nulibtable_inv_energies(j_prime)
                 
                 ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) + &
                      ies_temp*M1flux_energy(j)
              enddo
        
              !set cell fluxes, add to existing value
              flux_M1_scatter(k,i,j,1) = ies_sourceterms(j)
              flux_M1_scatter(k,i,j,2) = ies_sourceterms(j+number_groups)        
              
              ies_sourceterm(k,i,j,1) = ies_sourceterms(j)/(implicit_factor*dts*alp2)
              ies_sourceterm(k,i,j,2) = ies_sourceterms(j+number_groups)/(implicit_factor*dts*X2)

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO! end do
  endif

end subroutine M1_explicitterms
