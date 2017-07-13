!-*-f90-*-
subroutine output_all(modeflag)

  use GR1D_module
  use nulibtable
  implicit none

  character*1024 filename
  character*256 basename
  integer modeflag

  integer i,tempghosts
  real*8, allocatable :: cs(:)
  real*8 alp_ana(n1),rho_ana(n1),vel_ana(n1),X_ana(n1)
  real*8 maxr
  
  integer, parameter :: nscalars0 = 256
  real*8 scalars(nscalars0)
  integer nscalars
  integer ishock_radius(1)

  real*8 luminosity(number_species)
  real*8 luminosity_fluid(number_species)
  real*8 luminosity_rad(n1,number_species,2)
  
  real*8 enden_rad(n1,number_species,2)
  real*8 fluxden_rad(n1,number_species,2)
  
  real*8 num_luminosity(number_species)
  real*8 num_luminosity_fluid(number_species)
  real*8 num_luminosity_rad(n1,number_species,2)

  real*8 average_energy(number_species)
  real*8 average_energy_fluid(number_species)
  real*8 average_energy_rad(n1,number_species,2)

  real*8 rms_energy(number_species)
  real*8 rms_energy_fluid(number_species)
  real*8 rms_energy_rad(n1,number_species,2)

  real*8 spectrum(number_groups)

  real*8 fluxfactor_enweighted(n1,number_species)
  real*8 fluxfactor_fluxweighted(n1,number_species)
  real*8 eddingtonfactor_enweighted(n1,number_species)
  real*8 eddingtonfactor_fluxweighted(n1,number_species)
  real*8 total_nu_energy,total_matter_energy,total_energy
  
  integer k,j

  real*8 tau(n1,2) !scattering+absorption/absorption
  real*8 ave_tau(n1,number_species)
  logical nusphere_not_found(number_species,number_groups,2)
  logical ave_nusphere_not_found(number_species)
  real*8 nusphere(number_species,number_groups,2)
  real*8 ave_nusphere(number_species,5)

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then
     
     filename = trim(adjustl(outdir))//"/v1.xg"
     if (.not.small_output) call output_single(v1*clite,filename)
     
     if(do_rotation) then
        filename = trim(adjustl(outdir))//"/omega.xg"
        call output_single(omega*time_gf,filename)
        filename = trim(adjustl(outdir))//"/ToverW.xg"
        call output_single(ToverW,filename)        
        
        if(GR) then
           filename = trim(adjustl(outdir))//"/vphi.xg"
           call output_single(vphi,filename)
        else
           filename = trim(adjustl(outdir))//"/vphi1.xg"
           call output_single(vphi1,filename)
        endif
     endif
     
     filename = trim(adjustl(outdir))//"/rho.xg"
     call output_single(rho/rho_gf,filename)
     
     if (do_nupress.or.do_M1) then
        filename = trim(adjustl(outdir))//"/nuchem.xg"
        if (.not.small_output) call output_single(nuchem,filename)
        filename = trim(adjustl(outdir))//"/press_nu.xg"
        if (.not.small_output) call output_single(press_nu/press_gf,filename)
        filename = trim(adjustl(outdir))//"/energy_nu.xg"
        if (.not.small_output) call output_single(energy_nu/press_gf,filename) !use press_gf b/c neutrino energy is not specific
        filename = trim(adjustl(outdir))//"/dnupdr.xg"
        if (.not.small_output) call output_single(dnupdr/press_gf*length_gf,filename)
     endif

     filename = trim(adjustl(outdir))//"/ye.xg"
     call output_single(ye,filename)
     if (do_M1) then
        filename = trim(adjustl(outdir))//"/dyedt_hydro.xg"
        if (.not.small_output) call output_single(dyedt_hydro*time_gf,filename)
        filename = trim(adjustl(outdir))//"/depsdt.xg"
        if (.not.small_output) call output_single(depsdt,filename)
        filename = trim(adjustl(outdir))//"/ynu.xg"
        if (.not.small_output) call output_single(ynu,filename)
     endif
     
     
     filename = trim(adjustl(outdir))//"/press.xg"
     call output_single(press/press_gf,filename)
     
     filename = trim(adjustl(outdir))//"/eps.xg"
     if (.not.small_output) call output_single(eps/eps_gf,filename)
     
     if (GR) then
        filename = trim(adjustl(outdir))//"/mass_grav.xg"
        call output_single(mgrav/mass_gf,filename)
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call output_single(mass/mass_gf,filename)
     else
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call output_single(mass/mass_gf,filename)
     endif
     
     if (eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/xn.xg"
        if (.not.small_output) call output_single(massfrac_n,filename)
        
        filename = trim(adjustl(outdir))//"/xp.xg"
        if (.not.small_output) call output_single(massfrac_p,filename)
        
        filename = trim(adjustl(outdir))//"/xa.xg"
        if (.not.small_output) call output_single(massfrac_a,filename)
        
        filename = trim(adjustl(outdir))//"/xh.xg"
        if (.not.small_output) call output_single(massfrac_h,filename)
        
        filename = trim(adjustl(outdir))//"/xabar.xg"
        if (.not.small_output) call output_single(massfrac_abar,filename)
        
        filename = trim(adjustl(outdir))//"/xzbar.xg"
        if (.not.small_output) call output_single(massfrac_zbar,filename)
     endif

     if (eoskey.eq.1) then
        filename = trim(adjustl(outdir))//"/pressth.xg"
        if (.not.small_output) call output_single(pressth/press_gf,filename)
     endif

     filename = trim(adjustl(outdir))//"/eps_kin.xg"
     if (.not.small_output) call output_single(eps_kin/eps_gf,filename)

     filename = trim(adjustl(outdir))//"/cs.xg"
     allocate(cs(n1))
     cs(:) = sqrt(cs2(:))*clite
     if (.not.small_output) call output_single(cs,filename)
     deallocate(cs)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           filename = trim(adjustl(outdir))//"/alpha_analytic.xg"
           call output_singlemod(alp_ana,filename,maxr)
           filename = trim(adjustl(outdir))//"/rho_analytic.xg"
           call output_single(rho_ana/rho_gf,filename)
           filename = trim(adjustl(outdir))//"/vel_analytic.xg"
           call output_singlemod(vel_ana*clite,filename,maxr)
           filename = trim(adjustl(outdir))//"/alphamod.xg"
           call output_singlemod(alp,filename,maxr)
           filename = trim(adjustl(outdir))//"/X_analytic.xg"
           call output_single(X_ana,filename)
        endif
        filename = trim(adjustl(outdir))//"/alpha.xg"
        if (.not.small_output) call output_single(alp,filename)
        filename = trim(adjustl(outdir))//"/X.xg"
        if (.not.small_output) call output_single(X,filename)
        filename = trim(adjustl(outdir))//"/W.xg"
        if (.not.small_output) call output_single(W,filename)
        filename = trim(adjustl(outdir))//"/v.xg"
        call output_single(v*clite,filename)
     endif

     if (do_effectivepotential) then
        filename = trim(adjustl(outdir))//"/alpha.xg"
        if (.not.small_output) call output_single(alp,filename)
     endif

     if (do_M1) then
        luminosity_rad = 0.0d0
        enden_rad = 0.0d0
        fluxden_rad = 0.0d0
        average_energy_rad = 0.0d0
        rms_energy_rad = 0.0d0

        do k=1,number_species
           do i=ghosts1+1,M1_imaxradii
              luminosity_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2))*(4.0d0*pi)**2*x1(i)**2*(clite**5/ggrav)
              enden_rad(i,k,1) = sum(q_M1_fluid(i,k,:,1))*4.0d0*pi*length_gf**3/energy_gf
              fluxden_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2))*4.0d0*pi*length_gf**3/energy_gf
              num_luminosity_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2)/nulibtable_energies(:))* &
                   (4.0d0*pi)**2*x1(i)**2*(clite**3/ggrav*mass_gf)
              average_energy_rad(i,k,1) = sum(q_M1_fluid(i,k,:,1))/ &
                   sum(q_M1_fluid(i,k,:,1)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
              rms_energy_rad(i,k,1) = sqrt(sum(q_M1_fluid(i,k,:,1)*nulibtable_energies(:))/ &
                   sum(q_M1_fluid(i,k,:,1)/nulibtable_energies(:))) / (mev_to_erg*energy_gf)

              luminosity_rad(i,k,2) = sum(q_M1(i,k,:,2))*(4.0d0*pi)**2*x1(i)**2*(clite**5/ggrav)
              enden_rad(i,k,2) = sum(q_M1(i,k,:,1))*4.0d0*pi*length_gf**3/energy_gf
              fluxden_rad(i,k,2) = sum(q_M1(i,k,:,2))*4.0d0*pi*length_gf**3/energy_gf
              num_luminosity_rad(i,k,2) = sum(q_M1(i,k,:,2)/nulibtable_energies(:))* &
                   (4.0d0*pi)**2*x1(i)**2*(clite**3/ggrav*mass_gf)
              average_energy_rad(i,k,2) = sum(q_M1(i,k,:,1))/ &
                   sum(q_M1(i,k,:,1)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
              rms_energy_rad(i,k,2) = sqrt(sum(q_M1(i,k,:,1)*nulibtable_energies(:))/ &
                   sum(q_M1(i,k,:,1)/nulibtable_energies(:))) / (mev_to_erg*energy_gf)
           enddo
        enddo

        tau = 0.0d0
        ave_tau = 0.0d0
        do k=1,number_species
           do i=M1_imaxradii,ghosts1+1,-1
              ave_tau(i,k) = ave_tau(i+1,k) + (x1i(i+1)-x1i(i))/ &
                   (sum(q_M1(i,k,:,1)/(eas(i,k,:,2)+eas(i,k,:,3)))/sum(q_M1(i,k,:,1)))
           enddo
        enddo

        fluxfactor_enweighted = 0.0d0
        fluxfactor_fluxweighted = 0.0d0
        eddingtonfactor_enweighted = 0.0d0
        eddingtonfactor_fluxweighted = 0.0d0
        do k=1,number_species
           do i=ghosts1+1,M1_imaxradii
              do j=1,number_groups
                 fluxfactor_enweighted(i,k) = fluxfactor_enweighted(i,k) &
                      + q_M1(i,k,j,1)*q_M1(i,k,j,2)/q_M1(i,k,j,1)
                 fluxfactor_fluxweighted(i,k) = fluxfactor_fluxweighted(i,k) &
                      + q_M1(i,k,j,2)*q_M1(i,k,j,2)/q_M1(i,k,j,1)
                 eddingtonfactor_enweighted(i,k) = eddingtonfactor_enweighted(i,k) &
                      + q_M1(i,k,j,1)*q_M1(i,k,j,3)
                 eddingtonfactor_fluxweighted(i,k) = eddingtonfactor_fluxweighted(i,k) &
                      + q_M1(i,k,j,2)*q_M1(i,k,j,3)
              enddo
              fluxfactor_enweighted(i,k) = fluxfactor_enweighted(i,k) & 
                   / sum(q_M1(i,k,:,1))
              fluxfactor_fluxweighted(i,k) = fluxfactor_fluxweighted(i,k) &
                   / sum(q_M1(i,k,:,2))
              eddingtonfactor_enweighted(i,k) = eddingtonfactor_enweighted(i,k) & 
                   / sum(q_M1(i,k,:,1)) / X(k)**2
              eddingtonfactor_fluxweighted(i,k) = eddingtonfactor_fluxweighted(i,k) &
                   / sum(q_M1(i,k,:,2)) / X(k)**2
           enddo
        enddo

        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nue.xg"
        if (.not.small_output) call output_single(fluxfactor_enweighted(:,1),filename)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_anue.xg"
        if (.not.small_output) call output_single(fluxfactor_enweighted(:,2),filename)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nux.xg"
        if (.not.small_output) call output_single(fluxfactor_enweighted(:,3),filename)
        
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nue.xg"
        if (.not.small_output) call output_single(fluxfactor_fluxweighted(:,1),filename)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_anue.xg"
        if (.not.small_output) call output_single(fluxfactor_fluxweighted(:,2),filename)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nux.xg"
        if (.not.small_output) call output_single(fluxfactor_fluxweighted(:,3),filename)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nue.xg"
        if (.not.small_output) call output_single(eddingtonfactor_enweighted(:,1),filename)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_anue.xg"
        if (.not.small_output) call output_single(eddingtonfactor_enweighted(:,2),filename)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nux.xg"
        if (.not.small_output) call output_single(eddingtonfactor_enweighted(:,3),filename)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nue.xg"
        if (.not.small_output) call output_single(eddingtonfactor_fluxweighted(:,1),filename)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_anue.xg"
        if (.not.small_output) call output_single(eddingtonfactor_fluxweighted(:,2),filename)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nux.xg"
        if (.not.small_output) call output_single(eddingtonfactor_fluxweighted(:,3),filename)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino.xg"
        if (.not.small_output) call output_single(dyedt_neutrino*time_gf,filename)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_fluid_rad.xg"
        call output_single(luminosity_rad(:,1,1),filename)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_lab_rad.xg"
        call output_single(luminosity_rad(:,1,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nue_enden_lab_rad.xg"
        if (.not.small_output) call output_single(enden_rad(:,1,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxden_lab_rad.xg"
        if (.not.small_output) call output_single(fluxden_rad(:,1,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call output_single(num_luminosity_rad(:,1,1),filename)
        
        filename = trim(adjustl(outdir))//"/M1_nue_aveenergy_fluid_rad.xg"
        call output_single(average_energy_rad(:,1,1),filename)

        filename = trim(adjustl(outdir))//"/M1_nue_rmsenergy_fluid_rad.xg"
        call output_single(rms_energy_rad(:,1,1),filename)


        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_fluid_rad.xg"
        call output_single(luminosity_rad(:,2,1),filename)

        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_lab_rad.xg"
        call output_single(luminosity_rad(:,2,2),filename)

        filename = trim(adjustl(outdir))//"/M1_anue_enden_lab_rad.xg"
        if (.not.small_output) call output_single(enden_rad(:,2,2),filename)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxden_lab_rad.xg"
        if (.not.small_output) call output_single(fluxden_rad(:,2,2),filename)

        filename = trim(adjustl(outdir))//"/M1_anue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call output_single(num_luminosity_rad(:,2,1),filename)
        
        filename = trim(adjustl(outdir))//"/M1_anue_aveenergy_fluid_rad.xg"
        call output_single(average_energy_rad(:,2,1),filename)

        filename = trim(adjustl(outdir))//"/M1_anue_rmsenergy_fluid_rad.xg"
        call output_single(rms_energy_rad(:,2,1),filename)


        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_fluid_rad.xg"
        call output_single(luminosity_rad(:,3,1),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_lab_rad.xg"
        call output_single(luminosity_rad(:,3,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_enden_lab_rad.xg"
        if (.not.small_output) call output_single(enden_rad(:,3,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxden_lab_rad.xg"
        if (.not.small_output) call output_single(fluxden_rad(:,3,2),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_numluminosity_fluid_rad.xg"
        if (.not.small_output) call output_single(num_luminosity_rad(:,3,1),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_aveenergy_fluid_rad.xg"
        call output_single(average_energy_rad(:,3,1),filename)

        filename = trim(adjustl(outdir))//"/M1_nux_rmsenergy_fluid_rad.xg"
        call output_single(rms_energy_rad(:,3,1),filename)

            
        filename = trim(adjustl(outdir))//"/M1_nue_ng1_rad.xg"
        if (.not.small_output) call output_single(q_M1(:,1,1,1),filename)

     endif

     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy.xg"
        if (.not.small_output) call output_single(entropy,filename)
        filename = trim(adjustl(outdir))//"/temperature.xg"
        call output_single(temp,filename)
     endif
     
  else if(modeflag.eq.2) then
     
     if (initial_data.eq.'Collapse') then
        !Shock radius
        if (bounce) then
           filename = trim(adjustl(outdir))//"/shock_radius_t.dat"
           call output_scalar(shock_radius/length_gf,filename)
           filename = trim(adjustl(outdir))//"/binding_energy_total.dat"
           call output_scalar(binding_energy_total/energy_gf,filename)
        endif
        
        ! Mass values
        filename = trim(adjustl(outdir))//"/mgrav_Xmax.dat"
        call output_scalar(mgravX,filename)        
        
        filename = trim(adjustl(outdir))//"/mbary_shock.dat"
        call output_scalar(mass(ishock(1)),filename)       
        
        filename = trim(adjustl(outdir))//"/mgrav_shock.dat"
        call output_scalar(mgrav(ishock(1)),filename)                
        
        filename = trim(adjustl(outdir))//"/mgrav_rho1e12.dat"
        call output_scalar(mgrav12,filename)      

        filename = trim(adjustl(outdir))//"/mbary_Xmax.dat"
        call output_scalar(mbaryX,filename)      
        
        filename = trim(adjustl(outdir))//"/mbary_rho1e12.dat"
        call output_scalar(mbary12,filename)     
        
        filename = trim(adjustl(outdir))//"/r_Xmax.dat"
        call output_scalar(rXmax/length_gf,filename)      
        
        filename = trim(adjustl(outdir))//"/r_rho1e12.dat"
        call output_scalar(r12max/length_gf,filename)    

        filename = trim(adjustl(outdir))//"/r_rho1e11.dat"
        call output_scalar(r11max/length_gf,filename)    
        
        filename = trim(adjustl(outdir))//"/M_innercore.dat"
        call output_scalar(mass_inner_core,filename)
        
        !rotation scalars
        if (do_rotation) then
           filename = trim(adjustl(outdir))//"/total_angular_momentum.dat"
           call output_scalar(angular_momentum/mass_gf*time_gf/length_gf**2,filename)
           filename = trim(adjustl(outdir))//"/ToverW_edge.dat"
           call output_scalar(ToverW(n1-ghosts1-1),filename)
        endif
     endif


     if (do_M1) then
        luminosity = 0.0d0 
        num_luminosity = 0.0d0 
        average_energy = 0.0d0
        rms_energy = 0.0d0
        luminosity_fluid = 0.0d0 
        num_luminosity_fluid = 0.0d0 
        average_energy_fluid = 0.0d0
        rms_energy_fluid = 0.0d0

         do k=1,number_species
           luminosity(k)     = sum(q_M1(M1_iextractradii,k,:,2))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**5/ggrav)
           luminosity_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**5/ggrav)
           num_luminosity_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**3/ggrav*mass_gf)
           average_energy_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2))&
                / sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
           rms_energy_fluid(k)     = sqrt( sum(q_M1_fluid(M1_iextractradii,k,:,2)*nulibtable_energies(:))&
                / sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:)) ) / (mev_to_erg*energy_gf)
           average_energy(k) = average_energy_fluid(k)/sqrt(1.0d0-v(M1_iextractradii)**2)*(1.0d0+v(M1_iextractradii))
           rms_energy(k)     = rms_energy_fluid(k)/sqrt(1.0d0-v(M1_iextractradii)**2)*(1.0d0+v(M1_iextractradii))
           num_luminosity(k) = luminosity(k)/(average_energy(k)*mev_to_erg)
        enddo

        total_nu_energy = 0.0d0
        do i=ghosts1+1,M1_imaxradii
           total_nu_energy = total_nu_energy + sum(q_M1(i,1,:,1)*4.0d0*pi*volume(i)/energy_gf) + &
                sum(q_M1(i,2,:,1)*4.0d0*pi*volume(i)/energy_gf) + &
                sum(q_M1(i,3,:,1)*4.0d0*pi*volume(i)/energy_gf)
        enddo

        tau = 0.0d0
        ave_tau = 0.0d0
        nusphere_not_found = .true.
        ave_nusphere_not_found = .true.
        do k=1,number_species
           do i=M1_imaxradii,ghosts1+1,-1
              ave_tau(i,k) = ave_tau(i+1,k) + (x1i(i+1)-x1i(i))/ &
                   (sum(q_M1(i,k,:,1)/(eas(i,k,:,2)+eas(i,k,:,3)))/sum(q_M1(i,k,:,1)))
              if (ave_nusphere_not_found(k).and.ave_tau(i,k).gt.0.66666666d0) then
                 ave_nusphere_not_found(k) = .false.
                 ave_nusphere(k,1) = x1i(i)/length_gf
                 ave_nusphere(k,2) = rho(i)/rho_gf
                 ave_nusphere(k,3) = temp(i)
                 ave_nusphere(k,4) = ye(i)
                 ave_nusphere(k,5) = entropy(i)
              endif
           enddo

           do j=1,number_groups
              tau = 0.0d0
              do i=M1_imaxradii,ghosts1+1,-1
                 tau(i,1) = tau(i+1,1) + (x1i(i+1)-x1i(i))*(eas(i,k,j,2)+eas(i,k,j,3))
                 tau(i,2) = tau(i+1,2) + (x1i(i+1)-x1i(i))*eas(i,k,j,2)
                 if (nusphere_not_found(k,j,1).and.tau(i,1).gt.0.66666666d0) then
                    nusphere_not_found(k,j,1) = .false.
                    nusphere(k,j,1) = x1i(i)
                 endif
                 if (nusphere_not_found(k,j,2).and.tau(i,2).gt.0.66666666d0) then
                    nusphere_not_found(k,j,2) = .false.
                    nusphere(k,j,2) = x1i(i)
                 endif
              enddo
           enddo
        enddo

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = total_nu_energy
        scalars(2) = total_energy_radiated
        scalars(3) = total_energy_absorped
        filename = trim(adjustl(outdir))//"/M1_energies.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = luminosity(1)
        scalars(2) = luminosity(2)
        scalars(3) = luminosity(3)
        filename = trim(adjustl(outdir))//"/M1_flux_lum.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = luminosity_fluid(1)
        scalars(2) = luminosity_fluid(2)
        scalars(3) = luminosity_fluid(3)
        filename = trim(adjustl(outdir))//"/M1_flux_lum_fluid.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = num_luminosity(1)
        scalars(2) = num_luminosity(2)
        scalars(3) = num_luminosity(3)
        filename = trim(adjustl(outdir))//"/M1_flux_numlum.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = num_luminosity_fluid(1)
        scalars(2) = num_luminosity_fluid(2)
        scalars(3) = num_luminosity_fluid(3)
        filename = trim(adjustl(outdir))//"/M1_flux_numlum_fluid.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)


        scalars(1:nscalars0) = 0.0d0
        nscalars = 4
        scalars(1) = total_net_heating/(luminosity(1)+luminosity(2))            
        scalars(2) = luminosity(1)+luminosity(2)
        scalars(3) = total_net_heating
        scalars(4) = total_mass_gain
        filename = trim(adjustl(outdir))//"/M1_net_heating.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = number_groups*number_species*2
        if (nscalars.gt.256) stop "increase scalar count"
        do k=1,number_species
           do j=1,number_groups
              scalars((k-1)*number_groups*2+(j-1)*2+1) = nusphere(k,j,1)/length_gf
              scalars((k-1)*number_groups*2+(j-1)*2+2) = nusphere(k,j,2)/length_gf
           enddo
        enddo
        filename = trim(adjustl(outdir))//"/M1_nusphere_allE.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = number_species*5
        if (nscalars.gt.256) stop "increase scalar count"
        do k=1,number_species
           scalars((k-1)*5+1) = ave_nusphere(k,1)
           scalars((k-1)*5+2) = ave_nusphere(k,2)
           scalars((k-1)*5+3) = ave_nusphere(k,3)
           scalars((k-1)*5+4) = ave_nusphere(k,4)
           scalars((k-1)*5+5) = ave_nusphere(k,5)
        enddo
        filename = trim(adjustl(outdir))//"/M1_nusphere_ave.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = average_energy(1)
        scalars(2) = average_energy(2)
        scalars(3) = average_energy(3)
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_lab.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = rms_energy_fluid(1)
        scalars(2) = rms_energy_fluid(2)
        scalars(3) = rms_energy_fluid(3)
        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_fluid.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = average_energy_fluid(1)
        scalars(2) = average_energy_fluid(2)
        scalars(3) = average_energy_fluid(3)
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_fluid.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = rms_energy(1)
        scalars(2) = rms_energy(2)
        scalars(3) = rms_energy(3)
        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_lab.dat"
        call output_many_scalars(scalars,nscalars0,nscalars,filename)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_out.xg"
        spectrum = q_M1(M1_iextractradii,1,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_out.xg"
        spectrum = q_M1(M1_iextractradii,2,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_out.xg"
        spectrum = q_M1(M1_iextractradii,3,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_out.xg"
        spectrum = q_M1_fluid(M1_iextractradii,1,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_out.xg"
        spectrum = q_M1_fluid(M1_iextractradii,2,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_out.xg"
        spectrum = q_M1_fluid(M1_iextractradii,3,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_cen.xg"
        spectrum = q_M1(ghosts1+1,1,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_cen.xg"
        spectrum = q_M1(ghosts1+1,2,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_cen.xg"
        spectrum = q_M1(ghosts1+1,3,:,2)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_cen.xg"
        spectrum = q_M1_fluid(ghosts1+1,1,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_cen.xg"
        spectrum = q_M1_fluid(ghosts1+1,2,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)
        
        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_cen.xg"
        spectrum = q_M1_fluid(ghosts1+1,3,:,1)*M1_moment_to_distro(:)
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/capturing_factors.xg"
        spectrum = alp(ghosts1+1)*(1.0d0-eas(ghosts1+1,1,:,2)* &
             q_M1_fluid(ghosts1+1,1,:,1)/eas(ghosts1+1,1,:,1))
        call output_spectra(spectrum,filename)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino_c_t.dat"
        if (.not.small_output)  call output_central(dyedt_neutrino*time_gf,filename)

     endif

     ! central values
     filename = trim(adjustl(outdir))//"/rho_c_t.dat"
     call output_central(rho/rho_gf,filename)
        
     filename = trim(adjustl(outdir))//"/ye_c_t.dat"
     call output_central(ye,filename)

     filename = trim(adjustl(outdir))//"/dyedt_hydro_c_t.dat"
     if (.not.small_output) call output_central(dyedt_hydro*time_gf,filename)

     filename = trim(adjustl(outdir))//"/ynu_c_t.dat"
     call output_central(ynu,filename)        

     filename = trim(adjustl(outdir))//"/csound_c_t.dat"
     call output_central(sqrt(cs2),filename)
     
     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy_c_t.dat"
        call output_central(entropy,filename)
        filename = trim(adjustl(outdir))//"/temperature_c_t.dat"
        call output_central(temp,filename)
     endif
     
     filename = trim(adjustl(outdir))//"/totalmass.dat"
     call output_scalar(totalmass/mass_gf,filename)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           filename = trim(adjustl(outdir))//"/alpha_analytic_c_t.dat"
           call output_central(alp_ana,filename)
           filename = trim(adjustl(outdir))//"/rho_analytic_c_t.dat"
           call output_central(rho_ana/rho_gf,filename)
           filename = trim(adjustl(outdir))//"/vel_analytic_c_t.dat"
           call output_central(vel_ana*clite,filename)
        endif
        filename = trim(adjustl(outdir))//"/alpha_c_t.dat"
        call output_central(alp,filename)
        filename = trim(adjustl(outdir))//"/time_c.dat"
        call output_scalar(time_c,filename)
     endif
     
     if(initial_data.eq."Sedov") then
        ishock_radius = maxloc(abs(v1))
        filename = trim(adjustl(outdir))//"/Sedov_radius.dat"
        call output_scalar(x1(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_velocity.dat"
        call output_scalar(v1(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_press.dat"
        call output_scalar(press(ishock_radius(1)),filename)
        filename = trim(adjustl(outdir))//"/Sedov_density.dat"
        call output_scalar(rho(ishock_radius(1)),filename)
        
     endif
     
     if (initial_data.eq.'Collapse') then
        filename = trim(adjustl(outdir))//"/accretion_rates.dat"
        call output_accretion(accretion_rates,filename)
        filename = trim(adjustl(outdir))//"/accreted_mass.dat"
        call output_accretion(accreted_mass,filename)
     endif
     
  endif
  
end subroutine output_all

! *******************************************************************
subroutine output_accretion(var,filename)
      
  use GR1D_module, only: time
  
  implicit none
  real*8 var(*)
  character(*) filename
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time, var(1), var(2), var(3), var(4), &
       var(5), var(6), var(7), var(8), var(9), var(10), var(11)
  
  close(666)
  
end subroutine output_accretion
! ******************************************************************
subroutine output_single(var,filename)
  
  use GR1D_module, only: vs_mass,x1,mass,time,n1,length_gf
  
  implicit none
  real*8 var(n1)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  write(666,*) '"Time = ',time
  
  if(vs_mass) then
     do i=1,n1
        if (abs(var(i)).lt.1.0d-90) then
           write(666,"(1P20E18.9)") mass(i),0.0d0
        else
           write(666,"(1P20E18.9)") mass(i),var(i)
        endif
     enddo
  else
     do i=1,n1
        if (abs(var(i)).lt.1.0d-90) then
           write(666,"(1P20E18.9)") x1(i)/length_gf,0.0d0
        else
           write(666,"(1P20E18.9)") x1(i)/length_gf,var(i)
        endif
     enddo
  endif

  write(666,*) " "
  write(666,*) " "
  close(666)
  
end subroutine output_single

! ******************************************************************
    subroutine output_spectra(var,filename)
      
      use GR1D_module, only: time,number_groups,nulib_energy_gf
      use nulibtable, only: nulibtable_energies

      implicit none
      real*8 var(number_groups)
      character(len=100) filename
      integer nt
      integer i

      open(unit=666,file=trim(adjustl(filename)),status="unknown",&
           form='formatted',position="append")
      write(666,*) '"Time = ',time

      do i=1,number_groups
         if (abs(var(i)).lt.1.0d-90) then
            write(666,"(1P20E18.9)") nulibtable_energies(i)/(nulib_energy_gf),0.0d0
         else
            write(666,"(1P20E18.9)") nulibtable_energies(i)/(nulib_energy_gf),var(i)
         endif
      enddo
      write(666,*) " "
      write(666,*) " "
      close(666)



    end subroutine output_spectra

! ******************************************************************

subroutine output_central(var,filename)
  
  use GR1D_module, only: x1,time,n1,ghosts1
  
  implicit none
  real*8 var(n1)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time,var(ghosts1+1)
  
  close(666)
  
end subroutine output_central

! *******************************************************************
subroutine output_scalar(var,filename)
  
  use GR1D_module, only: time
  implicit none
  real*8 var
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P20E18.9)") time,var
  
  close(666)

end subroutine output_scalar
! *******************************************************************
subroutine output_many_scalars(var,n0,n,filename)
  
  use GR1D_module, only: time
  implicit none
  integer n,n0
  real*8 var(n0)
  character(len=100) filename
  integer nt
  integer i
  
  open(unit=666,file=filename,status="unknown",form='formatted',position="append")
  
  write(666,"(1P64E18.9)") time,var(1:n)
  
  close(666)
  
end subroutine output_many_scalars
! *******************************************************************
subroutine generate_filename(varname,outdir,time,nt,suffix,fname)
  
  implicit none
  
  real*8 time
  integer nt
  character(*) varname
  character(len=128) outdir
  character*(*) suffix
  character*(*) fname
  character*(400) aa
  character(len=100) outtime
  character(len=20) outnt
  integer i,ii
  
  aa=" "
  fname=" "
  write(outnt,"(i10.10)") nt
  
  fname = trim(adjustl(outdir))//"/"//trim(adjustl(varname))//"_nt_"//outnt
  write(outtime,"(f9.7)") time
  
  fname = trim(adjustl(fname))//"_time_"//trim(adjustl(outtime))//".dat"
  
end subroutine generate_filename
! ******************************************************************
subroutine output_singlemod(var,filename,maxr)
  
  use GR1D_module, only: x1,time,n1,length_gf,ghosts1
  
  implicit none
  real*8 var(*)
  character(*) filename
  integer nt
  real*8 maxr
  integer i
  
  open(unit=666,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  write(666,*) '"Time = ',time
  
  do i=ghosts1,n1
     if (x1(i).lt.maxr) then
        write(666,"(1P20E18.9)") x1(i),var(i)
     endif
  enddo
  write(666,*) " "
  write(666,*) " "
  close(666)
  
end subroutine output_singlemod

! ******************************************************************
