! -*-f90-*-
subroutine allocate_vars

  use GR1D_module

  allocate(x1(n1))
  allocate(x1i(n1))

  allocate(rho(n1))
  allocate(rhop(n1))
  allocate(rhom(n1))

  allocate(v1(n1))
  allocate(v1p(n1))
  allocate(v1m(n1))

  if(do_rotation) then
     allocate(vphi1(n1))
     allocate(vphi1p(n1))
     allocate(vphi1m(n1))
     allocate(vphi(n1))
     allocate(vphip(n1))
     allocate(vphim(n1))
     allocate(omega(n1))
  endif
  allocate(ToverW(n1))

  allocate(v(n1))
  allocate(vp(n1))
  allocate(vm(n1))
  allocate(v_prev(n1))

  allocate(eps(n1))
  allocate(epsp(n1))
  allocate(epsm(n1))

  allocate(energy_nu(n1))
  allocate(mom_nu(n1))

  allocate(eps_kin(n1))
  allocate(binding_energy(n1))

  allocate(mass(n1))
  allocate(mass1(n1))
  allocate(volume(n1))

  allocate(press(n1))
  allocate(pressth(n1))
  allocate(pressp(n1))
  allocate(pressm(n1))
  allocate(press_nu(n1))
  allocate(dnupdr(n1))
  
  allocate(cs2(n1))
  allocate(cs2p(n1))
  allocate(cs2m(n1))
  
  allocate(temp(n1))

  allocate(entropy(n1))
  allocate(nuchem(n1))
  allocate(elechem(n1))
  allocate(massfrac_p(n1))
  allocate(massfrac_n(n1))
  allocate(massfrac_h(n1))
  allocate(massfrac_a(n1))
  allocate(massfrac_abar(n1))
  allocate(massfrac_zbar(n1))

  allocate(ye(n1),ye_prev(n1))
  allocate(dyedt_hydro(n1),dyedt_neutrino(n1))
  allocate(depsdt(n1),dyedt(n1))
  allocate(yep(n1))
  allocate(yem(n1))
  allocate(ynu(n1))

  allocate(q(n1,n_cons))
  allocate(qold(n1,n_cons))
  allocate(qp(n1,n_cons))
  allocate(qm(n1,n_cons))

  allocate(q_hat(n1,n_cons))
  allocate(q_hat_old(n1,n_cons))

  !M1 arrays
  allocate(eas(n1,number_species,number_groups,number_eas)) !emissivity,absorpive crosssection, scattering crosssection
  allocate(ies(n1,number_species,number_groups,number_groups,2)) ! scattering kernel \phis
  allocate(ies_sourceterm(n1,number_species,number_groups,2)) !source terms for matter
  allocate(epannihil(n1,number_species,number_groups,number_groups,4)) ! scattering kernel \phis
  allocate(epannihil_sourceterm(n1,number_species,number_groups,2)) !source terms for matter
  allocate(q_M1(n1,number_species,number_groups,3)) !conserved variables
  allocate(q_M1_prev(n1,number_species,number_groups,3)) !conserved variables
  allocate(q_M1p(n1,number_species,number_groups,3,2)) !conserved variables
  allocate(q_M1m(n1,number_species,number_groups,3,2)) !conserved variables
  allocate(q_M1_extra(n1,number_species,number_groups,4)) !closure variables
  allocate(q_M1_extrap(n1,number_species,number_groups,1,2)) !closure variables
  allocate(q_M1_extram(n1,number_species,number_groups,1,2)) !closure variables
  allocate(q_M1_old(n1,number_species,number_groups,3)) !conserved variables
  allocate(B_M1(n1,number_species,number_groups,3)) !conserved variables
  allocate(C_M1(n1,number_species,number_groups,3)) !conserved variables
  allocate(D_M1(n1,number_species,number_groups,3)) !conserved variables
  allocate(flux_M1(n1,number_species,number_groups,3)) !conserved variables
  allocate(flux_M1_energy(n1,number_species,number_groups,3)) !conserved variables
  allocate(flux_M1_scatter(n1,number_species,number_groups,3)) !conserved variables
  allocate(q_M1_fluid(n1,number_species,number_groups,3)) !neutrino variables in the fluid frame
  allocate(M1_matter_source(n1,4)) !matter source terms
  allocate(M1_moment_to_distro(number_groups))
  allocate(M1_moment_to_distro_inverse(number_groups))

  allocate(sqrt_gamma(n1))

  allocate(flux_diff(n1,n_cons))
  allocate(gravsource(n1,n_cons))
  allocate(presssource(n1,n_cons))
  allocate(coolingsource(n1,n_cons))
  allocate(denergyloss(n1))
  
  allocate(atmo(n1))
  
! #############################################
! GR VARIABLES

  allocate(alp(n1))
  allocate(alpp(n1))
  allocate(alpm(n1))
  allocate(phi(n1))
  allocate(phii(n1))
  allocate(dphidr(n1))
  allocate(X(n1))
  allocate(Xp(n1))
  allocate(Xm(n1))
  allocate(W(n1))
  allocate(Wp(n1))
  allocate(Wm(n1))
  allocate(mgrav(n1))
  allocate(mgravi(n1))
  
end subroutine allocate_vars
