!-*-f90-*-
!Initialize the M1 module
subroutine M1_init
  
  use GR1D_module, only : opacity_table,nulib_opacity_gf,nulib_emissivity_gf, &
       energy_gf, mev_to_erg,M1_maxradii,M1_imaxradii,length_gf,n1,x1,ghosts1, &
       ye,temp,rho_gf,x1i,number_groups,rho,eas,ghosts1,q_M1,nulib_energy_gf, &
       temp_mev_to_kelvin,number_species,volume,pi,M1_moment_to_distro,clite, &
       hbarc_mevcm,M1_testcase_number,v_order,include_nes_kernels, &
       M1_moment_to_distro_inverse,nulib_kernel_gf,number_species_to_evolve, &
       include_epannihil_kernels
  use nulibtable

  implicit none

  integer i

  call nulibtable_reader(opacity_table,include_nes_kernels,include_epannihil_kernels)
  
  !change units of emissivities, opacities, energies and inverse
  !energies to code units, then we only have to do it once these are
  !located in GR1D_module.F90

  !emissivities are in units of erg/cm^3/s/srad (I've multiplied by bin width (in MeV) in the table reader)
  !nulib_emissivity_gf = energy_gf/(length_gf**3*time_gf)

  !opacities are in units of cm^-1
  !nulib_opacity_gf = 1.0d0/length_gf

  !energies are in units of MeV
  !nulib_energy_gf = mev_to_erg*energy_gf

  !kernels are in units of cm^3 s^-1
  !nulib_kernel_gf = length_gf**3/time_gf

  nulibtable_emissivities = max(-200.0d0,log10(10.0d0**(nulibtable_emissivities)*nulib_emissivity_gf))
  nulibtable_absopacity = log10(10.0d0**nulibtable_absopacity*nulib_opacity_gf)
  nulibtable_scatopacity = log10(10.0d0**nulibtable_scatopacity*nulib_opacity_gf)
  if (include_nes_kernels) then
     nulibtable_Itable_Phi0 = log10(10.0d0**nulibtable_Itable_Phi0*nulib_kernel_gf)
     !Phi1 is stored as the ratio of Phi1/Phi0, so no log, no units
  endif
  if (include_epannihil_kernels) then
     nulibtable_epannihiltable_Phi0 = log10(10.0d0**nulibtable_epannihiltable_Phi0*nulib_kernel_gf)
     !Phi1 is stored as the ratio of Phi1/Phi0, so no log, no units
  endif

  !convert energies to reduced units to save time, note nulib_ewidth is NOT converted
  nulibtable_energies = nulibtable_energies*nulib_energy_gf
  nulibtable_ebottom = nulibtable_ebottom*nulib_energy_gf
  nulibtable_etop = nulibtable_etop*nulib_energy_gf
  nulibtable_inv_energies = nulibtable_inv_energies/nulib_energy_gf

  !set the log based version for time savings
  nulibtable_logenergies = log(nulibtable_energies)
  nulibtable_logetop = log(nulibtable_etop)

  !find zone that matches maximum radii
  do i=1,n1
     if (x1(i)/length_gf.lt.M1_maxradii) M1_imaxradii = i
  enddo

  if (v_order.eq.0.or.v_order.eq.-1) then
     write(*,*) "Velocity order is:",v_order," (-1 for all orders)"
  else
     stop "implement v_order"
  endif
  write(*,*) "M1_init: extract radii at", x1(M1_imaxradii)/length_gf, "index:", M1_imaxradii, "of", n1

  !conversion from energy (momentum) density to angle integrated distribution function (*\mu)
  M1_moment_to_distro(:) =  (2.0d0*pi*hbarc_mevcm)**3 / &
       (clite*nulib_emissivity_gf/nulib_opacity_gf * &
       nulibtable_ewidths(:)*mev_to_erg*(nulibtable_energies(:)/nulib_energy_gf)**3)
  M1_moment_to_distro_inverse(:) = 1.0d0/M1_moment_to_distro(:)

  if (M1_imaxradii.gt.n1-ghosts1) stop "M1_init: Your extraction radii is too big"

  if (number_species_to_evolve.eq.-1) then
     number_species_to_evolve = number_species
  endif

  call M1_updateeas

end subroutine M1_init
