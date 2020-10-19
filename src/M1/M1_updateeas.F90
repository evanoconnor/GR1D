!-*-f90-*-
!this subroutine updates the neutrino-matter interaction coefficents
subroutine M1_updateeas
  
  use GR1D_module
  use nulibtable
  implicit none

  real*8 :: xrho,xtemp,xye,xeta,xxp,xxn
  real*8,dimension(3) :: xn_N
  integer :: k,i,j,spec

  real*8 :: tempspectrum(number_species,number_groups,number_eas)
  real*8 :: singlespecies_tempspectrum(number_groups,number_eas)

  real*8 :: inelastic_tempspectrum(number_species,number_groups,number_groups,2)
  real*8 :: singlespecies_inelastic_tempspectrum(number_groups,number_groups,2)

  real*8 :: epannihil_tempspectrum(number_species,number_groups,number_groups,4)
  real*8 :: singlespecies_epannihil_tempspectrum(number_groups,number_groups,4)
  
  real*8 :: bremsstrahlung_tempspectrum(number_species,number_groups,number_groups,2)
  real*8 :: singlespecies_bremsstrahlung_tempspectrum(number_groups,number_groups,2)

  real*8 :: blackbody_spectra(number_species,number_groups)
  real*8 :: energy_x
  
  real*8 :: integrator,integrator2
  integer :: myid
  real,parameter :: m_amu = 931.494061d0 !MeV
  real,parameter :: mev_to_gram = 1.782661758d-27 !one MeV is # grams

  integer :: keytemp,keyerr
  real*8 :: eosdummy(17)

  if (M1_testcase_number.eq.0 .or.M1_testcase_number.eq.1) then

     !$OMP PARALLEL DO PRIVATE(xrho,xtemp,xye,xxn,xxp,xn_N,tempspectrum,singlespecies_tempspectrum, &
     !$OMP keytemp,keyerr,eosdummy,xeta,inelastic_tempspectrum,singlespecies_inelastic_tempspectrum, &
     !$OMP epannihil_tempspectrum,singlespecies_epannihil_tempspectrum,integrator2,integrator, &
     !$OMP bremsstrahlung_tempspectrum,singlespecies_bremsstrahlung_tempspectrum,blackbody_spectra, & 
     !$omp energy_x,i,j,spec)


     
     do k=2,M1_imaxradii+ghosts1-1
    
   !     write(*,*) nulibtable_logtemp_max
    
        xrho = rho(k)/rho_gf
        xtemp = temp(k)
        xye = ye(k)

        tempspectrum = 0.0d0
        if (log10(xrho).lt.nulibtable_logrho_min) then
           tempspectrum = 0.0d0
        else
           if (log10(xtemp).lt.nulibtable_logtemp_min) then
              stop "M1_update_eas: temp too low"
           endif
           if (xye.lt.nulibtable_ye_min) stop "M1_update_eas: ye too low"
           if (log10(xtemp).gt.nulibtable_logtemp_max) stop "M1_update_eas: temp too high"
           if (xye.gt.nulibtable_ye_max) stop "M1_update_eas: ye too high"

           if (number_species_to_evolve.eq.1) then
              call nulibtable_single_species_range_energy(xrho,xtemp,xye,1, &
                   singlespecies_tempspectrum,number_groups,number_eas)
              tempspectrum(1,:,:) = singlespecies_tempspectrum(:,:)
           else if (number_species_to_evolve.eq.3) then
              call nulibtable_range_species_range_energy(xrho,xtemp,xye,tempspectrum, &
                   number_species,number_groups,number_eas)
           else
              stop "set up eas interpolation for this number of species"
           endif
        endif
        
        !set new eas variables
!        if (include_epannihil_kernels) then
!           tempspectrum(3,:,2) = 0.0d0
!           tempspectrum(3,:,1) = 0.0d0
!        endif

        if (.true.) then
           eas(k,:,:,2) = tempspectrum(:,:,2)
           eas(k,:,:,3) = tempspectrum(:,:,3)

           !re calculate emmisivity from black body.
           keytemp = 1 !keep temperature
           keyerr = 0       
           call nuc_eos_full(xrho,xtemp,xye,eosdummy(1),eosdummy(2),eosdummy(3), &
                eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),eosdummy(8), &
                eosdummy(9),eosdummy(10),eosdummy(11),eosdummy(12), &
                eosdummy(13),elechem(k),eosdummy(15),eosdummy(16),eosdummy(17), &
                keytemp,keyerr,eos_rf_prec)
           if(keyerr.ne.0) then
              write(6,*) "############################################"
              write(6,*) "EOS PROBLEM in inelastic scatter :"
              write(6,*) "timestep number: ",nt
              write(6,"(i5,1P3E18.9)") k,xrho,xtemp,xye
              stop "This is bad!"
           endif

           xeta = (elechem(k)-eosdummy(17))/xtemp
           do j=1,number_groups
              energy_x = nulibtable_energies(j)/(nulib_energy_gf*xtemp)
              blackbody_spectra(1,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                   (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x-xeta))
              blackbody_spectra(2,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                   (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x+xeta))
              blackbody_spectra(3,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                   (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x))
           enddo

           eas(k,1,:,1) = (blackbody_spectra(1,:)*tempspectrum(1,:,2)* &
                nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf
           eas(k,2,:,1) = (blackbody_spectra(2,:)*tempspectrum(2,:,2)* &
                nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf
           eas(k,3,:,1) = 4.0d0*(blackbody_spectra(3,:)*tempspectrum(3,:,2)* &
                nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf

        else
           eas(k,:,:,:) = tempspectrum(:,:,:)
        endif

        !get electron chemical potential if scattering/epannihil is turned on &
        ! get number density if bremsstrahlung is turned on
        
        if (include_Ielectron_exp.or.include_Ielectron_imp.or.include_epannihil_kernels &
		.or. include_bremsstrahlung_kernels.or.include_gang_kernels) then
           keytemp = 1 !keep temperature, not resetting any hydro variables
           keyerr = 0       
           call nuc_eos_full(xrho,xtemp,xye,eosdummy(1),eosdummy(2),eosdummy(3), &
                eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),eosdummy(8), &
                eosdummy(9),xxn,xxp,eosdummy(12), &
                eosdummy(13),elechem(k),eosdummy(15),eosdummy(16),eosdummy(17), &
                keytemp,keyerr,eos_rf_prec)

                
           if(keyerr.ne.0) then
              write(6,*) "############################################"
              write(6,*) "EOS PROBLEM in M1_updateeas.F90:"
              write(6,*) "timestep number: ",nt
              write(6,"(i5,1P3E18.9)") k,xrho,xtemp,xye
              stop "This is bad!"
           endif

           xtemp = temp(k)
           xeta = elechem(k)/temp(k)
           
           if (include_bremsstrahlung_kernels) then
	           xn_N(1) = xxn*xrho/(m_amu*mev_to_gram)
			   xn_N(2) = xxp*xrho/(m_amu*mev_to_gram)
               xn_N(3) = SQRT(xxn*xxp)*xrho/(m_amu*mev_to_gram)
           endif
           
        endif
        
        if (include_Ielectron_imp.or.include_Ielectron_exp) then
           inelastic_tempspectrum = 0.0d0
           if (log10(xrho).lt.nulibtable_logrho_min) then
              inelastic_tempspectrum = 0.0d0
           else
              if (log10(xtemp).lt.nulibtable_logItemp_min) stop "M1_update_eas: Itemp too low"
              if (log10(xeta).lt.nulibtable_logIeta_min) then
                 write(*,*) xrho,xtemp,xye,xeta
                 stop "M1_update_eas: Ieta too low"
              endif
              if (log10(xtemp).gt.nulibtable_logItemp_max) then 
                        write(*,*) xtemp, log10(xtemp),nulibtable_logItemp_max
                        stop "M1_update_eas: temp too high"
              endif
              if (log10(xeta).gt.nulibtable_logIeta_max) then
                 write(*,*) xrho,xtemp,xye,xeta,nulibtable_logIeta_max,k
                 stop "M1_update_eas: Ieta too high"
              endif
                
              if (number_species_to_evolve.eq.1) then
                 call nulibtable_inelastic_single_species_range_energy2(xtemp,xeta,1, &
                      singlespecies_inelastic_tempspectrum,number_groups,number_groups,2)
                 inelastic_tempspectrum(1,:,:,:) = singlespecies_inelastic_tempspectrum(:,:,:)
              else if (number_species_to_evolve.eq.3) then
                 call nulibtable_inelastic_range_species_range_energy2(xtemp,xeta, &
                      inelastic_tempspectrum,number_species,number_groups,number_groups,2)
              else
                 stop "set up eas interpolation for this number of species"
              endif


           endif
           
           !set new ies variables
           ies(k,:,:,:,:) = inelastic_tempspectrum(:,:,:,:)

        endif

        !get ep-annihilation kernels for i.eq.3 only, if turned on
        if (include_epannihil_kernels) then
           epannihil_tempspectrum = 0.0d0
           if (log10(xrho).lt.nulibtable_logrho_min) then
              epannihil_tempspectrum = 0.0d0
           else
              if (log10(xtemp).lt.nulibtable_logItemp_min) stop "M1_update_eas: Itemp too low"
              if (log10(xeta).lt.nulibtable_logIeta_min) then
                 write(*,*) xrho,xtemp,xye,xeta
                 stop "M1_update_eas: Ieta too low"
              endif
              if (log10(xtemp).gt.nulibtable_logItemp_max) stop "M1_update_eas: temp too high"
              if (log10(xeta).gt.nulibtable_logIeta_max) then
                 write(*,*) xrho,xtemp,xye,xeta,nulibtable_logIeta_max,k
                 stop "M1_update_eas: Ieta too high"
              endif
                
              !only do this for nux
              i = 3
              call nulibtable_epannihil_single_species_range_energy2(xtemp,xeta,i, &
                   singlespecies_epannihil_tempspectrum,number_groups,number_groups,4)

              epannihil_tempspectrum(i,:,:,:) = singlespecies_epannihil_tempspectrum(:,:,:)

           endif
           
           !set new ep annihilation variables
           epannihil(k,:,:,:,:) = epannihil_tempspectrum(:,:,:,:)

        endif
        
        
        !get bremsstrahlung kernels for i.eq.3 only, if turned on
        if (include_bremsstrahlung_kernels) then
           bremsstrahlung_tempspectrum = 0.0d0
           !if (log10(xrho).lt.nulibtable_logrho_min) then
           if (log10(xrho) .lt. 11.22d0 .or. log10(xrho) .gt. 15.22d0 &
                .or. xtemp .LT. 2.0d0 .or. xtemp .GT. 50.0d0 .or.  &
                xye .GT. 0.5d0) then  
                  bremsstrahlung_tempspectrum = 0.0d0
           else
              if (log10(xtemp).lt.nulibtable_logItemp_min) stop "M1_update_eas: Itemp too low"
              if (ANY(log10(xn_N).lt.nulibtable_logn_N_min)) then
                 write(*,*) xrho,xtemp,xye,xn_N
                 stop "M1_update_eas: n_N too low"
              endif
              if (log10(xtemp).gt.nulibtable_logItemp_max) then 
                 write(*,*) xtemp, log10(xtemp),nulibtable_logItemp_max
                 stop "M1_update_eas: temp too high"
              endif
              if (ANY(log10(xn_N).gt.nulibtable_logn_N_max)) then
                 write(*,*) xrho,xtemp,xye,xn_N,nulibtable_logn_N_max,k
                 stop "M1_update_eas: n_N too high"
              endif
!~                 write(*,*) log10(xrho)
              !only do this for nux
              i = 3
!~               xtemp = 26.0298d0
!~               xn_N=1.54314e+28
              call nulibtable_bremsstrahlung_single_species_range_energy2(xtemp,xn_N,i, &
                   singlespecies_bremsstrahlung_tempspectrum,number_groups,number_groups,2)
                   
                if ( ANY(singlespecies_bremsstrahlung_tempspectrum < 0.0d0) ) then
                    write(*,*) xrho,xtemp,xye,xn_N,nulibtable_logn_N_max,k
                    stop "brem <0"
                endif
!~                 write(*,*)  singlespecies_bremsstrahlung_tempspectrum(11,11,1)
!~                 write(*,*) nulibtable_energies(10),nulibtable_energies(11)
!~                 write(*,*) nulibtable_ewidths(11)
!~                 stop
!~ 			  do spec = 1,3
	              bremsstrahlung_tempspectrum(1,:,:,:) = singlespecies_bremsstrahlung_tempspectrum(:,:,:)
	              bremsstrahlung_tempspectrum(2,:,:,:) = singlespecies_bremsstrahlung_tempspectrum(:,:,:)
	              bremsstrahlung_tempspectrum(3,:,:,:) = singlespecies_bremsstrahlung_tempspectrum(:,:,:)
!~               enddo
           endif
           
           !set new ep annihilation variables
           bremsstrahlung(k,:,:,:,:) = bremsstrahlung_tempspectrum(:,:,:,:)

        endif
        
               
        if (include_gang_kernels) then 
              bremsstrahlung_tempspectrum = 0.0d0
              singlespecies_bremsstrahlung_tempspectrum = 0.0d0
              if (log10(xtemp).lt.nulibtable_logItemp_min) then 
                    write(*,*) xtemp,log10(xtemp)   
                    stop "M1_update_eas: Itemp too low"
              endif 
              if (log10(xye).gt.nulibtable_logYe_min .or. &
                  log10(xye).lt.nulibtable_logYe_max .or. &
                log10(xrho/(m_amu*mev_to_gram)).lt.nulibtable_logn_N_min .or. &
                log10(xrho/(m_amu*mev_to_gram)).gt.nulibtable_logn_N_max .or. &
                log10(xtemp).gt.nulibtable_logItemp_max .or. &
                log10(xtemp).lt.nulibtable_logItemp_min) then 
              else
              
                  if (log10(xrho/(m_amu*mev_to_gram)).lt.nulibtable_logn_N_min) then
                        write(*,*) xrho,xtemp,xye,xrho/(m_amu*mev_to_gram),nulibtable_logn_N_min,k
                         stop "M1_update_eas: n_N too low"
                  endif
                 if (log10(xtemp).gt.nulibtable_logItemp_max) then 
                           write(*,*) xtemp,log10(xtemp)
                        stop "M1_update_eas: temp too high"
                endif 

                if (log10(xrho/(m_amu*mev_to_gram)).gt.nulibtable_logn_N_max) then
                         write(*,*) xrho,xtemp,xye,xrho/(m_amu*mev_to_gram),nulibtable_logn_N_max,k
                        stop "M1_update_eas: n_N too high"
                endif
                i=3
                call nulibtable_bremsstrahlung_gang_single_species_range_energy2(xtemp,&
                                        xrho/(m_amu*mev_to_gram),xYe, &
                            i,singlespecies_bremsstrahlung_tempspectrum,number_groups,number_groups,2)

                bremsstrahlung_tempspectrum(1,:,:,:) =singlespecies_bremsstrahlung_tempspectrum(:,:,:)
                bremsstrahlung_tempspectrum(2,:,:,:) =singlespecies_bremsstrahlung_tempspectrum(:,:,:)
                bremsstrahlung_tempspectrum(3,:,:,:) =singlespecies_bremsstrahlung_tempspectrum(:,:,:)

              endif  
              bremsstrahlung(k,:,:,:,:) =bremsstrahlung_tempspectrum(:,:,:,:)
              
        

              if ( ANY(singlespecies_bremsstrahlung_tempspectrum < 0.0d0) ) then
                    write(*,*)  singlespecies_bremsstrahlung_tempspectrum
                    write(*,*) xrho,xtemp,xye,xn_N,nulibtable_logn_N_max,k
                    stop "brem gang <0"
              endif
        endif
    enddo
     !$OMP END PARALLEL DO! end do
      
  else if (M1_testcase_number.ge.2.and.M1_testcase_number.le.9) then
     !this is taken care of
  else
     stop "add in eas updating code for this test case"
  endif

end subroutine M1_updateeas
