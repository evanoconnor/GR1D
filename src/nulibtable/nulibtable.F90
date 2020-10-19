!-*-f90-*-
module nulibtable

  implicit none

  integer :: nulibtable_number_species
  integer :: nulibtable_number_groups
  
  real*8, allocatable,save :: nulibtable_logrho(:)
  real*8, allocatable,save :: nulibtable_logtemp(:)
  real*8, allocatable,save :: nulibtable_ye(:)
  real*8, allocatable,save :: nulibtable_logItemp(:)
  real*8, allocatable,save :: nulibtable_logIeta(:)
  real*8, allocatable,save :: nulibtable_logn_N(:)
  real*8, allocatable,save :: nulibtable_logYe(:)
  
  
  real*8, allocatable,save :: nulibtable_energies(:)
  real*8, allocatable,save :: nulibtable_inv_energies(:)
  real*8, allocatable,save :: nulibtable_ewidths(:)
  real*8, allocatable,save :: nulibtable_ebottom(:)
  real*8, allocatable,save :: nulibtable_etop(:)

  real*8, allocatable,save :: nulibtable_logenergies(:)
  real*8, allocatable,save :: nulibtable_logetop(:)

  real*8, allocatable,save :: nulibtable_emissivities(:,:,:,:)
  real*8, allocatable,save :: nulibtable_absopacity(:,:,:,:)
  real*8, allocatable,save :: nulibtable_scatopacity(:,:,:,:)

  real*8, allocatable,save :: nulibtable_Itable_Phi0(:,:,:)
  real*8, allocatable,save :: nulibtable_Itable_Phi1(:,:,:)

  real*8, allocatable,save :: nulibtable_epannihiltable_Phi0(:,:,:)
  real*8, allocatable,save :: nulibtable_epannihiltable_Phi1(:,:,:)
  
  real*8, allocatable,save :: nulibtable_bremsstrahlung_Phi0(:,:,:)
  real*8, allocatable,save :: nulibtable_bremsstrahlung_gang_Phi0(:,:,:,:)
!~   real*8, allocatable,save :: gang_table(:,:,:,:)

  integer :: nulibtable_nrho
  integer :: nulibtable_ntemp
  integer :: nulibtable_nye

  integer :: nulibtable_nItemp
  integer :: nulibtable_nIeta
  integer :: nulibtable_nn_N
  integer :: nulibtable_nIYe

  real*8 :: nulibtable_logrho_min
  real*8 :: nulibtable_logrho_max

  real*8 :: nulibtable_logtemp_min
  real*8 :: nulibtable_logtemp_max

  real*8 :: nulibtable_ye_min
  real*8 :: nulibtable_ye_max

  real*8 :: nulibtable_logItemp_min
  real*8 :: nulibtable_logItemp_max

  real*8 :: nulibtable_logIeta_min
  real*8 :: nulibtable_logIeta_max
  
  real*8 :: nulibtable_logn_N_max
  real*8 :: nulibtable_logn_N_min

  real*8 :: nulibtable_logYe_max
  real*8 :: nulibtable_logYe_min

  integer :: nulibtable_number_easvariables

end module nulibtable

!this takes rho,temp,ye,species and energy and return eas
subroutine nulibtable_single_species_single_energy(xrho,xtemp,xye,lns,lng,eas,eas_n1)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: lns, lng
  integer, intent(in) :: eas_n1
  real*8, intent(out) :: eas(eas_n1)
  integer :: startindex,endindex

  if (eas_n1.ne.nulibtable_number_easvariables) stop "supplied array dimensions (1) is not commensurate with table"
  
  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  startindex = (lns-1)*nulibtable_number_groups+(lng-1)+1
  endindex = startindex

  call intp3d_many_mod(xlrho,xltemp,xye,eas(1), &
       nulibtable_emissivities(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  call intp3d_many_mod(xlrho,xltemp,xye,eas(2), &
       nulibtable_absopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  call intp3d_many_mod(xlrho,xltemp,xye,eas(3), &
       nulibtable_scatopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  eas(:) = 10.0d0**eas(:)

end subroutine nulibtable_single_species_single_energy

!this takes rho,temp,ye,species and return eas over energy range
subroutine nulibtable_single_species_range_energy(xrho,xtemp,xye,lns,eas,eas_n1,eas_n2)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: lns
  integer, intent(in) :: eas_n1,eas_n2
  real*8, intent(out) :: eas(eas_n1,eas_n2)
  integer :: ing
  real*8 :: xeas(eas_n1)
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_easvariables) then
     stop "nulibtable_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  startindex = (lns-1)*nulibtable_number_groups+1
  endindex = startindex + nulibtable_number_groups - 1

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_emissivities(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,1) = 10.0d0**xeas(:)

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_absopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,2) = 10.0d0**xeas(:)

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas, &
       nulibtable_scatopacity(:,:,:,startindex:endindex),nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)
  
  eas(:,3) = 10.0d0**xeas(:)

end subroutine nulibtable_single_species_range_energy

!this takes rho,temp,ye and return eas over energy and species range
subroutine nulibtable_range_species_range_energy(xrho,xtemp,xye,eas,eas_n1,eas_n2,eas_n3)
  
  use nulibtable
  implicit none

  real*8, intent(in) :: xrho, xtemp, xye !inputs
  real*8 :: xlrho, xltemp !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ins,ing
  real*8 :: xeas(eas_n1*eas_n2)
  integer :: index

  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_easvariables) then
     stop "nulibtable_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xlrho = log10(xrho)
  xltemp = log10(xtemp)

  if (xlrho.lt.nulibtable_logrho_min) stop "density below nulib table minimum rho"
  if (xlrho.gt.nulibtable_logrho_max) stop "density above nulib table maximum rho"
  if (xltemp.lt.nulibtable_logtemp_min) stop "temperature below nulib table minimum temp"
  if (xltemp.gt.nulibtable_logtemp_max) stop "temperature above nulib table maximum temp"
  if (xye.lt.nulibtable_ye_min) stop "ye below nulib table minimum ye"
  if (xye.gt.nulibtable_ye_max) stop "ye above nulib table maximum ye"

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_emissivities,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,1) = 10.0d0**xeas(index)
     enddo
  enddo

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_absopacity,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,2) = 10.0d0**xeas(index)
     enddo
  enddo

  xeas = 0.0d0
  call intp3d_many_mod(xlrho,xltemp,xye,xeas,nulibtable_scatopacity,nulibtable_nrho, &
       nulibtable_ntemp,nulibtable_nye,eas_n1*eas_n2,nulibtable_logrho, &
       nulibtable_logtemp,nulibtable_ye)

  do ins=1,nulibtable_number_species
     do ing=1,nulibtable_number_groups
        index = (ins-1)*nulibtable_number_groups + (ing-1) + 1
        eas(ins,ing,3) = 10.0d0**xeas(index)
     enddo
  enddo

end subroutine nulibtable_range_species_range_energy

!this takes temp,eta, and return phi0/1 over energy (both in and out) and species range
subroutine nulibtable_inelastic_range_species_range_energy2(xtemp,xeta, &
     eas,eas_n1,eas_n2,eas_n3,eas_n4)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta !inputs
  real*8 :: xltemp, xleta !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3,eas_n4
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3,eas_n4)
  integer :: ins,ing_in,ing_out
  real*8 :: xeas(eas_n1*eas_n2*(eas_n2+1)/2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55


  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_inelastic_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_inelastic_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_groups) then
     stop "nulibtable_inelastic_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif
  if(size(eas,4).ne.2) then
     stop "nulibtable_inelastic_range_species_range_energy: supplied array dimensions (4) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  xleta = log10(xeta)

  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib inelastic table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib inelastic table maximum temp"
  if (xleta.lt.nulibtable_logIeta_min) stop "eta below nulib inelastic table minimum eta"
  if (xleta.gt.nulibtable_logIeta_max) stop "eta above nulib inelastic table maximum eta"

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_Itable_Phi0,nulibtable_nItemp, &
       nulibtable_nIeta,eas_n1*eas_n2*(eas_n2+1)/2,nulibtable_logItemp, &
       nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_in=1,nulibtable_number_groups
        do ing_out=1,ing_in
           index = index + 1
           eas(ins,ing_in,ing_out,1) = 10.0d0**xeas(index)
        enddo
     enddo
     do ing_in=1,nulibtable_number_groups
        do ing_out=ing_in+1,nulibtable_number_groups
           eas(ins,ing_in,ing_out,1) =  &
                exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in))/ &
                (xtemp*energy_conversion))*eas(ins,ing_out,ing_in,1) !cernohorsky 94
        enddo
     enddo
  enddo

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_Itable_Phi1,nulibtable_nItemp, &
       nulibtable_nIeta,eas_n1*eas_n2*(eas_n2+1)/2,nulibtable_logItemp, &
       nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_in=1,nulibtable_number_groups
        do ing_out=1,ing_in
           index = index + 1
           !this way because we interpolate \phi1/\phi0
           eas(ins,ing_in,ing_out,2) = xeas(index)*eas(ins,ing_in,ing_out,1)
        enddo
     enddo
     do ing_in=1,nulibtable_number_groups
        do ing_out=ing_in+1,nulibtable_number_groups
           eas(ins,ing_in,ing_out,2) =  &
                exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in))/ &
                (xtemp*energy_conversion))*eas(ins,ing_out,ing_in,2) !cernohorsky 94
        enddo
     enddo
  enddo  

  
  
end subroutine nulibtable_inelastic_range_species_range_energy2

!this takes temp,eta, and return phi0/1 over energy (both in and out) range for a single species
subroutine nulibtable_inelastic_single_species_range_energy2(xtemp,xeta, &
     lns,eas,eas_n1,eas_n2,eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta !inputs
  real*8 :: xltemp, xleta !log versions
  integer, intent(in) :: lns,eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ing_in,ing_out
  real*8 :: xeas(eas_n1*(eas_n1+1)/2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_inelastic_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif  
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_inelastic_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,3).ne.2) then
     stop "nulibtable_inelastic_single_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  xleta = log10(xeta)

  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib inelastic table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib inelastic table maximum temp"
  if (xleta.lt.nulibtable_logIeta_min) stop "eta below nulib inelastic table minimum eta"
  if (xleta.gt.nulibtable_logIeta_max) stop "eta above nulib inelastic table maximum eta"

  startindex = (lns-1)*nulibtable_number_groups*(nulibtable_number_groups+1)/2+1
  endindex = startindex + nulibtable_number_groups*(nulibtable_number_groups+1)/2 - 1

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_Itable_Phi0(:,:,startindex:endindex), &
       nulibtable_nItemp,nulibtable_nIeta,eas_n1*(eas_n1+1)/2, &
       nulibtable_logItemp,nulibtable_logIeta)

  index = 0
  do ing_in=1,nulibtable_number_groups
     do ing_out=1,ing_in
        index = index + 1
        eas(ing_in,ing_out,1) = 10.0d0**xeas(index)
     enddo
  enddo
  do ing_in=1,nulibtable_number_groups
     do ing_out=ing_in+1,nulibtable_number_groups
        eas(ing_in,ing_out,1) =  &
             exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in))/ &
             (xtemp*energy_conversion))*eas(ing_out,ing_in,1) !cernohorsky 94
     enddo
  enddo

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_Itable_Phi1(:,:,startindex:endindex), &
       nulibtable_nItemp,nulibtable_nIeta,eas_n1*(eas_n1+1)/2, &
       nulibtable_logItemp,nulibtable_logIeta)

  index = 0
  do ing_in=1,nulibtable_number_groups
     do ing_out=1,ing_in
        index = index + 1
        !this way because we interpolate \phi1/\phi0
        eas(ing_in,ing_out,2) = xeas(index)*eas(ing_in,ing_out,1)
     enddo
  enddo
  do ing_in=1,nulibtable_number_groups
     do ing_out=ing_in+1,nulibtable_number_groups
        eas(ing_in,ing_out,2) =  &
             exp(-(nulibtable_energies(ing_out)-nulibtable_energies(ing_in))/ &
             (xtemp*energy_conversion))*eas(ing_out,ing_in,2) !cernohorsky 94
     enddo
  enddo

  
  
end subroutine nulibtable_inelastic_single_species_range_energy2

!this takes temp,eta, and return phi0/1 for ep-annihilation over energy (both neutrinos) and species range
subroutine nulibtable_epannihil_range_species_range_energy2(xtemp,xeta, &
     eas,eas_n1,eas_n2,eas_n3,eas_n4)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta !inputs
  real*8 :: xltemp, xleta !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3,eas_n4
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3,eas_n4)
  integer :: ins,ing_this,ing_that
  real*8 :: xeas(eas_n1*eas_n2*eas_n3*2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55


  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_epannihil_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_epannihil_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_groups) then
     stop "nulibtable_enannihil_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif
  if(size(eas,4).ne.4) then
     stop "nulibtable_epannihil_range_species_range_energy: supplied array dimensions (4) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  xleta = log10(xeta)

  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib epannihil table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib epannihil table maximum temp"
  if (xleta.lt.nulibtable_logIeta_min) stop "eta below nulib epannihil table minimum eta"
  if (xleta.gt.nulibtable_logIeta_max) stop "eta above nulib epannihil table maximum eta"

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_epannihiltable_Phi0,nulibtable_nItemp, &
       nulibtable_nIeta,eas_n1*eas_n2*eas_n2*2,nulibtable_logItemp, &
       nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_this=1,nulibtable_number_groups
        do ing_that=1,nulibtable_number_groups
           index = index + 1
           eas(ins,ing_this,ing_that,1) = 10.0d0**xeas(index)
           index = index + 1
           eas(ins,ing_this,ing_that,2) = 10.0d0**xeas(index)
!           eas(ins,ing_this,ing_that,1) = eas(ins,ing_this,ing_that,2)* &
!                exp(-(nulibtable_energies(ing_this)+nulibtable_energies(ing_that))/(xtemp*energy_conversion))

        enddo
     enddo
  enddo

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_epannihiltable_Phi1,nulibtable_nItemp, &
       nulibtable_nIeta,eas_n1*eas_n2*eas_n2*2,nulibtable_logItemp, &
       nulibtable_logIeta)

  index = 0
  do ins=1,nulibtable_number_species
     do ing_this=1,nulibtable_number_groups
        do ing_that=1,nulibtable_number_groups
           !this way because we interpolate \phi1/\phi0
           index = index + 1
           eas(ins,ing_this,ing_that,3) = xeas(index)*eas(ins,ing_this,ing_that,1)
           index = index + 1
           eas(ins,ing_this,ing_that,4) = xeas(index)*eas(ins,ing_this,ing_that,2)
!           eas(ins,ing_this,ing_that,3) = eas(ins,ing_this,ing_that,4)* &
!                exp(-(nulibtable_energies(ing_this)+nulibtable_energies(ing_that))/(xtemp*energy_conversion))
        enddo
     enddo
  enddo  

  
  
end subroutine nulibtable_epannihil_range_species_range_energy2

!this takes temp,eta, and return epannihilation phi0/1 over energy (both nus) range for a single species
subroutine nulibtable_epannihil_single_species_range_energy2(xtemp,xeta, &
     lns,eas,eas_n1,eas_n2,eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xeta !inputs
  real*8 :: xltemp, xleta !log versions
  integer, intent(in) :: lns,eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ing_this,ing_that
  real*8 :: xeas(eas_n1*eas_n2*2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_epannihil_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif  
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_epannihil_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,3).ne.4) then
     stop "nulibtable_epannihil_single_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  xleta = log10(xeta)

  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib epannihil table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib epannihil table maximum temp"
  if (xleta.lt.nulibtable_logIeta_min) stop "eta below nulib epannihil table minimum eta"
  if (xleta.gt.nulibtable_logIeta_max) stop "eta above nulib epannihil table maximum eta"

  startindex = (lns-1)*nulibtable_number_groups*nulibtable_number_groups*2+1
  endindex = startindex + nulibtable_number_groups*nulibtable_number_groups*2 - 1

  xeas = 0.0d0
  
!~   write(*,*) MAXVAL(nulibtable_epannihiltable_Phi0),MINVAL(nulibtable_epannihiltable_Phi0)
  
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_epannihiltable_Phi0(:,:,startindex:endindex), &
       nulibtable_nItemp,nulibtable_nIeta,eas_n1*eas_n2*2,nulibtable_logItemp,nulibtable_logIeta)

  index = 0
  do ing_this=1,nulibtable_number_groups
     do ing_that=1,nulibtable_number_groups
        index = index + 1
        eas(ing_this,ing_that,1) = 10.0d0**xeas(index)
        index = index + 1
        eas(ing_this,ing_that,2) = 10.0d0**xeas(index)
!        eas(ing_this,ing_that,1) = eas(ing_this,ing_that,2)* &
!             (exp(-(nulibtable_energies(ing_this)+nulibtable_energies(ing_that))/(xtemp*energy_conversion)))
     enddo
  enddo

  xeas = 0.0d0
  call intp2d_many_mod(xltemp,xleta,xeas,nulibtable_epannihiltable_Phi1(:,:,startindex:endindex), &
       nulibtable_nItemp,nulibtable_nIeta,eas_n1*eas_n2*2,nulibtable_logItemp,nulibtable_logIeta)

  index = 0
  do ing_this=1,nulibtable_number_groups
     do ing_that=1,nulibtable_number_groups
        index = index + 1
        !this way because we interpolate \phi1/\phi0
        eas(ing_this,ing_that,3) = xeas(index)*eas(ing_this,ing_that,1)
        index = index + 1
        eas(ing_this,ing_that,4) = xeas(index)*eas(ing_this,ing_that,2)
!        eas(ing_this,ing_that,3) = eas(ing_this,ing_that,4)*&
!             (exp(-(nulibtable_energies(ing_this)+nulibtable_energies(ing_that))/(xtemp*energy_conversion)))
     enddo
  enddo
  
end subroutine nulibtable_epannihil_single_species_range_energy2

!this takes temp,eta, and return phi0/1 for ep-annihilation over energy (both neutrinos) and species range
subroutine nulibtable_bremsstrahlung_range_species_range_energy2(xtemp,xn_N, &
     eas,eas_n1,eas_n2,eas_n3,eas_n4)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xn_N(3) !inputs
  real*8 :: xltemp, xln_N !log versions
  integer, intent(in) :: eas_n1,eas_n2,eas_n3,eas_n4
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3,eas_n4)
  integer :: ins,ing_this,ing_that,i
  real*8 :: xeas(eas_n1*eas_n2*eas_n3*2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55


  if(size(eas,1).ne.nulibtable_number_species) then
     stop "nulibtable_bremsstrahlung_range_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_bremsstrahlung_range_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif  
  if(size(eas,3).ne.nulibtable_number_groups) then
     stop "nulibtable_bermsstrahlung_range_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif
  if(size(eas,4).ne.2) then
     stop "nulibtable_bremsstrahlung_range_species_range_energy: supplied array dimensions (4) is not commensurate with table"
  endif
  
  xltemp = log10(xtemp)
  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib bremsstrahlung table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib bremsstrahlung table maximum temp"  
  
  do i=1,3
  

	  xln_N = log10(xn_N(i))
	

	  if (xln_N.lt.nulibtable_logn_N_min) stop "n_N below nulib bremsstrahlung table minimum n_N"
	  if (xln_N.gt.nulibtable_logn_N_max) stop "n_N above nulib bremsstrahlung table maximum n_N"
	
	  xeas = 0.0d0
	  call intp2d_many_mod(xltemp,xln_N,xeas,nulibtable_bremsstrahlung_Phi0,nulibtable_nItemp, &
	       nulibtable_nn_N,eas_n1*eas_n2*eas_n2*2,nulibtable_logItemp, &
	       nulibtable_logn_N)
	
	  index = 0
	  do ins=1,nulibtable_number_species
	     do ing_this=1,nulibtable_number_groups
	        do ing_that=1,nulibtable_number_groups
				if ( i .EQ. 3) then
		           index = index + 1
		           eas(ins,ing_this,ing_that,1) = eas(ins,ing_this,ing_that,1) + &
													28.0d0/3.0d0*10.0d0**xeas(index)
		           index = index + 1
		           eas(ins,ing_this,ing_that,2) = eas(ins,ing_this,ing_that,2) + &
													28.0d0/3.0d0*10.0d0**xeas(index)
		        else 
		           index = index + 1
		           eas(ins,ing_this,ing_that,1) = eas(ins,ing_this,ing_that,1) + &
										           10.0d0**xeas(index)
		           index = index + 1
		           eas(ins,ing_this,ing_that,2) = eas(ins,ing_this,ing_that,2) + &
										           10.0d0**xeas(index)
		        endif
	        enddo
	     enddo
	  enddo
  enddo
  
  
end subroutine nulibtable_bremsstrahlung_range_species_range_energy2

!this takes temp,n_N, and return bremsstrahlungation phi0 over energy (both nus) range for a single species
subroutine nulibtable_bremsstrahlung_single_species_range_energy2(xtemp,xn_N, &
     lns,eas,eas_n1,eas_n2,eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xn_N(3) !inputs
  real*8 :: xltemp, xln_N !log versions
  integer, intent(in) :: lns,eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ing_this,ing_that,i
  real*8 :: xeas(eas_n1*eas_n2*2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif  
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,3).ne.2) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib bremsstrahlung table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib bremsstrahlung table maximum temp"  
  
  eas = 0.0d0
  
  do i =1,3
	  xln_N = log10(xn_N(i))
	

	  if (xln_N.lt.nulibtable_logn_N_min) stop "n_N below nulib bremsstrahlung table minimum n_N"
	  if (xln_N.gt.nulibtable_logn_N_max) stop "n_N above nulib bremsstrahlung table maximum n_N"
	
	  startindex = (lns-1)*nulibtable_number_groups*nulibtable_number_groups*2+1
	  endindex = startindex + nulibtable_number_groups*nulibtable_number_groups*2 - 1
!~         write(*,*) MAXVAL(nulibtable_bremsstrahlung_Phi0),MINVAL(nulibtable_bremsstrahlung_Phi0)
	  xeas = 0.0d0
	  call intp2d_many_mod(xltemp,xln_N,xeas,nulibtable_bremsstrahlung_Phi0(:,:,startindex:endindex), &
	       nulibtable_nItemp,nulibtable_nn_N,eas_n1*eas_n2*2,nulibtable_logItemp,nulibtable_logn_N)
!~         write(*,*) MAXVAL(xeas),MINVAL(xeas)
!~     stop
	  index = 0
	  do ing_this=1,nulibtable_number_groups
	     do ing_that=1,nulibtable_number_groups
		     if (i .eq. 3) then 
		        index = index + 1
		        eas(ing_this,ing_that,1) = eas(ing_this,ing_that,1) + &
											28.0d0/3.0d0*10.0d0**xeas(index)
		        index = index + 1
		        eas(ing_this,ing_that,2) = eas(ing_this,ing_that,2) + &
									        28.0d0/3.0d0*10.0d0**xeas(index)
		     else
		        index = index + 1
		        eas(ing_this,ing_that,1) = eas(ing_this,ing_that,1) + &
									        10.0d0**xeas(index)
		        index = index + 1
		        eas(ing_this,ing_that,2) = eas(ing_this,ing_that,2) + &
                                           10.0d0**xeas(index)
		     endif
	     enddo
	  enddo
  enddo

  if (ANY(eas <0.0d0)) then
    stop "nulib : eas_brem <0"
  endif

   
  
end subroutine nulibtable_bremsstrahlung_single_species_range_energy2
  
!this takes temp,n_N, and return gang bremsstrahlungation phi0 over energy (both nus) range for a single species
subroutine nulibtable_bremsstrahlung_gang_single_species_range_energy2(xtemp,xn_N,xYe, &
     lns,eas,eas_n1,eas_n2,eas_n3)

  use nulibtable
  implicit none

  real*8, intent(in) :: xtemp, xn_N,xYe !inputs
  real*8 :: xltemp, xln_N,xlYe !log versions
  integer, intent(in) :: lns,eas_n1,eas_n2,eas_n3
  real*8, intent(out) :: eas(eas_n1,eas_n2,eas_n3)
  integer :: ing_this,ing_that
  real*8 :: xeas(eas_n1*eas_n2*2)
  integer :: index
  real*8 :: energy_conversion = 1.60217733d-6*5.59424238d-55
  integer :: startindex,endindex

  if(size(eas,1).ne.nulibtable_number_groups) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (1) is not commensurate with table"
  endif  
  if(size(eas,2).ne.nulibtable_number_groups) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (2) is not commensurate with table"
  endif
  if(size(eas,3).ne.2) then
     stop "nulibtable_bremsstrahlung_single_species_range_energy: supplied array dimensions (3) is not commensurate with table"
  endif

  xltemp = log10(xtemp)
  if (xltemp.lt.nulibtable_logItemp_min) stop "temp below nulib bremsstrahlung table minimum temp"
  if (xltemp.gt.nulibtable_logItemp_max) stop "temp above nulib bremsstrahlung table maximum temp"  
  
  eas = 0.0d0
  xln_N = log10(xn_N)
  xlYe = log10(xYe)
	

  if (xln_N.lt.nulibtable_logn_N_min) stop "n_N below nulib bremsstrahlung table minimum n_N"
  if (xln_N.gt.nulibtable_logn_N_max) stop "n_N above nulib bremsstrahlung table maximum n_N"

  if (xlYe.gt.nulibtable_logYe_min) then 
       write(*,*) xlYe,nulibtable_logYe_min  
       stop "Ye below nulib bremsstrahlung table minimum Ye"
  endif
  if (xlYe.lt.nulibtable_logYe_max) then 
       write(*,*) xlYe,nulibtable_logYe_max
         stop "Ye above nulib bremsstrahlung table maximum Ye"
  endif
  startindex = (lns-1)*nulibtable_number_groups*nulibtable_number_groups*2+1
  endindex = startindex + nulibtable_number_groups*nulibtable_number_groups*2 - 1
!~         write(*,*) MAXVAL(nulibtable_bremsstrahlung_Phi0),MINVAL(nulibtable_bremsstrahlung_Phi0)
  xeas = 0.0d0
  call intp3d_many_mod(xltemp,xln_N,xlYe,xeas,nulibtable_bremsstrahlung_gang_Phi0(:,:,:,startindex:endindex), &
       nulibtable_nItemp,nulibtable_nn_N,nulibtable_nIYe,eas_n1*eas_n2*2,&
       nulibtable_logItemp,nulibtable_logn_N,nulibtable_logYe)
       
  index = 0
  do ing_this=1,nulibtable_number_groups
     do ing_that=1,nulibtable_number_groups
                index = index + 1
                eas(ing_this,ing_that,1) = eas(ing_this,ing_that,1) + &
                                           10.0d0**xeas(index)
                index = index + 1
                eas(ing_this,ing_that,2) = eas(ing_this,ing_that,2) + &
                                           10.0d0**xeas(index)

     enddo
  enddo


  if (ANY(eas <0.0d0)) then
    stop "nulib : eas_brem <0"
  endif

   
  
end subroutine nulibtable_bremsstrahlung_gang_single_species_range_energy2
  
!~ subroutine gang_table_range_energy(xtemp,xn_N,xYe, &
!~  brem_array)
	
!~ 	use nulibtable
	
!~ 	implicit none
	
	
!~ 	real*8, intent(in) :: xtemp,xn_N,xYe
!~ 	real*8, intent(out) :: brem_array(nulibtable_number_groups,nulibtable_number_groups,2)
	
!~ 	integer :: i,j,indx
!~ 	real*8 :: temp_array(25) = (/(i, i=2,50, 2)/)
!~ 	real*8 :: Ye_array(26) = (/(i, i=0,50,2)/)/100.0d0
!~ 	real*8 :: n_array(37)
!~ 	real*8 :: n_fem,dx
!~ 	real*8 :: om_array(40) = (/(i,i = 1,40,1)/)/10.0d0 -1.4d0
!~ 	real*8 :: omega
!~ 	real*8 :: eas(40)
!~     real*8,parameter :: nulib_energy_gf = 1.60217733d-6*5.59424238d-55
!~ 	real*8, parameter :: G_f = 1.16637d-11*(1.97326966d-11)  ! Mev-1  cm 
!~     real*8,parameter :: nulib_kernel_gf = 6.77140812d-06**3/2.03001708d+05 
!~ 	real*8, parameter :: hbarc_mevcm = 1.97326966d-11
!~     real*8, parameter :: c_light = 29979245800.0d0 
!~ 	n_array(1:9) = (/(i, i=1,9, 1)/)/1.0d4
!~ 	n_array(10:18) = (/(i, i=1,9, 1)/)/1.0d3
!~ 	n_array(19:27) = (/(i, i=1,9, 1)/)/1.0d2
!~ 	n_array(28:37) = (/(i, i=1,10, 1)/)/1.0d1
	
!~ 	eas = 0.0d0
!~ 	n_fem = xn_N /(1.0d13)**3 
!~ 	if ( n_fem .GT. maxval(n_array) .OR. n_fem .LT. minval(n_array)) then 

!~ 	else
!~ 		call intp3d_gang_table( xtemp, n_fem, xYe, eas, gang_table, 25, 37, 26, 40, temp_array, n_array, Ye_array)
!~ 	endif
!~ 	brem_array = 0.0d0
!~ 	do i = 1,nulibtable_number_groups
!~ 		do j = 1, nulibtable_number_groups
!~ 			omega = log10((nulibtable_energies(i)+nulibtable_energies(j))/nulib_energy_gf)
!~ 			if ( omega .LT. maxval(om_array) .AND. omega .GT. minval(om_array)) then
!~ 				dx = 10.0d0
!~ 				indx= 1+ int((omega-om_array(1))*dx)
!~ 				brem_array(i,j,1) = eas(indx) * xn_N*(1.26d0/2.0d0)**2*G_f**2*3.0d0 &  ! production
!~ 									* exp(-10**(omega)/xtemp)*nulib_kernel_gf*hbarc_mevcm**3*c_light
!~ 				brem_array(i,j,2) = eas(indx) * xn_N*(1.26d0/2.0d0)**2*G_f**2*3.0d0*nulib_kernel_gf&
!~ 									*hbarc_mevcm**3*c_light!annihilation					
!~ 			endif 
!~ 		enddo
!~ 	enddo
!~ 	stop
!~ end subroutine 
