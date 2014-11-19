!-*-f90-*-
module leakage_rosswog

  use GR1D_module, only : GR,n1,ghosts1,alp,X,W,rho,temp,ye,cs2,entropy,eps,press, &
       dt,coolingsource,v,v1,length_gf,x1,x1i,rho_gf,eps_gf,press_gf,time_gf, &
       nt,radial_zones, time,t_bounce, massfrac_n,massfrac_p,massfrac_a, &
       massfrac_h, massfrac_abar,massfrac_zbar,pi,clite,small_output,outdir, &
       heat_fac,erg_to_mev,mev_to_erg,do_heating,massn_cgs,eos_rf_prec,avo, &
       do_leak_ros,ishock,shock_radius,sqrt_gamma,do_NNBrem
       
  use eosmodule, only : energy_shift,eos_yemin,eos_yemax
  implicit none


!!!!!!!!!!!!terminology!!!!!!!!!!!!
!
! neutrino species ---- electron type (e)     - 1
!                       antielectron type (a) - 2
!                       all others (x)        - 3
!
! matter species   ---- neutrons (n)          - 1
!                       protones (p)          - 2
!                       heavy nuclei (h)      - 3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !number emission rate used in determining change in ye
  real*8,allocatable :: R_eff(:,:) !effective number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3
  real*8,allocatable :: R_loc(:,:) !local number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3 
  real*8,allocatable :: R_diff(:,:) !diffusive number emission rate per volume, indices <radial zones:neutrino species> - number / sec / cm^3 
  !energy emission rate used in determining change in energy
  real*8,allocatable :: Q_eff(:,:) !effective energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3 
  real*8,allocatable :: Q_loc(:,:) !local energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3 
  real*8,allocatable :: Q_diff(:,:) !diffusive energy emission rate per volume, indices <radial zones:neutrino species> - MeV / sec / cm^3 
  !optical depth
  real*8,allocatable :: chi(:,:) !optical depth with E^2 factored out, indices <radial zones:neutrino species> - 1 / MeV^2
  real*8,allocatable :: zeta(:,:) !Eq. A21, mean free path with E^2 factored out, indices <radial zones:neutrino species> - 1 / MeV^2 / cm
  real*8,allocatable,save :: leak_tau(:,:) !optical depth, , indices <radial zones:neutrino species> - dimensionless
  real*8,allocatable :: kappa_tilde_nu_scat(:,:,:) !opacity without energy dependence for neutrino-matter scattering
                                                    !indices <radial zones:neutrino species:matter species> - 1 / MeV^2 / cm
  real*8,allocatable :: kappa_tilde_nu_abs(:,:,:) !opacity without energy dependence for neutrino absorption
                                                    !indices <radial zones:neutrino species:matter species> - 1 / MeV^2 / cm
  !local EOS variables, so we can use cgs units
  real*8,allocatable,save :: lx1(:) !spatial positions of cell center in cm, indices <radial position>
  real*8,allocatable,save :: lx1i(:) !spatial positions of inner cell interface in cm, indices <radial position>
  real*8,allocatable,save :: lvol(:) !cell volume in cm^3, indices <radial position>
  real*8,allocatable,save :: lIarea(:) !inverse cell area in cm^-2, indices <radial position>
  real*8,allocatable :: lrho(:) !density in g/cm^3, indices <radial position>
  real*8,allocatable :: ltemp(:) !temp in MeV, indices <radial position>
  real*8,allocatable :: lye(:) !number fraction, indices <radial position>
  real*8,allocatable :: lpress(:) !press in dynes/cm^2, indices <radial position>
  real*8,allocatable :: leps(:) !internal specific energy in erg/g, indices <radial position>
  real*8,allocatable :: lent(:) !entropy in k_b / baryon, indices <radial position>
  real*8,allocatable :: lcs2(:) !speed of sound in cm/s, indices <radial position>
  real*8,allocatable :: lmass(:) !cell mass in g, indices <radial position>

  !degeneracy parameters
  real*8,allocatable :: eta_e(:) !e degeneracy, indices <radial position>, including rest mass - dimensionless
  real*8,allocatable :: eta_n(:) !n degeneracy, indices <radial position>, including rest mass, relative to mass of EOS - dimensionless
  real*8,allocatable :: eta_p(:) !p degeneracy, indices <radial position>, including rest mass, relative to mass of EOS - dimensionless
  real*8,allocatable :: eta_hat(:) !n degeneracy - p dengeneracy - Qnp/temp, removes mass difference, indices <radial position> - dimensionless
  real*8,allocatable :: eta_nue(:) !nu_e degeneracy, from eta_e - eta_n + eta_p, indices <radial position> - dimensionless
  real*8,allocatable :: eta_nua(:) !nu_a degeneracy = -eta_nue, indices <radial position> - dimensionless
  real*8,allocatable :: eta_nux(:) !nu_x degeneracy = 0.d00, indices <radial position> - dimensionless
  real*8,allocatable :: lepton_blocking(:,:) !blocking terms for electrons and positrons, indices <radial position,lepton type> - dimensionless
  
  !mass fractions
  real*8,allocatable :: xxn(:) !neutron mass fraction, indices <radial position>
  real*8,allocatable :: xxp(:) !proton mass fraction, indices <radial position>
  real*8,allocatable :: xxa(:) !alpha mass fraction, indices <radial position>
  real*8,allocatable :: xxh(:) !heavy mass fraction, indices <radial position>
  real*8,allocatable :: xabar(:) !abar fraction, indices <radial position>
  real*8,allocatable :: xzbar(:) !zbar fraction, indices <radial position>
  real*8,allocatable :: eta_pn(:) !proton number density corrected for neutron final state blocking, indices <radial position> - dimensionless
  real*8,allocatable :: eta_np(:) !neutron number density corrected for proton final state blocking, indices <radial position> - dimensionless

  !neutrino sphere
  real*8 :: ns_energy(3) !neutrino sphere energies, indices <neutrino species> - in MeV
  integer :: ns_location(3) !neutrino sphere location, indices <neutrino species> - radial index
  logical :: have_ns_energies !flag that if true uses old energies as starting point, otherwise 18.0MeV is used

  !outputs
  real*8,allocatable :: lum(:,:) !luminosity of neutrinos streaming though cell, indices <radial position:neutrino species> in erg/s
  real*8,allocatable :: depsdt(:) !energy change in cell due to neutrinos, indices <radial position> in erg/g/s
  real*8,allocatable :: dyedt(:) !change in electron number fraction, indices <radial position> in number fraction / s
  real*8,allocatable :: ave_enr_nu(:,:) !average neutrino energy, indices <radial position:neutrino species> in MeV
  real*8 :: heat_erms(3) !root mean squared energy at neutrinosphere (F5/F3), indices <neutrino species> in MeV
  real*8 :: heat_em(3) !mean energy at neutrinosphere (F4/F3), indices <neutrino species> in MeV

  !constants
  real*8,parameter :: me_mev = 0.510998910d0 !mass of electron in MeV
  real*8,parameter :: sigma_0 = 1.76d-44 ! in units of cm^2
  real*8,parameter :: alpha = 1.23d0 ! dimensionless
  real*8,parameter :: Qnp = 1.293333d0 !neutron proton mass difference in MeV
  real*8,parameter :: hc_mevcm = 1.23984172d-10 !hc in units of MeV*cm
  real*8,parameter :: rho_min_leak = 1.0d4 !in g/cm^3
  real*8,parameter :: Cv = 0.5d0 + 2.0d0*0.23d0 ! dimensionless
  real*8,parameter :: Ca = 0.5d0 !dimensionless
  real*8,parameter :: gamma_0 = 5.565d-2 ! dimensionless
  real*8,parameter :: fsc = 1.0d0/137.036d0 ! fine structure constant, dimensionless

  real*8 :: ldt !dt in seconds

  integer :: leak_global_imax
  logical :: output
  logical,save :: have_old_tau

contains
!######################################################################
  subroutine initialize_leakage_rosswog

    implicit none

    integer i

    !allocate leakage variables
    allocate(R_eff(n1,3))
    allocate(R_loc(n1,3))
    allocate(R_diff(n1,3))
    allocate(Q_eff(n1,3))
    allocate(Q_loc(n1,3))
    allocate(Q_diff(n1,3))

    allocate(chi(n1,3))
    allocate(zeta(n1,3))
    allocate(leak_tau(n1,3))
    allocate(kappa_tilde_nu_scat(n1,3,3))
    allocate(kappa_tilde_nu_abs(n1,3,3))

    allocate(lx1(n1))
    allocate(lx1i(n1))
    allocate(lvol(n1))
    allocate(lIarea(n1))
    allocate(lrho(n1))
    allocate(ltemp(n1))
    allocate(lye(n1))
    allocate(lpress(n1))
    allocate(leps(n1))
    allocate(lent(n1))
    allocate(lcs2(n1))
    allocate(lmass(n1))

    allocate(eta_e(n1))
    allocate(eta_n(n1))
    allocate(eta_p(n1))
    allocate(eta_hat(n1))
    allocate(eta_nue(n1))
    allocate(eta_nua(n1))
    allocate(eta_nux(n1))
    allocate(lepton_blocking(n1,2))

    allocate(lum(n1,3))
    allocate(depsdt(n1))
    allocate(dyedt(n1))
    allocate(ave_enr_nu(n1,3))

    allocate(xxn(n1))
    allocate(xxp(n1))
    allocate(xxa(n1))
    allocate(xxh(n1))
    allocate(xzbar(n1))
    allocate(xabar(n1))

    allocate(eta_pn(n1))
    allocate(eta_np(n1))

    R_loc(:,:) = 0.0d0
    R_diff(:,:) = 0.0d0
    Q_loc(:,:) = 0.0d0
    Q_diff(:,:) = 0.0d0

    kappa_tilde_nu_scat(:,:,:) = 0.0d0
    kappa_tilde_nu_abs(:,:,:) = 0.0d0

    lx1(:) = 0.0d0
    lx1i(:) = 0.0d0
    lvol(:) = 0.0d0
    lIarea(:) = 0.0d0
    lrho(:) = 0.0d0
    ltemp(:) = 0.0d0
    lye(:) = 0.0d0
    lpress(:) = 0.0d0
    leps(:) = 0.0d0
    lent(:) = 0.0d0
    lcs2(:) = 0.0d0
    lmass(:) = 0.0d0

    eta_e(:) = 0.0d0
    eta_n(:) = 0.0d0
    eta_p(:) = 0.0d0
    eta_hat(:) = 0.0d0
    eta_nue(:) = 0.0d0
    eta_nua(:) = 0.0d0
    eta_nux(:) = 0.0d0
    lepton_blocking(:,:) = 0.0d0

    xxn(:) = 0.0d0
    xxp(:) = 0.0d0
    xxa(:) = 0.0d0
    xxh(:) = 0.0d0
    xzbar(:) = 0.0d0
    xabar(:) = 0.0d0

    eta_pn(:) = 0.0d0
    eta_np(:) = 0.0d0

    output = .true.
    have_old_tau = .false.

    lx1 = x1/length_gf
    lx1i= x1i/length_gf

    do i=ghosts1+1,n1-2
       lvol(i) = 4.0d0*pi/3.0d0*(lx1i(i+1)**3-lx1i(i)**3)
       lIarea(i) = 1.0d0/(4.0d0*pi*lx1(i)**2)
    enddo

  end subroutine initialize_leakage_rosswog
!######################################################################
  subroutine deallocate_leakage_rosswog

    implicit none

    !deallocate leakage variables
    deallocate(R_eff)
    deallocate(R_loc)
    deallocate(R_diff)
    deallocate(Q_eff)
    deallocate(Q_loc)
    deallocate(Q_diff)

    deallocate(chi)
    deallocate(zeta)
    deallocate(leak_tau)
    deallocate(kappa_tilde_nu_scat)
    deallocate(kappa_tilde_nu_abs)

    deallocate(lx1)
    deallocate(lx1i)
    deallocate(lvol)
    deallocate(lIarea)
    deallocate(lrho)
    deallocate(ltemp)
    deallocate(lye)
    deallocate(lpress)
    deallocate(leps)
    deallocate(lent)
    deallocate(lmass)
    deallocate(lcs2)

    deallocate(eta_e)
    deallocate(eta_n)
    deallocate(eta_p)
    deallocate(eta_hat)
    deallocate(eta_nue)
    deallocate(eta_nua)
    deallocate(eta_nux)
    deallocate(lepton_blocking)

    deallocate(lum)
    deallocate(depsdt)
    deallocate(dyedt)
    deallocate(ave_enr_nu)

    deallocate(xxn)
    deallocate(xxp)
    deallocate(xxa)
    deallocate(xxh)
    deallocate(xzbar)
    deallocate(xabar)

    deallocate(eta_pn)
    deallocate(eta_np)

  end subroutine deallocate_leakage_rosswog
!######################################################################
  subroutine leak_rosswog

    implicit none
    integer i

    ! copy over density, energy, temperature, ye
    ! convert units
    lrho  = rho/rho_gf
    ltemp = temp
    lye   = ye    
    leps  = eps/eps_gf
    lent  = entropy
    lcs2 = cs2*clite**2
    lpress = press/press_gf
    ldt = dt/time_gf
    do i=ghosts1+1,n1-2
       lmass(i) = lvol(i)*lrho(i)
    enddo

    i=ghosts1+1
    do while(lrho(i).gt.rho_min_leak.and.i.lt.n1-ghosts1)
       i=i+1
    enddo
    leak_global_imax = i

    !set up degeneracy/mass factors 
    call get_EOS_variables

    !find tau and interpolate degeneracy factors
    call find_ruf_tau

    !find Q_diff and R_diff
    call compute_diffusion_ros

    !find Q_loc and R_loc
    call compute_emission_ros

    !compute the rate of energy and ye loss/gain in each cell for RK
    call compute_update_ros

  end subroutine leak_rosswog
!######################################################################
  subroutine compute_update_ros

    implicit none

    real*8 :: eosdummy(14)
    real*8 :: heat_rad(n1,3),heat_net(n1)
    real*8 :: heat_net_total
    real*8 :: F(2),heat_const
    integer :: keytemp,keyerr
    real*8 get_fermi_integral
    integer :: arraymax

    character(len=100) filename

    integer gain_radius_nue(1)
    integer gain_radius_nua(1)

    integer :: i,j

    lum(:,:) = 0.0d0
    depsdt(:) = 0.0d0
    dyedt(:) = 0.0d0
    R_eff(:,:) = 0.0d0
    Q_eff(:,:) = 0.0d0
    heat_rad(:,:) = 0.0d0
    ave_enr_nu(:,:) = 0.0d0
    heat_net(:) = 0.0d0
    heat_net_total = 0.0d0

    heat_const = (1.0d0+3.0d0*alpha**2) * sigma_0 / (me_mev**2 * massn_cgs) * 0.25d0
    
    do i=leak_global_imax,ghosts1+1,-1
       !determine effective rates
       R_eff(i,:) = R_loc(i,:)/(1.0d0+R_loc(i,:)/R_diff(i,:))
       Q_eff(i,:) = Q_loc(i,:)/(1.0d0+Q_loc(i,:)/Q_diff(i,:))
 
       do j=1,3 
          if (Q_diff(i,j).eq.0.0d0) then
             Q_eff(i,j) = 0.0d0
             ave_enr_nu(i,j) = 0.0d0
          else
             ave_enr_nu(i,j) = Q_eff(i,j)/R_eff(i,j)
          endif

          if (R_diff(i,j).eq.0.0d0) then
             R_eff(i,j) = 0.0d0
             ave_enr_nu(i,j) = 0.0d0
          else
             ave_enr_nu(i,j) = Q_eff(i,j)/R_eff(i,j)
          endif
       enddo

    enddo

    !leakage
    do i=ghosts1+1,leak_global_imax
       if (do_heating.and.i.lt.ishock(1)) then
          !here we sacrifice computational efficiency for clarity...
          lepton_blocking(i,1) = 1.0d0/(1.0d0 + exp(eta_e(i) - &
               get_fermi_integral(5,eta_nue(ns_location(1)))/ &
               get_fermi_integral(4,eta_nue(ns_location(1)))))
          lepton_blocking(i,2) = 1.0d0/(1.0d0 + exp(-eta_e(i) - &
               get_fermi_integral(5,eta_nua(ns_location(2)))/ &
               get_fermi_integral(4,eta_nua(ns_location(2)))))

          F(1:2) = (4.275d0*leak_tau(i,1:2)+1.15d0)*exp(-2.0d0*leak_tau(i,1:2))*lepton_blocking(i,1:2)
          
          if (GR) then
             ! contains GR terms
             ! in erg/s, divide by alp to get luminosity value at r=x1(i), W(i)*alp((1+v) to go to fluid frame
             heat_rad(i,1) = heat_fac*heat_const * lrho(i) * xxn(i) * lum(i-1,1) * &
                  lIarea(i) * heat_erms(1)**2 * F(1) * lvol(i) / (alp(i)**2*W(i)*(1.0d0+v(i)))
             heat_rad(i,2) = heat_fac*heat_const * lrho(i) * xxp(i) * lum(i-1,2) * &
                  lIarea(i) * heat_erms(2)**2 * F(2) * lvol(i) / (alp(i)**2*W(i)*(1.0d0+v(i)))
          
             !don't take more energy then available
             heat_rad(i,:) = min(lum(i-1,:)/alp(i)**2,heat_rad(i,:))
          else
             !newtonian
             ! in erg/s, divide by to get luminosity value at r=x1(i)
             heat_rad(i,1) = heat_fac*heat_const * lrho(i) * xxn(i) * lum(i-1,1) * &
                  lIarea(i) * heat_erms(1)**2 * F(1) * lvol(i)
             heat_rad(i,2) = heat_fac*heat_const * lrho(i) * xxp(i) * lum(i-1,2) * &
                  lIarea(i) * heat_erms(2)**2 * F(2) * lvol(i)
          
             !don't take more energy then available
             heat_rad(i,:) = min(lum(i-1,:),heat_rad(i,:))
          endif

          !in MeV/cm^3/s (same units as Q)
          heat_rad(i,1) = heat_rad(i,1) / lvol(i) * erg_to_mev
          heat_rad(i,2) = heat_rad(i,2) / lvol(i) * erg_to_mev

          if (GR) then
             !red shifted value at infinity
             lum(i,:) = lum(i-1,:) + (Q_eff(i,:)-heat_rad(i,:))* &
                  mev_to_erg*lvol(i)*alp(i)**2*W(i)*(1.0d0+v(i))*X(i)
          else
             lum(i,:) = lum(i-1,:) + (Q_eff(i,:)-heat_rad(i,:))*mev_to_erg*lvol(i)
          endif

          depsdt(i) = (sum(-Q_eff(i,:)+heat_rad(i,:)))/lrho(i)*mev_to_erg
          dyedt(i) = -(R_eff(i,1)-R_eff(i,2)-heat_rad(i,1)/heat_em(1)+heat_rad(i,2)/heat_em(2))* &
               massn_cgs/lrho(i)
          
          heat_net(i) =  (-Q_eff(i,1)-Q_eff(i,2)+heat_rad(i,1)+heat_rad(i,2))*mev_to_erg*lvol(i)

          heat_net_total = heat_net_total + max(heat_net(i),0.0d0)

          !in ergs/g/s (for output)
          heat_rad(i,1) = heat_rad(i,1) / lrho(i) * mev_to_erg
          heat_rad(i,2) = heat_rad(i,2) / lrho(i) * mev_to_erg

       else
          !just cooling
          if (i.lt.ishock(1)) then
             if (GR) then
                lum(i,:) = lum(i-1,:) + Q_eff(i,:)* &
                     mev_to_erg*lvol(i)*X(i)*alp(i)**2*W(i)*(1.0d0+v(i))
             else
                lum(i,:) = lum(i-1,:) + Q_eff(i,:)*mev_to_erg*lvol(i)
             endif
             depsdt(i) = -sum(Q_eff(i,1:3))*mev_to_erg/lrho(i)
             dyedt(i) = -(R_eff(i,1)-R_eff(i,2))*massn_cgs/lrho(i)
          else
             lum(i,:) = lum(i-1,:)
          endif
       endif


       if (lye(i).le.eos_yemin*1.01d0) then
          !need to surpress any cooling/heating if dyedt.lt.0 near boundary
          if (dyedt(i).lt.0.0d0) then
             dyedt(i) = 0.0d0
             depsdt(i) = 0.0d0
          endif
       endif

      if (lye(i).ge.eos_yemax*0.99d0) then
          !need to surpress any cooling/heating if dyedt.lt.0 near boundary
          if (dyedt(i).gt.0.0d0) then
             dyedt(i) = 0.0d0
             depsdt(i) = 0.0d0
          endif
       endif

       !don't actually use these, done in RK, just a check on the energy
       lye(i) = lye(i) + dyedt(i)*ldt
       leps(i) = leps(i) + depsdt(i)*ldt
       if (leps(i).lt.-energy_shift) then
          write(*,*) "You want an energy less then the table...", &
               i,leps(i),depsdt(i)*ldt,leps(i)-depsdt(i)*ldt,Q_eff(i,:),heat_rad(i,:),&
               "this means there's definetly an error now, maybe already be one though...", &
               xxp(i),lum(i-1,2),lx1(i),F(2),ns_location(2),R_eff(i,:),&
               R_loc(i,:),R_diff(i,:)
          stop
       endif
       leps(i) = leps(i) - depsdt(i)*ldt

       !set variables for RK, get units right
       if (GR) then
          coolingsource(i,2) = W(i)*alp(i)*v(i)*rho(i)*depsdt(i)*eps_gf/time_gf
          coolingsource(i,3) = W(i)*alp(i)*rho(i)*depsdt(i)*eps_gf/time_gf
          coolingsource(i,4) = alp(i)*X(i)*rho(i)*dyedt(i)/time_gf
       else
          coolingsource(i,2) = v1(i)*rho(i)*depsdt(i)*eps_gf/time_gf*sqrt_gamma(i)
          coolingsource(i,3) = rho(i)*depsdt(i)*eps_gf/time_gf*sqrt_gamma(i)
          coolingsource(i,4) = rho(i)*dyedt(i)/time_gf*sqrt_gamma(i)
       endif
          

!!!!! code used to update EOS in operator split format,  i.e. NOT done in RK !!!!!

!       keytemp = 0 !keep eps constant
!       keyerr = 0       
!       call nuc_eos_full(lrho(i),ltemp(i),lye(i),leps(i),lpress(i),lent(i), &
!            lcs2(i),eosdummy(2),eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
!            eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10),eosdummy(11), &
!            eosdummy(12),eosdummy(13),eosdummy(14),keytemp,keyerr,eos_rf_prec)
!       if(keyerr.ne.0) then
!          write(6,*) "############################################"
!          write(6,*) "EOS PROBLEM in leakage:"
!          write(6,*) "timestep number: ",nt
!          write(6,"(i4,1P10E15.6)") i,lx1(i),lrho(i),ltemp(i),leps(i),lye(i)
!          stop "This is bad!"
!       endif
!       if(GR) then
!          lcs2(i) = lcs2(i)/(1.0d0+lpress(i)/lrho(i)/clite**2+leps(i)/clite**2)
!       endif

!!!!! done EOS update !!!!!
       
    enddo
    
    arraymax = min(max(500,ishock(1)+10),n1)
    output = .not.output !alternates output so RK doesn't lead to double output with RK=2
    
    if (output) then
       if (mod(nt,2000).eq.0) then

          if (.not.small_output) then
             filename = trim(adjustl(outdir))//"/blocking.xg"
             open(667,file=filename,status='unknown',position='append')
             write(667,*) '"Time = ',time
             do i=ghosts1+1,arraymax
                write(667,"(1P20E18.9)") lx1(i),lepton_blocking(i,1), &
                     lepton_blocking(i,2)
             enddo
             write(667,*) " "
             write(667,*) " "
             close(667)
          endif
          
          filename = trim(adjustl(outdir))//"/lum.xg"
          open(667,file=filename,status='unknown',position='append')
          write(667,*) '"Time = ',time
          do i=ghosts1+1,arraymax
             write(667,"(1P20E18.9)") lx1(i),lum(i,1),lum(i,2),lum(i,3)
          enddo
          write(667,*) " "
          write(667,*) " "
          close(667)

          filename = trim(adjustl(outdir))//"/tau.xg"
          open(667,file=filename,status='unknown',position='append')
          write(667,*) '"Time = ',time
          do i=ghosts1+1,arraymax
             write(667,"(1P20E18.9)") lx1(i),leak_tau(i,1), &
                  leak_tau(i,2),leak_tau(i,3)
          enddo
          write(667,*) " "
          write(667,*) " "
          close(667)
          
          if (do_heating) then
             
             filename = trim(adjustl(outdir))//"/heat_rad.xg"
             open(667,file=filename,status='unknown',position='append')
             write(667,*) '"Time = ',time
             do i=ghosts1+1,arraymax
                write(667,"(1P20E18.9)") lx1(i),heat_rad(i,1),heat_rad(i,2)
             enddo
             write(667,*) " "
             write(667,*) " "
             close(667)
             
          endif
          
          filename = trim(adjustl(outdir))//"/depsdt.xg"
          open(667,file=filename,status='unknown',position='append')
          write(667,*) '"Time = ',time
          do i=ghosts1+1,arraymax
             write(667,"(1P20E18.9)") lx1(i),depsdt(i)
          enddo
          write(667,*) " "
          write(667,*) " "
          close(667)
         
          if (.not.small_output) then
             filename = trim(adjustl(outdir))//"/dyedt.xg"
             open(667,file=filename,status='unknown',position='append')
             write(667,*) '"Time = ',time
             do i=ghosts1+1,arraymax
                write(667,"(1P20E18.9)") lx1(i),dyedt(i)
             enddo
             write(667,*) " "
             write(667,*) " "
             close(667)
             
             filename = trim(adjustl(outdir))//"/leak_Q.xg"
             open(667,file=filename,status='unknown',position='append')
             write(667,*) '"Time = ',time
             do i=ghosts1+1,arraymax
                write(667,"(1P20E18.9)") lx1(i),Q_loc(i,1),Q_diff(i,1)
             enddo
             write(667,*) " "
             write(667,*) " "
             close(667)
             
             filename = trim(adjustl(outdir))//"/leak_R.xg"
             open(667,file=filename,status='unknown',position='append')
             write(667,*) '"Time = ',time
             do i=ghosts1+1,arraymax
                write(667,"(1P20E18.9)") lx1(i),R_loc(i,1),R_diff(i,1)
             enddo
             write(667,*) " "
             write(667,*) " "
             close(667)
          endif
          
       endif

       if (mod(nt,100).eq.0) then
          
          filename = trim(adjustl(outdir))//"/lum_nu.dat"
          open(667,file=filename,status='unknown',position='append')
          write(667,"(1P20E18.9)") time-t_bounce,lum(leak_global_imax,1),lum(leak_global_imax,2), &
               lum(leak_global_imax,3)
          close(667) 
          
          filename = trim(adjustl(outdir))//"/nu_spheres.dat"
          open(667,file=filename,status='unknown',position='append')
          write(667,"(1P20E18.9)") time-t_bounce,lx1(ns_location(1)),lx1(ns_location(2)), &
               lx1(ns_location(3))
          close(667)

          if (do_heating) then

             gain_radius_nue = maxloc(lum(ghosts1:ishock(1),1))
             gain_radius_nua = maxloc(lum(ghosts1:ishock(1),2))
             
             filename = trim(adjustl(outdir))//"/heat_eff.dat"
             open(667,file=filename,status='unknown',position='append')
             if (GR) then
                write(667,"(1P20E18.9)") time-t_bounce,heat_net_total/ &
                     (lum(gain_radius_nue(1),1)/alp(gain_radius_nue(1)) + & 
                     lum(gain_radius_nua(1),2)/alp(gain_radius_nua(1)))
             else
                write(667,"(1P20E18.9)") time-t_bounce,heat_net_total/ &
                     (lum(gain_radius_nue(1),1) + lum(gain_radius_nua(1),2))
             endif
             close(667)
             
             filename = trim(adjustl(outdir))//"/heat_erms.dat"
             open(667,file=filename,status='unknown',position='append')
             write(667,"(1P20E18.9)") time-t_bounce,heat_erms(1),heat_erms(2),heat_erms(3)
             close(667)

             filename = trim(adjustl(outdir))//"/heat_em.dat"
             open(667,file=filename,status='unknown',position='append')
             write(667,"(1P20E18.9)") time-t_bounce,heat_em(1),heat_em(2),heat_em(3)
             close(667)
             
          endif

       endif
    endif


  end subroutine compute_update_ros

!######################################################################
  subroutine compute_emission_ros

    implicit none
    
    real*8 :: beta
    real*8 :: pair_const, R_pair, Q_pair
    real*8 :: gamma, gamma_const, R_gamma
    real*8 :: block_factor_e, block_factor_a, block_factor_x
    real*8 :: enr_m, enr_p, enr_tilde_m, enr_tilde_p
    real*8 get_fermi_integral


    integer :: i

    !electron & positron capture
    
    beta = pi*clite*(1.0d0+3.0d0*alpha**2)*sigma_0/(hc_mevcm**3*me_mev**2)
    
    R_loc(:,:) = 0.0d0
    Q_loc(:,:) = 0.0d0
    
    do i=ghosts1+1,leak_global_imax

    !electron & positron capture


       R_loc(i,1) = beta*eta_pn(i)*ltemp(i)**5*get_fermi_integral(4,eta_e(i))
       Q_loc(i,1) = beta*eta_pn(i)*ltemp(i)**6*get_fermi_integral(5,eta_e(i))
       R_loc(i,2) = beta*eta_np(i)*ltemp(i)**5*get_fermi_integral(4,-eta_e(i))
       Q_loc(i,2) = beta*eta_np(i)*ltemp(i)**6*get_fermi_integral(5,-eta_e(i))

       !e-e+ pair processes from Ruffert et al.
       block_factor_e = 1.0d0+exp(eta_nue(i)-0.5d0*( &
            get_fermi_integral(4,eta_e(i))/get_fermi_integral(3,eta_e(i)) + &
            get_fermi_integral(4,-eta_e(i))/get_fermi_integral(3,-eta_e(i)) &
            ))
       block_factor_a = 1.0d0+exp(eta_nua(i)-0.5d0*( &
            get_fermi_integral(4,eta_e(i))/get_fermi_integral(3,eta_e(i)) + &
            get_fermi_integral(4,-eta_e(i))/get_fermi_integral(3,-eta_e(i)) &
            ))
       block_factor_x = 1.0d0+exp(eta_nux(i)-0.5d0*( &
            get_fermi_integral(4,eta_e(i))/get_fermi_integral(3,eta_e(i)) + &
            get_fermi_integral(4,-eta_e(i))/get_fermi_integral(3,-eta_e(i)) &
            ))
       
       enr_m = 8.0d0*pi/hc_mevcm**3*ltemp(i)**4*get_fermi_integral(3,eta_e(i))
       enr_p = 8.0d0*pi/hc_mevcm**3*ltemp(i)**4*get_fermi_integral(3,-eta_e(i))

       enr_tilde_m = 8.0d0*pi/hc_mevcm**3*ltemp(i)**5*get_fermi_integral(4,eta_e(i))
       enr_tilde_p = 8.0d0*pi/hc_mevcm**3*ltemp(i)**5*get_fermi_integral(4,-eta_e(i))

       pair_const = sigma_0*clite/me_mev**2*enr_m*enr_p
       
       R_pair =  pair_const/(36.0d0*block_factor_e*block_factor_a)* &
            ((Cv-Ca)**2+(Cv+Ca)**2)

       R_loc(i,1) = R_loc(i,1) + R_pair
       Q_loc(i,1) = Q_loc(i,1) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
       R_loc(i,2) = R_loc(i,2) + R_pair
       Q_loc(i,2) = Q_loc(i,2) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)
       
       R_pair =  pair_const/(9.0d0*block_factor_x**2)*((Cv-Ca)**2+(Cv+Ca-2.0d0)**2)

       R_loc(i,3) = R_loc(i,3) + R_pair
       Q_loc(i,3) = Q_loc(i,3) + R_pair*0.5d0*(enr_tilde_m*enr_p+enr_m*enr_tilde_p)/(enr_m*enr_p)

       !plasmon decay from Ruffert et al.
       gamma = gamma_0*sqrt((pi**2+3.0d0*eta_e(i)**2)/3.0d0)
       block_factor_e = 1.0d0 + exp(eta_nue(i)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
       block_factor_a = 1.0d0 + exp(eta_nua(i)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))
       block_factor_x = 1.0d0 + exp(eta_nux(i)-(1.0d0+0.5d0*gamma**2/(1.0d0+gamma)))

       gamma_const = pi**3*sigma_0*clite*ltemp(i)**8/(me_mev**2*3.0d0*fsc*hc_mevcm**6)* &
            gamma**6*exp(-gamma)*(1.0d0+gamma)
       
       R_gamma = Cv**2*gamma_const/(block_factor_e*block_factor_a)
       R_loc(i,1) = R_loc(i,1) + R_gamma
       Q_loc(i,1) = Q_loc(i,1) + R_gamma*0.5d0*ltemp(i)*(2.0d0+gamma**2/(1.0d0+gamma))
       R_loc(i,2) = R_loc(i,2) + R_gamma
       Q_loc(i,2) = Q_loc(i,2) + R_gamma*0.5d0*ltemp(i)*(2.0d0+gamma**2/(1.0d0+gamma))

       R_gamma = (Cv-1.0d0)**2*4.0d0*gamma_const/block_factor_x**2
       R_loc(i,3) = R_loc(i,3) + R_gamma 
       Q_loc(i,3) = Q_loc(i,3) + R_gamma*0.5d0*ltemp(i)*(2.0d0+gamma**2/(1.0d0+gamma))
       
       !NN Bremsstrahlung (non degenerate limit, BRT06 (with fix to constant out front 1.04 -> 2.0778, Burrows)

       if (do_NNBrem) then
          R_pair = 0.231d0*(2.0778d2*erg_to_mev)*0.5d0*(xxn(i)**2+xxp(i)**2+28.0d0/3.0d0*xxn(i)*xxp(i))* &
               lrho(i)**2*ltemp(i)**(4.5d0)
          Q_pair = R_pair*ltemp(i)/0.231d0*0.504d0
          
          R_loc(i,1) = R_loc(i,1) + R_pair
          Q_loc(i,1) = Q_loc(i,1) + Q_pair
          
          R_loc(i,2) = R_loc(i,2) + R_pair
          Q_loc(i,2) = Q_loc(i,2) + Q_pair
          R_loc(i,3) = R_loc(i,3) + 4.0d0*R_pair
          Q_loc(i,3) = Q_loc(i,3) + 4.0d0*Q_pair
       endif

    enddo

  end subroutine compute_emission_ros
!######################################################################
  subroutine get_EOS_variables
    
    implicit none
    
    integer :: i
    integer :: keytemp, keyerr
    real*8 :: eosdummy(9)

    real*8 :: xmu_e, xmu_n, xmu_p

    do i=ghosts1+1,leak_global_imax

       !call eos to get quantities we need
       keytemp = 1 !keep temperature
       keyerr = 0       
       call nuc_eos_full(lrho(i),ltemp(i),lye(i),eosdummy(1),eosdummy(2),eosdummy(3), &
            eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),xxa(i),xxh(i), &
            xxn(i),xxp(i),xabar(i),xzbar(i),xmu_e,xmu_n,xmu_p,eosdummy(9),keytemp,keyerr,eos_rf_prec)
       if(keyerr.ne.0) then
          write(6,*) "############################################"
          write(6,*) "EOS PROBLEM in leakage:"
          write(6,*) "timestep number: ",nt
          write(6,"(i4,1P10E15.6)") i,lx1(i),lrho(i),ltemp(i),leps(i),lye(i)
          stop "This is bad!"
       endif
       
       ! in our EOS the rest mass difference is in the chemical potentials of the neucleons
       eta_e(i) = xmu_e/ltemp(i)
       eta_p(i) = xmu_p/ltemp(i) 
       eta_n(i) = xmu_n/ltemp(i)

       massfrac_a(i) = xxa(i)
       massfrac_h(i) = xxh(i)
       massfrac_p(i) = xxp(i)
       massfrac_n(i) = xxn(i)
       massfrac_abar(i) = xabar(i)
       massfrac_zbar(i) = xzbar(i)
       
       eta_hat(i) = eta_n(i)-eta_p(i) - Qnp/ltemp(i)
       eta_nue(i) = eta_e(i) - eta_n(i) + eta_p(i) !fully includes effects of rest masses
       eta_nua(i) = -eta_nue(i)
       eta_nux(i) = 0.0d0

    end do

  end subroutine get_EOS_variables
!######################################################################
  subroutine compute_diffusion_ros

    implicit none

    real*8 :: scattering_kappa,abs_kappa,rate_const
    real*8 :: block_factor
    real*8 get_fermi_integral

    integer :: i

    
    zeta(:,:) = 0.0d0
    chi(:,:) = 0.0d0
    
    
    do i=ghosts1+1,leak_global_imax

       !scattering
       scattering_kappa = lrho(i)*avo*0.25d0*sigma_0/me_mev**2
       kappa_tilde_nu_scat(i,1,1) = xxn(i)*scattering_kappa
       kappa_tilde_nu_scat(i,1,2) = xxp(i)*scattering_kappa
       kappa_tilde_nu_scat(i,2,1) = xxn(i)*scattering_kappa
       kappa_tilde_nu_scat(i,2,2) = xxp(i)*scattering_kappa
       kappa_tilde_nu_scat(i,3,1) = xxn(i)*scattering_kappa
       kappa_tilde_nu_scat(i,3,2) = xxp(i)*scattering_kappa
       
       scattering_kappa = lrho(i)*avo*0.0625d0*sigma_0/me_mev**2* &
            xabar(i)*(1.0d0-xzbar(i)/xabar(i))**2 ! only have 1 factor of A because kappa multiples the number fraction, not mass fractions
       kappa_tilde_nu_scat(i,1,3) = xxh(i)*scattering_kappa
       kappa_tilde_nu_scat(i,2,3) = xxh(i)*scattering_kappa
       kappa_tilde_nu_scat(i,3,3) = xxh(i)*scattering_kappa

       eta_pn(i) = avo*lrho(i)*(xxn(i)-xxp(i))/(exp(eta_hat(i))-1.0d0)
       eta_pn(i) = max(eta_pn(i),0.0d0)
       eta_np(i) = avo*lrho(i)*(xxp(i)-xxn(i))/(exp(-eta_hat(i))-1.0d0)
       eta_np(i) = max(eta_np(i),0.0d0)

       if (lrho(i).lt.1.0d11) then
          !non degenerate here, use mass fractions as chemical potentials fail at low densities
          eta_pn(i) = avo*lrho(i)*xxp(i)
          eta_np(i) = avo*lrho(i)*xxn(i)
       endif

       !absorption
!       abs_kappa = lrho(i)*avo*
       abs_kappa = (1.0d0+3.0d0*alpha**2)*0.25d0*sigma_0/me_mev**2
       block_factor = 1.0d0 + exp(eta_e(i)-get_fermi_integral(5,eta_nue(i))/ &
            get_fermi_integral(4,eta_nue(i)))
       kappa_tilde_nu_abs(i,1,1) = eta_np(i)*abs_kappa/block_factor
       kappa_tilde_nu_abs(i,2,1) = 0.0d0 !no absorption of a-type on neutrons
       kappa_tilde_nu_abs(i,3,1) = 0.0d0 !no absorption of x-type neutrinos
       kappa_tilde_nu_abs(i,1,2) = 0.0d0 !no absorption of e-type on protons
       block_factor = 1.0d0 + exp(-eta_e(i)-get_fermi_integral(5,eta_nua(i))/ &
            get_fermi_integral(4,eta_nua(i)))
       kappa_tilde_nu_abs(i,2,2) = eta_pn(i)*abs_kappa/block_factor
       kappa_tilde_nu_abs(i,3,2) = 0.0d0 !no absorption of x-type neutrinos
       kappa_tilde_nu_abs(i,1,3) = 0.0d0 !no absorption on nuclei
       kappa_tilde_nu_abs(i,2,3) = 0.0d0 !no absorption on nuclei
       kappa_tilde_nu_abs(i,3,3) = 0.0d0 !no absorption on nuclei

       !sum up opacities to get zeta (again, factoring out energy dependence)
       zeta(i,1) = kappa_tilde_nu_scat(i,1,1) + kappa_tilde_nu_scat(i,1,2) + &
            kappa_tilde_nu_scat(i,1,3) + kappa_tilde_nu_abs(i,1,1) + &
            kappa_tilde_nu_abs(i,1,2) + kappa_tilde_nu_abs(i,1,3)

       zeta(i,2) = kappa_tilde_nu_scat(i,2,1) + kappa_tilde_nu_scat(i,2,2) + &
            kappa_tilde_nu_scat(i,2,3) + kappa_tilde_nu_abs(i,2,1) + &
            kappa_tilde_nu_abs(i,2,2) + kappa_tilde_nu_abs(i,2,3)

       zeta(i,3) = kappa_tilde_nu_scat(i,3,1) + kappa_tilde_nu_scat(i,3,2) + &
            kappa_tilde_nu_scat(i,3,3) + kappa_tilde_nu_abs(i,3,1) + &
            kappa_tilde_nu_abs(i,3,2) + kappa_tilde_nu_abs(i,3,3)
    enddo

    do i=n1-ghosts1,ghosts1+1,-1
       !integrate zeta to get chi, tau with energy dependence factored out
       chi(i,1) = chi(i+1,1) + zeta(i,1)*(lx1i(i+1)-lx1i(i))
       chi(i,2) = chi(i+1,2) + zeta(i,2)*(lx1i(i+1)-lx1i(i))
       chi(i,3) = chi(i+1,3) + zeta(i,3)*(lx1i(i+1)-lx1i(i))
    enddo
    
    do i=ghosts1+1,leak_global_imax
       !now we can determine diffusion rates
       rate_const = 4.0d0*pi*clite*zeta(i,1)/(hc_mevcm**3*6.0d0*chi(i,1)**2)
       R_diff(i,1) = rate_const*ltemp(i)*get_fermi_integral(0,eta_nue(i))
       Q_diff(i,1) = rate_const*ltemp(i)**2*get_fermi_integral(1,eta_nue(i))
       
       rate_const = 4.0d0*pi*clite*zeta(i,2)/(hc_mevcm**3*6.0d0*chi(i,2)**2)
       R_diff(i,2) = rate_const*ltemp(i)*get_fermi_integral(0,eta_nua(i))
       Q_diff(i,2) = rate_const*ltemp(i)**2*get_fermi_integral(1,eta_nua(i))
       
       rate_const = 16.0d0*pi*clite*zeta(i,3)/(hc_mevcm**3*6.0d0*chi(i,3)**2)
       R_diff(i,3) = rate_const*ltemp(i)*get_fermi_integral(0,eta_nux(i))
       Q_diff(i,3) = rate_const*ltemp(i)**2*get_fermi_integral(1,eta_nux(i))

    enddo

  end subroutine compute_diffusion_ros
!######################################################################
  subroutine find_ruf_tau

    implicit none

        ! local variables:
    ! **** loop control variables
    integer i,jm1
    integer icount
    integer,parameter :: icount_max = 200
    ! **** EOS ****
    real*8 kappa_const_scat_n(n1)
    real*8 kappa_const_scat_p(n1)
!    real*8 kappa_const_scat_h(n1)
    real*8 kappa_const_abs(n1)
    real*8 csn_0,csp_0,t1,t2
    real*8 xerr
    real*8,parameter :: xerr_out = 1.0d-10
    real*8 :: kappa_tot(n1,3)   
    real*8 :: kappa_tot_p(n1,3) 
    real*8 :: kappa_scat_n(n1,3) ! 1/cm
    real*8 :: kappa_scat_p(n1,3) ! 1/cm
!    real*8 :: kappa_scat_h(n1,3) ! 1/cm
    real*8 :: kappa_abs_n(n1)  ! 1/cm
    real*8 :: kappa_abs_p(n1)  ! 1/cm
    real*8 :: local_eta_nux(n1)
    real*8 :: local_eta_nue(n1)
    real*8 :: local_eta_nua(n1)
    real*8 :: eta_nue_eq(n1)

    real*8 dr
    real*8 xlye,xyn,xynp,xyp,xypn
    real*8 get_fermi_integral

    ! initialize some suff:
    kappa_tot(:,:)   = 1.0d0
    kappa_tot_p(:,:) = 1.0d0
    kappa_scat_n(:,:) = 1.0d-5 ! 1/cm
    kappa_scat_p(:,:) = 1.0d-5 ! 1/cm
    kappa_abs_n(:)  = 1.0d-5 ! 1/cm
    kappa_abs_p(:)  = 1.0d-5 ! 1/cm

    local_eta_nux(:) = 0.0d0
    local_eta_nue(:) = 0.0d0
    local_eta_nua(:) = 0.0d0

    ! Neutrino-nucleon scattering transport cross-section
    
    ! C_s,N in equation (A1) of Ruffert et al.
    !
    csn_0 = (1.0d0 + 5.0D0*alpha**2) / 24.0d0
    csp_0 = (4.0d0*(Cv-1.0d0)**2 + 5.0d0*alpha**2) / 24.0d0
    !
    do i=ghosts1+1,leak_global_imax
       ! constant parts of kappa (A6)
       t1 = sigma_0 * avo * lrho(i) * (ltemp(i)/me_mev)**2
       kappa_const_scat_n(i) = csn_0 * t1
       kappa_const_scat_p(i) = csp_0 * t1
       ! (A11) constant part
       kappa_const_abs(i) = (1.0d0+3.0d0*alpha**2)/4.0d0 * t1
    enddo

    ! Loop to get converged result for tau.
    ! This is discussed in the text between equations (A5) and
    ! (A6). Note that for the initial iteration the kappas are set to 1.0d-5
    ! unless we have tau from previous time, then use it as starting point
    icount = 1
    xerr = 1.0d0

    do while(xerr.gt.xerr_out .and. icount.lt.icount_max)
       ! copy over current into previous kappa
       kappa_tot_p = kappa_tot

       ! set up new kappa based on individual
       ! contributions
       ! nu_e; (A17)
       kappa_tot(:,1) = &
              kappa_scat_p(:,1) &
            + kappa_scat_n(:,1) &
            + kappa_abs_n(:)
       ! antis; (A18)
       kappa_tot(:,2) = &
              kappa_scat_p(:,2) &
            + kappa_scat_n(:,2) &
            + kappa_abs_p(:)

       ! nu_xl (A19)
       kappa_tot(:,3) = &
            + kappa_scat_p(:,3) &
            + kappa_scat_n(:,3)

       ! Integrate optical depths: Equation (A20)
       ! Note that this is done for energy transport
       if(icount.gt.2.or..not.have_old_tau) then
          leak_tau(:,:) = 0.0d0
          if (GR) then
             do i=leak_global_imax,ghosts1+1,-1
                dr = lx1i(i+1)-lx1i(i)
                leak_tau(i,1:3) = leak_tau(i+1,1:3) + &
                     X(i)*kappa_tot(i,1:3)*dr 
             enddo
          else
             do i=leak_global_imax,ghosts1+1,-1
                dr = lx1i(i+1)-lx1i(i)
                leak_tau(i,1:3) = leak_tau(i+1,1:3) + &
                     kappa_tot(i,1:3)*dr 
             enddo
          endif
       endif

       jm1=1 !energy optical depth, switch to 0 for number
       do i=ghosts1+1,leak_global_imax
          local_eta_nux(i) = 0.0d0   ! (A2)
          ! (A5) equilibrium eta, we have rest masses in our chemical potentials
          ! no need to include mass difference in eta
          eta_nue_eq(i) = eta_nue(i)
          ! (A3); note that the ^0 etas are set to 0.0d0
          local_eta_nue(i) =  eta_nue_eq(i) * (1.0d0-exp(-leak_tau(i,1))) 
          ! (A4)
          local_eta_nua(i) = -eta_nue_eq(i) * (1.0d0-exp(-leak_tau(i,2))) 

          !assuming completely dissociated, valid in side shock
          xlye = lye(i)
          ! (A8)
          xyn = (1.0d0-xlye) / (1.0d0 + 2.0d0/3.0d0 * max(eta_n(i),0.0d0))
          xyp = xlye / (1.0d0 + 2.0d0/3.0d0*max(eta_p(i),0.0d0))
          t1 = exp(-eta_hat(i))
          ! (A13)
          xynp = max((2.0d0*xlye-1.0d0)/ (t1-1.0d0),0.0d0)
          ! (A14)
          xypn = max(xynp * t1,0.0d0)
 
          ! electron neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nue(i)) / & 
               get_fermi_integral(2+jm1,local_eta_nue(i))
          ! (A15)
          t2 = 1.0d0 + exp(eta_e(i)-get_fermi_integral(5,local_eta_nue(i)) / &
               get_fermi_integral(4,local_eta_nue(i)))
          ! (A6)

          kappa_scat_n(i,1) = kappa_const_scat_n(i) * xyn  * t1
          kappa_scat_p(i,1) = kappa_const_scat_p(i) * xyp  * t1

          ! (A11)
          kappa_abs_n(i) = kappa_const_abs(i) * xynp * t1 / t2 

          ! anti-electron neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nua(i)) / & 
               get_fermi_integral(2+jm1,local_eta_nua(i))
          ! (A16)
          t2 = 1.0d0 + exp(-eta_e(i)-get_fermi_integral(5,local_eta_nua(i)) / &
               get_fermi_integral(4,local_eta_nua(i)))
          ! (A6)
          kappa_scat_n(i,2) = kappa_const_scat_n(i) * xyn  * t1
          kappa_scat_p(i,2) = kappa_const_scat_p(i) * xyp  * t1

          ! (A12)
          kappa_abs_p(i) = kappa_const_abs(i) * xypn * t1 / t2 
          ! nux neutrinos
          t1 = get_fermi_integral(4+jm1,local_eta_nux(i)) / & 
               get_fermi_integral(2+jm1,local_eta_nux(i))
          ! (A6)
          kappa_scat_n(i,3) = kappa_const_scat_n(i) * xyn * t1
          kappa_scat_p(i,3) = kappa_const_scat_p(i) * xyp * t1

       enddo
    
       ! compute relative change xerr
       xerr = 0.0d0
       do i=ghosts1+1,leak_global_imax
          xerr = max(xerr,abs(kappa_tot(i,1)/kappa_tot_p(i,1)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(i,2)/kappa_tot_p(i,2)-1.0d0))
          xerr = max(xerr,abs(kappa_tot(i,3)/kappa_tot_p(i,3)-1.0d0))
       enddo

       icount = icount + 1

    enddo
    
    if(icount.ge.icount_max) then
       write(6,"(i5,1P10E15.6)") icount,xerr,xerr_out
       stop "icount > icount_max in leakage; leak_rosswog.F90"
    endif

    ! Recompute tau based on the most recent kappa_tot
    leak_tau(:,:) = 0.0d0
    if (GR) then
       do i=leak_global_imax,ghosts1+1,-1
          dr = lx1i(i+1)-lx1i(i)
          leak_tau(i,1:3) = leak_tau(i+1,1:3) + &
               X(i)*kappa_tot(i,1:3)*dr
       enddo
    else
       do i=leak_global_imax,ghosts1+1,-1
          dr = lx1i(i+1)-lx1i(i)
          leak_tau(i,1:3) = leak_tau(i+1,1:3) + &
               kappa_tot(i,1:3)*dr
       enddo
    endif
    
    have_old_tau = .true.

    do i=ghosts1+1,leak_global_imax
       !nu-spheres for heating
       if (leak_tau(i,1).gt.0.66666d0) then
          ns_location(1) = i
       endif
       if (leak_tau(i,2).gt.0.66666d0) then
          ns_location(2) = i
       endif
       if (leak_tau(i,3).gt.0.66666d0) then
          ns_location(3) = i
       endif
    enddo

    heat_erms(1) = ltemp(ns_location(1))*sqrt(get_fermi_integral(5,local_eta_nue(ns_location(1)))/ &
         get_fermi_integral(3,local_eta_nue(ns_location(1))))
    heat_erms(2) = ltemp(ns_location(2))*sqrt(get_fermi_integral(5,local_eta_nua(ns_location(2)))/ &
         get_fermi_integral(3,local_eta_nua(ns_location(2))))
    heat_erms(3) = ltemp(ns_location(3))*sqrt(get_fermi_integral(5,0.0d0)/ &
         get_fermi_integral(3,0.0d0))

    heat_em(1) = ltemp(ns_location(1))*get_fermi_integral(5,local_eta_nue(ns_location(1)))/ &
         get_fermi_integral(4,local_eta_nue(ns_location(1)))
    heat_em(2) = ltemp(ns_location(2))*get_fermi_integral(5,local_eta_nua(ns_location(2)))/ &
         get_fermi_integral(4,local_eta_nua(ns_location(2)))
    heat_em(3) = ltemp(ns_location(3))*get_fermi_integral(5,0.0d0)/ &
         get_fermi_integral(4,0.0d0) !not used
    
    !set degeneracy factors to interpolated values
    eta_nue(:) = local_eta_nue(:)
    eta_nua(:) = local_eta_nua(:)
    eta_nux(:) = local_eta_nux(:)
    
    
  end subroutine find_ruf_tau
!######################################################################
end module leakage_rosswog


!######################################################################
function get_fermi_integral(ifermi,eta)
  implicit none
  integer ifermi
  real*8 get_fermi_integral
  real*8 eta
  real*8 fermi_integral_analytical
  
  fermi_integral_analytical = 0.0d0
  
  ! Expressions for Fermi integrals given in Takahashi et al. 1978 
  if (eta.gt.1.D-3) then  
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             (eta**2/2.0D0 + 1.6449d0)/(1.0D0+EXP(-1.6855d0*eta))
     case (2)
        fermi_integral_analytical = &
             (eta**3/3.0D0 + 3.2899d0*eta)/(1.0D0-EXP(-1.8246d0*eta))
     case (3)
        fermi_integral_analytical = & 
             (eta**4/4.0D0 + 4.9348d0*eta**2+11.3644d0) / &
             (1.0D0+EXP(-1.9039d0*eta))        
     case (4)
        fermi_integral_analytical = &
             (eta**5/5.0D0 + 6.5797d0*eta**3+45.4576d0*eta) / &
             (1.0D0-EXP(-1.9484d0*eta))        
     case (5)
        fermi_integral_analytical = &
             (eta**6/6.0D0 + 8.2247d0*eta**4 + 113.6439d0*eta**2 + &
             236.5323d0)/(1.0D0+EXP(-1.9727d0*eta))
     end select
     
  else
     select case (ifermi)
     case (0)
        fermi_integral_analytical = &
             log10(1.0d0+exp(eta))
     case (1)
        fermi_integral_analytical = &
             EXP(eta)/(1.0D0+0.2159d0*EXP(0.8857d0*eta))
     case (2)
        fermi_integral_analytical = & 
             2.0D0*EXP(eta)/(1.0D0+0.1092d0*EXP(0.8908d0*eta))
     case (3)
        fermi_integral_analytical = & 
             6.0D0*EXP(eta)/(1.0D0+0.0559d0*EXP(0.9069d0*eta))
     case (4)
        fermi_integral_analytical = & 
             24.0D0*EXP(eta)/(1.0D0+0.0287d0*EXP(0.9257d0*eta))
     case (5)
        fermi_integral_analytical = &
             120.0D0*EXP(eta) / (1.0D0 + 0.0147d0*EXP(0.9431d0*eta))
     end select
     
  endif
  get_fermi_integral = fermi_integral_analytical
  
  return
end function get_fermi_integral
  
