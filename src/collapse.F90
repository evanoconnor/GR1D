!-*-f90-*-
subroutine collapse

  use GR1D_module
  use ye_of_rho
  implicit none
  
  real*8 xmin,xmax,mindx
  real*8 discrim
  integer i
  
  ! initial timestep
  dt = 1.0d-7*time_gf

  if(fake_neutrinos) then
     write(*,*) "Initializing Fake Neutrinos!"
     if (do_nupress) then
        write(*,*) "Neutrino Pressure on"
        if (.not.do_ye_of_rho) then
           stop "Neutrino Press with No Ye_of_rho is a bad idea"
        endif
     endif
     if (do_ye_of_rho.and..not.do_yeofrhofit) then
        write(*,*) "Using ye profile"
        call read_yeprofile
     else if (do_ye_of_rho.and.do_yeofrhofit) then
        write(*,*) "Using ye(rho) fit formula"
     else
        write(*,*) "No ye(rho) prescription"
     endif
  endif
  
  write(6,"(A15,A60)") "Initial data: ",profile_name
  
  if (do_profile_rmax) then
     call map_limit(profile_name)
  endif
  
  xmin = 0.0d0*length_gf
  xmax = grid_rmax*length_gf
  
  ! minimum grid spacing; does only apply if geometry = 2 & grid = log
  mindx = grid_custom_dx1*length_gf !1 km
  
  write(*,*) "Setting up grid: ", trim(adjustl(gridtype))
  call grid(xmin,xmax,mindx)
  
  if(profile_type.eq.1) then
     call map_profile(profile_name)
  else
     stop "profile_type not recognized!"
  endif
  
  !find initial mass/gravitational mass
  call mass_interior
  call boundaries(0,0)
  if (GR) then
     if (gravity_active) then
        !GR Collapse & OSC
        mgrav(:) = mass(:)
        call GR_terms_initial
        call GR_boundaries
        !redetermine gravitational and baryonic mass
        call mass_interior
     else
        !planer & sedov
        X(:) = 1.0d0
        Xp(:) = 1.0d0
        Xm(:) = 1.0d0
        alp(:) = 1.0d0
        alpp(:) = 1.0d0
        alpm(:) = 1.0d0

        v(:) = v1(:)
        do i=1,n1
           discrim = 1.0d0 - v(i)**2
           if (discrim.lt.0.0d0) then
              stop "We have a problem, v>1"
           endif
           W(i) = 1.0d0/sqrt(discrim)
        enddo
     endif
  endif

  if (GR) then
     write(*,*) "Using GR, totalmass_grav: ", totalmass/mass_gf, " grams"
     write(*,*) "          totalmass_bary: ", mass(n1-ghosts1+1)/mass_gf, " grams"
  else
     write(*,*) "Using Newtonian, totalmass_bary: ", mass(n1-ghosts1+1)/mass_gf, " grams"
  endif

end subroutine collapse
