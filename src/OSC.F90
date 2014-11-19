!-*-f90-*-
subroutine OSC

  use GR1D_module
  use poly_eos_module
  implicit none

  real*8 rchange
  real*8 OSCmass, OSCdensity
  real*8 testmass
  logical done

  real*8 xmin,xmax,mindx

  integer i

  if (GR.eqv..false.) then
     stop "Why are you doing OSC with no GR?"
  endif

  dt = 5.0d-11*time_gf

  OSCmass = 1.0d0
  rchange = 10.0d0*OSCmass
  OSCdensity = OSCmass/(4.0d0/3.0d0*pi*rchange**3)

  eoskey = 2
  polygamma = 5.0d0/3.0d0 !should be independant of this value
  polyK =  1.0d-20 !a very small number

  xmin = 0.0d0*length_gf
  xmax = rchange*2.0d0
  mindx = 0.0d0
     
  ! minimum grid spacing; does only apply if geometry = 2
  if (gridtype.ne.'unigrid') stop "work this gridtype into OSC"

  !setup grid
  call grid(xmin,xmax,mindx)

  !fill up zones
  testmass = 0.0d0
  done = .false.
  do i=ghosts1,n1
     if (done) then
        rho(i) = 1.0d0*rho_gf   
     else
        if ((testmass + OSCdensity*volume(i)).lt.1.0d0) then
           rho(i) = OSCdensity
        else
           rho(i) = (1.0d0-testmass)/volume(i)
           done = .true.
        endif
     endif
     
     press(i) = polyK*rho(i)**polygamma
     eps(i) = press(i)/rho(i)/(polygamma-1.d0)
     testmass = testmass + rho(i)*volume(i)
     v(i) = 0.0d0
     v1(i) = 0.0d0
  enddo

  !test
  write(*,*) "OSC: This number should be 1: ", testmass

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
        write(*,*) "OSC with no gravity!?!"
     endif
  endif

  if (GR) then
     write(*,*) "OSC: Using GR, totalmass_grav: ", totalmass/mass_gf
     write(*,*) "               totalmass_bary: ", mass(n1-ghosts1+1)/mass_gf
  else
     write(*,*) "OSC: Using Newtonian, totalmass_bary: ", mass(n1-ghosts1+1)/mass_gf
  endif
  

end subroutine OSC
