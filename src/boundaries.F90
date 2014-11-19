!-*-f90-*-
subroutine boundaries(innerflag,outerflag)

  use GR1D_module
  implicit none

  integer innerflag, outerflag
  integer iflag, oflag
  integer i,gi
  real*8 centralvalue
  
  iflag = innerflag
  oflag = outerflag

! inner boundary first; we don't care about inner/outer flag at this point
! always use reflective boundary for geometry=2, flow boundary for geometry=1
        

  if(geometry.eq.2) then
     centralvalue = 0.0d0
     v1m(ghosts1+1) = centralvalue
     if (do_rotation.and..not.GR) then
        vphim(ghosts1+1) = centralvalue
     endif

     gi=ghosts1
     do i=ghosts1,1,-1
        gi=gi+1
        rho(i) = rho(gi)
        rhop(i) = rhom(gi)
        rhom(i) = rhop(gi)
        eps(i) = eps(gi)
        epsp(i) = epsm(gi)
        epsm(i) = epsp(gi)
        mass(i) = mass(gi)
        press(i) = press(gi)
        pressp(i) = pressm(gi)
        pressm(i) = pressp(gi)
        temp(i) = temp(gi)
        entropy(i) = entropy(gi)
        if (fake_neutrinos) then
           nuchem(i) = nuchem(gi)
           if (do_nupress) then
              press_nu(i) = press_nu(gi)
              dnupdr(i) = dnupdr(gi)
           endif
        endif
        v1(i) = -v1(gi)
        v1p(i) = -v1m(gi)
        v1m(i) = -v1p(gi)
        ye(i) = ye(gi)
        yep(i) = yem(gi)
        yem(i) = yep(gi)

        if(do_rotation.and..not.GR) then
           vphi1(i) = -vphi1(gi)
           vphi1p(i) = -vphi1m(gi)
           vphi1m(i) = -vphi1p(gi)
        endif

     enddo
  else
     gi=ghosts1+1
     rhom(gi) = rho(gi)
     epsm(gi) = eps(gi)
     pressm(gi) = press(gi)
     v1m(gi) = v1(gi)
     yem(gi) = ye(gi)
     do i=ghosts1,1,-1
        rho(i) = rho(gi)
        rhop(i) = rho(gi)
        rhom(i) = rho(gi)
        eps(i) = eps(gi)
        epsp(i) = eps(gi)
        epsm(i) = eps(gi)
        mass(i) = mass(gi)
        press(i) = press(gi)
        pressm(i) = press(gi)
        pressp(i) = press(gi)
        temp(i) = temp(gi)
        entropy(i) = entropy(gi)
        if (fake_neutrinos) then
           nuchem(i) = nuchem(gi)
           if (do_nupress) then
              press_nu(i) = press_nu(gi)
              dnupdr(i) = dnupdr(gi)
           endif
        endif
        v1(i) = v1(gi)
        v1m(i) = v1(gi)
        v1p(i) = v1(gi)
        ye(i) = ye(gi)
        yem(i) = ye(gi)
        yep(i) = ye(gi)

     enddo
  endif

! outer boundaries, if stuff is going out, let boundary be outflow, in stuff is coming in, set v=0
  do i=n1-ghosts1+1,n1
     gi=n1-ghosts1
     rhop(gi) = rho(gi)
     rho(i) = rho(gi)
     rhop(i) = rho(gi)
     rhom(i) = rho(gi)
     if(initial_data.eq."OSC") then
        rho(i) = 1.0d0*rho_gf
        rhom(i) = 1.0d0*rho_gf
        rhop(i) = 1.0d0*rho_gf
        rho(gi) = 1.0d0*rho_gf
        rhop(gi) = rho(gi)
     endif
     epsp(gi) = eps(gi)
     eps(i) = eps(gi)
     epsp(i) = eps(gi)
     epsm(i) = eps(gi)
     mass(i) = mass(gi)
     pressp(gi) = press(gi)
     press(i) = press(gi)
     pressm(i) = press(gi)
     pressp(i) = press(gi)
     temp(i) = temp(gi)
     entropy(i) = entropy(gi)
     if (fake_neutrinos) then
        nuchem(i) = nuchem(gi)
        if (do_nupress) then
           press_nu(i) = press_nu(gi)
           dnupdr(i) = dnupdr(gi)
        endif
     endif
     if(v1(gi).gt.0.0d0) then
        v1p(gi) = v1(gi)
        v1(i) = v1(gi)
        v1m(i) = v1(gi)
        v1p(i) = v1(gi)
     else
        v1p(gi) = 0.0d0
        v1(i) = 0.0d0
        v1p(i) = 0.0d0
        v1m(i) = 0.0d0
     endif
     yep(gi) = ye(gi)
     ye(i) = ye(gi)
     yem(i) = ye(gi)
     yep(i) = ye(gi)
     
     if(do_rotation.and..not.GR) then
        vphip(gi) = vphi(gi)
        vphi1(i) = vphi1(gi)
        vphi1m(i) = vphi1(gi)
        vphi1p(i) = vphi1(gi)
     endif

  enddo
  
end subroutine boundaries
