!-*-f90-*-
subroutine gravity

  use GR1D_module
  implicit none

  integer i

  if (GR) then
    do i=1,n1
      gravsource(i,1) = 0.0d0
      gravsource(i,2) =   (q(i,2)*v(i)-q(i,3)-q(i,1))*(8.0d0*alp(i)*X(i)* &
           pi*x1(i)*press(i) + alp(i)*X(i)*mgrav(i)/x1(i)**2) &
           + alp(i)*X(i)*press(i)*mgrav(i)/x1(i)**2
      gravsource(i,3) = 0.0d0
      gravsource(i,4) = 0.0d0
    enddo
  else
    gravsource(1:n1,1) = 0.0d0
    gravsource(1:n1,2) = sqrt_gamma(1:n1)* &
         (-rho(1:n1))*dphidr(1:n1)  ! Force
    gravsource(1:n1,3) = sqrt_gamma(1:n1)* &
         (- rho(1:n1))*v1(1:n1)*dphidr(1:n1)
    gravsource(1:n1,4) = 0.0d0
  endif

end subroutine gravity

subroutine con2GR

  use GR1D_module
  implicit none

  integer i
  real*8 mgrav1(n1)
  real*8 nuE(n1)
  integer iminb,imaxb
  real*8 discrim

  iminb = ghosts1+1
  imaxb = n1-1

  !neutrino contribution to energy
  if (do_M1) then
     nuE(:) = energy_nu(:) !it matters how energy_nu is defined, we
                          !defines ours in the lab frame for M1,
                          !radiation tensor is then 3.29 of Shibata et
                          !al. instead of 3.5
  else
     nuE(:) = (W(:)**2-1.0d0)*press_nu(:) + W(:)**2*energy_nu(:)
  endif

  !find new GR quantities
  mgrav1(iminb) = (q(iminb,3)+q(iminb,1)+nuE(iminb))*volume(iminb)
  mgrav(iminb) = mgrav1(iminb) - 4.0d0/3.0d0*pi*(q(iminb,3) &
       +q(iminb,1)+nuE(iminb)) * ( x1i(iminb+1)**3 - x1(iminb)**3 ) 
  mgravi(iminb) = 0.0d0
  discrim = 1.0d0 - 2.0d0*mgrav(iminb)/x1(iminb)
  if (discrim.lt.0.0d0) then
     stop "We have a black hole"
  endif
  X(iminb) = 1.0d0/sqrt(discrim)
  do i=iminb+1,imaxb
     mgrav1(i) = (q(i,3)+q(i,1)+nuE(i))*volume(i)
     mgrav(i) = mgrav(i-1) + mgrav1(i) - 4.0d0/3.0d0*pi*(q(i,3)+q(i,1)+ &
          nuE(i)) *( x1i(i+1)**3 - x1(i)**3 ) + &
          4.0d0/3.0d0*pi*(q(i-1,3)+q(i-1,1)+nuE(i-1)) &
          * ( x1i(i)**3 - x1(i-1)**3 )    
     mgravi(i) = mgravi(i-1) + (q(i-1,3)+q(i-1,1)+nuE(i-1))*volume(i-1)
     discrim = 1.0d0 - 2.0d0*mgrav(i)/x1(i)
     if (discrim.lt.0.0d0) then
        stop "We have a black hole"
     endif
     X(i) = 1.0d0/sqrt(discrim)
  enddo
  
  totalmass = mgravi(n1-ghosts1+1)

  discrim = 1.0d0-2.0d0*mgravi(ghosts1+2)/x1i(ghosts1+2)
  if (discrim.lt.0.0d0) then
     write(*,*) "We have a black hole in GR_terms a", &
          mgravi(ghosts1+2),x1i(ghosts1+2)
     stop
  endif
  Xp(ghosts1+1) = 1.0d0/sqrt(discrim)
  Xm(ghosts1+1) = 1.0d0
  do i=ghosts1+2,n1-1
     discrim = 1.0d0-2.0d0*mgravi(i+1)/x1i(i+1)
     if (discrim.lt.0.0d0) then
        write(*,*) "We have a black hole in GR_terms b", mgravi(i+1),x1i(i+1)
        stop
     endif
     Xp(i) = 1.0d0/sqrt(discrim)
     discrim = 1.0d0-2.0d0*mgravi(i)/x1i(i)
     if (discrim.lt.0.0d0) then
        write(*,*) "We have a black hole in GR_terms c", mgravi(i),x1i(i)
        stop
     endif
     Xm(i) = 1.0d0/sqrt(discrim)	  
  enddo
  discrim = 1.0d0 - 2.0d0*mgravi(n1)/x1i(n1)
  if (discrim.lt.0.0d0) then
     write(*,*) "We have a black hole in GR_terms d", mgravi(n1),x1i(n1)
     stop
  endif
  Xm(n1) = 1.0d0/sqrt(discrim)	

end subroutine con2GR

subroutine GR_alp

  use GR1D_module
  implicit none

  integer i
  real*8 phi_bound
  real*8 nuE(n1)

  !Neutrino contribution to gravitational potential other than mass, i.e. via pressure
  if (do_M1) then
     nuE(:) = press_nu(:) !it matters how press_nu is defined, we
                          !defines ours in the lab frame for M1,
                          !radiation tensor is then 3.29 of Shibata et
                          !al. instead of 3.5
  else
     if (do_rotation) then
        nuE(:) = (W(:)**2*(v(:)**2+twothirds*vphi(:)**2)+1.0d0)*press_nu(:) + &
             W(:)**2*(v(:)**2+vphi(:)**2)*energy_nu(:)
     else
        nuE(:) = ((W(:)*v(:))**2+1.0d0)*press_nu(:) + (W(:)*v(:))**2*energy_nu(:)
     endif
  endif

  if (do_rotation) then
     do i=ghosts1+1,n1-ghosts1
        dphidr(i) = X(i)**2 * ( &
             mgrav(i)/(x1(i))**2 + &
             4.0d0*pi*x1(i) * ( &
             nuE(i)+ press(i) + &
             rho(i) * &
             (1.0d0 + eps(i) + press(i)/rho(i)) * &
             (v(i)**2+twothirds*vphi(i)**2) * W(i)**2))
     enddo
  else
     do i=ghosts1+1,n1-ghosts1
        dphidr(i) = X(i)**2 * ( &
             mgrav(i)/(x1(i))**2 + &
             4.0d0*pi*x1(i) * ( &
             nuE(i)+ press(i) + &
             rho(i) * &
             (1.0d0 + eps(i) + press(i)/rho(i)) * &
             v(i)**2 * W(i)**2))
     enddo
  endif

  phi(ghosts1+1) = 0.0d0 + x1(ghosts1+1)*(dphidr(ghosts1+1))
  phii(ghosts1+1) = 0.0d0
  do i=ghosts1+2,n1-ghosts1
     phi(i) = phi(i-1) + (x1i(i)-x1(i-1))*dphidr(i-1) &
          + (x1(i)-x1i(i))*dphidr(i)
     phii(i) = phi(i-1) + (x1i(i)-x1(i-1))*dphidr(i-1)
  enddo

  phii(n1-ghosts1+1) = phi(n1-ghosts1) + &
       (x1i(n1-ghosts1+1)-x1(n1-ghosts1))*dphidr(n1-ghosts1)

  phi_bound = 0.5d0 * &
       log(1.0d0 - 2.0d0 * mgravi(n1-ghosts1+1)/x1i(n1-ghosts1+1))
  
  do i=ghosts1+1,n1-ghosts1
     phii(i) = phii(i) + (phi_bound-phii(n1-ghosts1+1))
     phi(i) = phi(i) + (phi_bound-phii(n1-ghosts1+1))
     alp(i) = exp(phi(i))
  enddo
  phii(n1-ghosts1+1) = phi_bound

  do i=ghosts1+1,n1-ghosts1
     alpp(i) = exp(phii(i+1))
     alpm(i) = exp(phii(i))
  enddo
  alpm(n1-ghosts1+1) = exp(phii(n1-ghosts1+1))

end subroutine GR_alp

subroutine GR_boundaries

  use GR1D_module
  implicit none

  integer i,gi
  real *8 centralvel

  gi = ghosts1
  centralvel = 0.0d0
  vm(ghosts1+1) = centralvel
  if(do_rotation) then
     vphim(ghosts1+1) = centralvel
  endif
  Wm(ghosts1+1) = 1.0d0
  do i=ghosts1,1,-1
     gi=gi+1
     X(i) = X(gi)
     Xp(i) = Xm(gi)
     Xm(i) = Xp(gi)
     W(i) = W(gi)
     Wp(i) = Wm(gi)
     Wm(i) = Wp(gi)
     v(i) = -v(gi)
     vp(i) = -vm(gi)
     vm(i) = -vp(gi)
     if (do_rotation) then
        vphi(i) = -vphi(gi)
        vphip(i) = -vphim(gi)
        vphim(i) = -vphip(gi)
     endif
     alp(i) = alp(gi)
     alpp(i) = alpm(gi)
     alpm(i) = alpp(gi)
     phi(i) = phi(gi)
     dphidr(i) = dphidr(gi)
     mgrav(i) = mgrav(gi)
  enddo

  !outer boundaries
  gi = n1-ghosts1
  do i=n1-ghosts1+1,n1
     mgrav(i) = mgrav(gi)
     X(i) = Xp(gi)
     Xp(i) = Xp(gi)
     Xm(i) = Xp(gi)
     alp(i) = alpp(gi)
     alpp(i) = alpp(gi)
     alpm(i) = alpp(gi)
     if(v(gi).gt.0.0d0) then
        vp(gi) = v(gi)
        v(i) = v(gi)
        vp(i) = v(gi)
        vm(i) = v(gi)
        Wp(gi) = W(gi)
        W(i) = W(gi)
        Wm(i) = W(gi)
        Wp(i) = W(gi)
     else
        vp(gi) = 0.0d0
        v(i) = 0.0d0
        vp(i) = 0.0d0
        vm(i) = 0.0d0
        Wp(gi) = 1.0d0
        W(i) = 1.0d0
        Wp(i) = 1.0d0
        Wm(i) = 1.0d0
     endif
  enddo

end subroutine GR_boundaries

subroutine GR_terms_initial

  use GR1D_module
  implicit none

  real*8 discrim
  real*8 err, tol
  real*8 Xo(n1), mgrav1(n1), phi_bound
  integer i,gi
  real*8 tau(n1)

  tol = 1.0d-12
  err = 1.0d0
  Xo(:) = 0.0d0

  !iterate to find grav mass from baryonic
  do while(err.gt.tol)
    err = 0.0d0
    do i=ghosts1+1,n1
      discrim = 1.0d0 - 2.0d0*mgrav(i)/x1(i)
      if (discrim.lt.0.0d0) then
         stop "We have a black hole!"
      endif
      X(i) = 1.0d0 / sqrt(discrim)
      if (do_rotation) then
         discrim = 1.0d0 - (X(i)*v1(i))**2 - twothirds*vphi(i)**2
      else
         discrim = 1.0d0 - (X(i)*v1(i))**2
      endif
      if (discrim.lt.0.0d0) then
        stop "We have a problem in GR_terms"
      endif
      W(i) = 1.0d0 / sqrt(discrim)
      ! set up tau + D:
      ! tau + D= rho h W^2 - p 
      ! call tau + D just tau for the time being
  
      tau(i) = rho(i)  * &
         (1.0d0 + eps(i) + press(i)/rho(i)) &
         * W(i)*W(i) &
         - press(i)

    enddo
  
    ! compute gravitational mass
    mgrav1(ghosts1+1) = tau(ghosts1+1)*volume(ghosts1+1)
    mgrav(ghosts1+1) = mgrav1(ghosts1+1) - 4.0d0/3.0d0*pi*tau(ghosts1+1) * &
       ( x1i(ghosts1+2)**3 - x1(ghosts1+1)**3 ) 
    mgravi(ghosts1+1) = 0.0d0
    do i=ghosts1+2,n1 
      mgrav(i) = mgrav(i-1) + 4.0d0/3.0d0*pi*tau(i-1) * &
           (x1i(i)**3 - x1(i-1)**3) &
           + 4.0d0/3.0d0*pi*tau(i) * (x1(i)**3 - x1i(i)**3) 
      mgravi(i) = mgrav(i-1) + &
           4.0d0/3.0d0*pi*tau(i-1) * (x1i(i)**3 - x1(i-1)**3)
      discrim = 1.0d0 - 2.0d0*mgrav(i)/x1(i)
      if (discrim.lt.0.0d0) then
        stop "We have a black hole"
      endif
      Xo(i) = X(i)
      X(i) = 1.0d0/sqrt(discrim)
      if (abs(Xo(i)-X(i)).gt.err) then 
         err = abs(Xo(i)-X(i))
      endif
    enddo
  enddo

  totalmass = mgrav(n1-ghosts1+1)

  do i=ghosts1+1,n1-ghosts1 
     v(i) = X(i)*v1(i)
  enddo
     
  discrim = 1.0d0-2.0d0*mgravi(ghosts1+2)/x1i(ghosts1+2)
  if (discrim.lt.0.0d0) then
    stop "We have a black hole in GR_terms"
  endif
  Xp(ghosts1+1) = 1.0d0/sqrt(discrim)
  Xm(ghosts1+1) = 1.0d0
  do i=ghosts1+2,n1-1
    discrim = 1.0d0-2.0d0*mgravi(i+1)/x1i(i+1)
    if (discrim.lt.0.0d0) then
      stop "We have a black hole in GR_terms"
    endif
    Xp(i) = 1.0d0/sqrt(discrim)
    discrim = 1.0d0-2.0d0*mgravi(i)/x1i(i)
    if (discrim.lt.0.0d0) then
      stop "We have a black hole in GR_terms"
    endif
    Xm(i) = 1.0d0/sqrt(discrim)	  
  enddo
  discrim = 1.0d0 - 2.0d0*mgravi(n1)/x1i(n1)
  if (discrim.lt.0.0d0) then
    stop "We have a black hole in GR_terms"
  endif
  Xm(n1) = 1.0d0/sqrt(discrim)	  

  !Now do lapse  
  call GR_alp

end subroutine GR_terms_initial

