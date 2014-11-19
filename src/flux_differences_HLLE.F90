! -*-f90-*-
subroutine flux_differences_hlle
    
  use GR1D_module
  use timers, only: t_fdhlle
  implicit none
  
  integer(kind=4) :: eosflag,keyerr,keytemp
  integer(kind=4) :: i1,i,mw,m
  
  real(kind=8), allocatable :: smin(:),smax(:), dspeed(:), acl(:), acr(:)
  real(kind=8), allocatable :: flux(:,:), fluxl(:,:), fluxr(:,:)
  real(kind=8), allocatable :: ul(:),ur(:)
  real(kind=8), allocatable :: dthx(:)
  real(kind=8), allocatable :: eigenvaluesl(:,:), eigenvaluesr(:,:)
  
  real(kind=8) :: rp,rm,dx,dx_mod
  
  !timers
  real*8 :: t1 = 0.0d0
  real*8 :: t2 = 0.0d0
  
  call cpu_time(t1)

  allocate(smin(n1))
  allocate(smax(n1))
  allocate(dspeed(n1))
  allocate(flux(n1,n_cons))
  allocate(fluxl(n1,n_cons))
  allocate(fluxr(n1,n_cons))
  allocate(ul(n1))
  allocate(ur(n1))
  allocate(acl(n1))
  allocate(acr(n1))
  allocate(dthx(n1))
  allocate(eigenvaluesl(n1,n_cons))
  allocate(eigenvaluesr(n1,n_cons))		
  
  ! Set left and right sound speeds, velocities and 
  ! eigenvalues (characteristic speeds) in Romero et al. (1996).
  ! In the case of rotation we make the approximation that
  ! vphi << v_r and disregard the contribution of vphi to the
  ! eigenvalues. This may lead to inconsistencies in very rapidly
  ! spinning cases.
  do i=1,n1-1
     acl(i) = sqrt(cs2p(i))
     acr(i) = sqrt(cs2m(i+1))
     ul(i)  = v1p(i)
     ur(i)  = v1m(i+1)
     if(GR) then
        ul(i) = vp(i)
        ur(i) = vm(i+1)
        eigenvaluesl(i,1) = ul(i)
        eigenvaluesl(i,2) = (ul(i)+acl(i))/(1.0d0+ul(i)*acl(i))
        eigenvaluesl(i,3) = (ul(i)-acl(i))/(1.0d0-ul(i)*acl(i))
        eigenvaluesr(i,1) = ur(i)
        eigenvaluesr(i,2) = (ur(i)+acr(i))/(1.0d0+ur(i)*acr(i))
        eigenvaluesr(i,3) = (ur(i)-acr(i))/(1.0d0-ur(i)*acr(i))
     endif
  enddo

  ! set wave speeds from eigenvalues
  ! only two wave speeds since we are using HLLE, need max and min
  if (GR) then
     do i=1,n1-1
        smin(i) = dmin1(eigenvaluesl(i,3), &
             eigenvaluesr(i,3), &
             eigenvaluesl(i,2), &
             eigenvaluesr(i,2), &
             eigenvaluesl(i,1), &
             eigenvaluesr(i,1), 0.0d0) 
        smax(i) = dmax1(eigenvaluesl(i,3), &
             eigenvaluesr(i,3), &
             eigenvaluesl(i,2), &
             eigenvaluesr(i,2), &
             eigenvaluesl(i,1), &
             eigenvaluesr(i,1), 0.0d0)
     enddo
  else 
     smin(1:n1-1) = dmin1(ul(1:n1-1) - acl(1:n1-1), &
          ur(1:n1-1) - acr(1:n1-1),0.0d0 )
     smax(1:n1-1) = dmax1(ul(1:n1-1) + acl(1:n1-1), &
          ur(1:n1-1) + acr(1:n1-1),0.0d0 )
  endif
  
  ! determine L,R fluxes at interfaces from conserved variables at interfaces			 
  do i=1,n1-1
     if(GR) then
        fluxl(i,1) = qp(i,1)*vp(i)
        fluxl(i,2) = qp(i,2)*vp(i) + pressp(i)
        fluxl(i,3) = (qp(i,3)+pressp(i))*vp(i)
        fluxl(i,4) = qp(i,4)*vp(i)
        
        fluxr(i,1) = qm(i+1,1)*vm(i+1)
        fluxr(i,2) = qm(i+1,2)*vm(i+1) + pressm(i+1)
        fluxr(i,3) = (qm(i+1,3)+pressm(i+1))*vm(i+1)
        fluxr(i,4) = qm(i+1,4)*vm(i+1)
        
        if(do_rotation) then
           fluxl(i,5) = qp(i,5)*vp(i)
           fluxr(i,5) = qm(i+1,5)*vm(i+1)
        endif
        
     else
        fluxl(i,1) = qp(i,1)*v1p(i)
        fluxl(i,2) = qp(i,2)*v1p(i) + pressp(i)
        fluxl(i,3) = (qp(i,3)+pressp(i))*v1p(i)
        fluxl(i,4) = qp(i,4)*v1p(i)
        
        fluxr(i,1) = qm(i+1,1)*v1m(i+1)
        fluxr(i,2) = qm(i+1,2)*v1m(i+1) + pressm(i+1)
        fluxr(i,3) = (qm(i+1,3)+pressm(i+1))*v1m(i+1)
        fluxr(i,4) = qm(i+1,4)*v1m(i+1)
        
        if(do_rotation) then
           fluxl(i,5) = qp(i,5)*v1p(i)
           fluxr(i,5) = qm(i+1,5)*v1m(i+1)
        endif
        
     endif
  enddo

  do i=1,n1
     dspeed(i) = smax(i)-smin(i)
     flux_diff(i,1:n_cons) = 0.0d0
  enddo

  do i=ghosts1,n1-ghosts1+1
     do m=1,n_cons
        ! this is the flux at i + 1/2
        flux(i,m) = ( smax(i)*fluxl(i,m) &
             - smin(i)*fluxr(i,m) &
             + smax(i) * smin(i) &
             * ( qm(i+1,m) - qp(i,m) ) ) &
             / ( dspeed(i) )
     enddo
  enddo

  ! update by flux differencing
  ! Note that we are not confusing indices here
  ! flux(i) is the flux at the interface located at  i + 1/2.
  ! According to Einfeldt, the flux difference needed for the update
  ! is g_i+1/2 - g_i-1/2 ^= g_i+1/2 - g_i-1+1/2
  if(geometry.eq.2) then
     !spherical geometry
     do i=ghosts1+1,n1-ghosts1
        rm = x1i(i)
        rp = x1i(i+1) 
        dx = rp - rm
        dx_mod = (rp**3-rm**3)/(3.0d0*x1(i)**2)
        do m=1,n_cons
           !approx divergence in spherical co-ordinate
           if (GR) then 
              if (m.eq.5) then
                 !since the quanitity goes as r^2 use improved diverence for small r
                 flux_diff(i,m)  =  ( flux(i,m)*alpp(i)/Xp(i)*rp**2 &
                      - flux(i-1,m)*alpm(i)/Xm(i)*rm**2 ) &
                      / x1(i)**2 / dx_mod
              else
                 flux_diff(i,m)  =  ( flux(i,m)*alpp(i)/Xp(i)*rp**2 &
                      - flux(i-1,m)*alpm(i)/Xm(i)*rm**2 ) &
                      / x1(i)**2 / dx
              endif
           else
              if (m.eq.5) then
                 !since the quanitity goes as r^2 use improved diverence for small r
                 flux_diff(i,m) = ( flux(i,m)*rp**2 - flux(i-1,m)*rm**2 ) / dx_mod
              else
                 flux_diff(i,m) = ( flux(i,m)*rp**2 - flux(i-1,m)*rm**2 ) / dx
              endif
           endif
        enddo
     enddo
  else
     !planar geometry
     do i=ghosts1+1,n1-ghosts1
        rm = x1i(i)
        rp = x1i(i+1) 
        dx = rp - rm
        do m=1,n_cons
           flux_diff(i,m) = (flux(i,m) - flux(i-1,m)) / dx
        enddo
     enddo
  endif
  
  deallocate(smin)
  deallocate(smax)
  deallocate(dspeed)
  deallocate(flux)
  deallocate(fluxl)
  deallocate(fluxr)
  deallocate(ul)
  deallocate(ur)
  deallocate(acl)
  deallocate(acr)
  deallocate(dthx)
  deallocate(eigenvaluesl)
  deallocate(eigenvaluesr)
  
  call cpu_time(t2)

  t_fdhlle = t_fdhlle + (t2-t1)
  
end subroutine flux_differences_hlle
