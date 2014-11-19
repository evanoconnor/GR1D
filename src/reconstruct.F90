! -*-f90-*-
subroutine reconstruct
  
  use GR1D_module
  implicit none
  
  call reconstruct_1
  
end subroutine reconstruct

subroutine reconstruct_1
  use GR1D_module
  implicit none
  
  integer i,m
  real*8 discrim
  real*8 temp_cell(n1)
  
  !piecewise constant...
  
  if(reconstruction_method.eq.'pc') then
     call reconstruct_with_pc
     
  else if(reconstruction_method.eq.'tvd') then
     call reconstruct_with_tvd

  else if(reconstruction_method.eq.'ppm') then
     call reconstruct_with_ppm
     
  else
     stop 'reconstruction method not implemented'
     
  endif

  !call the EOS to update the pressure and the speed of sound at the
  !interfaces, this uses keytemp = 0, need to keep temp at the cell
  !center value though

  temp_cell = temp
  call reconstruction_eos_call(rhop,temp,yep,epsp,pressp,cs2p,1)
  temp = temp_cell
  call reconstruction_eos_call(rhom,temp,yem,epsm,pressm,cs2m,0)
  temp = temp_cell

  !get new conservative variables at interfaces
  call prim2con_if

end subroutine reconstruct_1

subroutine reconstruct_with_pc
  
  use GR1D_module
  implicit none
  
  integer i
  real*8 discrim

  do i=1,n1-1
        
     cs2p(i) = cs2(i)
     cs2m(i) = cs2(i)
        
     rhop(i) = rho(i)
     rhom(i) = rho(i)
        
     epsp(i) = eps(i)
     epsm(i) = eps(i)
        
     pressp(i) = press(i)
     pressm(i) = press(i)
        
     yep(i) = ye(i)
     yem(i) = ye(i)
        
     v1p(i) = v1(i)
     v1m(i) = v1(i)

     if (GR) then
        vp(i) = v(i)
        vm(i) = v(i)
     endif

     if(do_rotation) then
        if (GR) then
           vphip(i) = vphi(i)
           vphim(i) = vphi(i)
        else
           vphi1p(i) = omega(i)*x1i(i+1)
           vphi1m(i) = omega(i)*x1i(i)           
        endif
     endif

     if (GR) then
        if (do_rotation) then
           discrim = 1.0d0-vp(i)**2-twothirds*vphip(i)**2
           if (discrim.lt.0.0d0) then
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
     
           discrim = 1.0d0-vm(i)**2-twothirds*vphim(i)**2
           if (discrim.lt.0.0d0) then
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        else
           discrim = 1.0d0-vp(i)**2
           if (discrim.lt.0.0d0) then
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
           
           discrim = 1.0d0-vm(i)**2
           if (discrim.lt.0.0d0) then
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        endif
     endif
  enddo

end subroutine reconstruct_with_pc

subroutine reconstruct_with_tvd

  use GR1D_module
  implicit none
  
  integer i
  real*8 discrim
  real*8 tvdomega(n1),tvdomegap(n1),tvdomegam(n1)
  
  call tvd_reconstruction(n1,ghosts1,rho,rhop,rhom,tvd_limiter)
  call tvd_reconstruction(n1,ghosts1,eps,epsp,epsm,tvd_limiter)
  call tvd_reconstruction(n1,ghosts1,ye,yep,yem,tvd_limiter)
  call tvd_reconstruction(n1,ghosts1,v1,v1p,v1m,tvd_limiter)
  call tvd_reconstruction(n1,ghosts1,v,vp,vm,tvd_limiter)

  if(do_rotation) then
     if (GR) then
        call tvd_reconstruction(n1,ghosts1,vphi,vphip,vphim,tvd_limiter)
     else
        call tvd_reconstruction(n1,ghosts1,vphi1,vphi1p,vphi1m,tvd_limiter)
     endif
  endif

  if (GR) then
     if (do_rotation) then
        do i=1,n1
           discrim = 1.0d0-vp(i)**2-twothirds*vphip(i)**2
           if (discrim.lt.0.0d0) then
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
        enddo
        
        do i=1,n1
           discrim = 1.0d0-vm(i)**2-twothirds*vphim(i)**2
           if (discrim.lt.0.0d0) then
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        enddo
     else
        do i=1,n1
           discrim = 1.0d0-vp(i)**2
           if (discrim.lt.0.0d0) then
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
        enddo
        
        do i=1,n1
           discrim = 1.0d0-vm(i)**2
           if (discrim.lt.0.0d0) then
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        enddo
     endif
  endif

  
end subroutine reconstruct_with_tvd

subroutine reconstruct_with_ppm

  use GR1D_module
  use ppm
  implicit none

  integer i,gi
  real*8 discrim, cv
  real*8 ppmomega(n1),ppmomegap(n1),ppmomegam(n1)

  call ppm_interpolate(rho,rhop,rhom)
  call ppm_steepen

  call ppm_interpolate(v1,v1p,v1m)
  call ppm_interpolate(eps,epsp,epsm)
  call ppm_interpolate(ye,yep,yem)
  call ppm_interpolate(v,vp,vm)

  if(do_rotation) then
     if (GR) then
        do i=1,ghosts1
           vphi(i) = -vphi(2*ghosts1-i+1)
           vphi(n1-i+1) = vphi(n1-ghosts1)
        enddo
        call ppm_interpolate(vphi,vphip,vphim)
     else
        ppmomega = vphi1/x1
        do i=1,ghosts1
           ppmomega(i) = ppmomega(2*ghosts1-i+1)
           ppmomega(n1-i+1) = ppmomega(n1-ghosts1)
        enddo
        call ppm_interpolate(ppmomega,ppmomegap,ppmomegam)
        do i=2,n1-1
           vphi1p(i) = ppmomegap(i)*x1i(i+1)
           vphi1m(i) = ppmomegam(i)*x1i(i)
        enddo
     endif
  endif

  call ppm_flatten

  call ppm_monotonize(rho,rhop,rhom)
  call ppm_monotonize(v1,v1p,v1m)
  call ppm_monotonize(eps,epsp,epsm)
  call ppm_monotonize(ye,yep,yem)
  call ppm_monotonize(v,vp,vm)

  if(do_rotation) then
     if (GR) then
        do i=1,ghosts1
           vphi(i) = -vphi(2*ghosts1-i+1)
           vphi(n1-i+1) = vphi(n1-ghosts1)
        enddo
        call ppm_monotonize(vphi,vphip,vphim)
     else
        ppmomega = vphi1/x1
        do i=1,ghosts1
           ppmomega(i) = ppmomega(2*ghosts1-i+1)
           ppmomega(n1-i+1) = ppmomega(n1-ghosts1)
        enddo
        call ppm_monotonize(ppmomega,ppmomegap,ppmomegam)
        do i=2,n1-1
           vphi1p(i) = ppmomegap(i)*x1i(i+1)
           vphi1m(i) = ppmomegam(i)*x1i(i)
        enddo
     endif
  endif

  cv = 0.0d0
  if (GR) then
     vm(ghosts1+1) = cv
  else
     v1m(ghosts1+1) = cv
  endif

  gi = ghosts1
  do i=ghosts1,1,-1
     gi=gi+1
     if (GR) then
        v(i) = -v(gi)
        vp(i) = -vm(gi)
        vm(i) = -vp(gi)
     else
        v1(i) = -v1(gi)
        v1p(i) = -v1m(gi)
        v1m(i) = -v1p(gi)
     endif
  enddo

  if(do_rotation) then
     if(GR) then
        gi = ghosts1
        do i=ghosts1,1,-1
           gi=gi+1
           vphi(i) = -vphi(gi)
           vphip(i) = -vphim(gi)
           vphim(i) = -vphip(gi)
        enddo
     else
        gi = ghosts1
        do i=ghosts1,1,-1
           gi=gi+1
           vphi1(i) = -vphi1(gi)
           vphi1p(i) = -vphi1m(gi)
           vphi1m(i) = -vphi1p(gi)
        enddo
     endif
  endif

  if (GR) then
     if(do_rotation) then
        do i=1,n1
           discrim = 1.0d0-vp(i)**2-twothirds*vphip(i)**2
           if (discrim.lt.0.0d0) then
              write (6,*) vphip(i), vphi(i), vp(i)
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
        enddo
        
        do i=1,n1
           discrim = 1.0d0-vm(i)**2-twothirds*vphim(i)**2
           if (discrim.lt.0.0d0) then
              write (6,*) vphim(i), vphi(i), vm(i), n1
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        enddo
     else
        do i=1,n1
           discrim = 1.0d0-vp(i)**2
           if (discrim.lt.0.0d0) then
              write (*,*) "Problem in reconstruct", i
              stop
           endif
           Wp(i) = 1.0d0/sqrt(discrim)
        enddo
        
        do i=1,n1
           discrim = 1.0d0-vm(i)**2
           if (discrim.lt.0.0d0) then
              write(*,*) "Problem in reconstruct", i
              stop
           endif
           Wm(i) = 1.0d0/sqrt(discrim)
        enddo
     endif
  endif
  
end subroutine reconstruct_with_ppm

subroutine reconstruction_eos_call(rhoin,tempin,yein,epsin,pressin,cs2in,idir)

  use GR1D_module
  implicit none
  
  integer keytemp,keyerr,eosflag
  integer i,idir,itc,j
  real*8 eosdummy(20)
  real*8 rhoin(n1),tempin(n1),yein(n1),epsin(n1),pressin(n1),cs2in(n1)
  real*8 epsin0
  real*8 rfeps
     
  character(len=256) warnline

  do i=ghosts1,n1-ghosts1+1

       keytemp = 0
       keyerr = 0
       call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
            eosdummy(19), &
            cs2in(i), & 
            eosdummy(2),&
            eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
            eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
            eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
            keytemp,keyerr,eoskey,eos_rf_prec)
       if(keyerr.ne.0) then
          ! -> Issues with the EOS, this can happen around bounce
          !    and is due to very large temperature gradients
          !    in the bouncing inner core in adiabatic collapse
          !    for which the EOS was not really designed. The
          !    problems seen here should not show up for leakage/ye_of_rho
          !    runs.
          write(6,*) "############################################"
          write(6,*) "EOS PROBLEM in reconstruction:"
          write(6,*) "timestep number: ",nt
          if(idir.eq.1) then
             write(6,*) "Reconstruction PLUS"
          else
             write(6,*) "Reconstruction MINUS"
          endif
          write(6,"(i4,1P10E15.6)") i,x1(i),rho(i)/rho_gf,temp(i),eps(i)/eps_gf,ye(i)
          write(6,"(i4,1P10E15.6)") i,x1(i),rhoin(i)/rho_gf,tempin(i),epsin(i)/eps_gf,yein(i)
          write(6,*) "debug:"
          write(6,"(i4,1P10E15.6)") i+1,x1(i+1),rhoin(i+1)/rho_gf,tempin(i+1),epsin(i+1)/eps_gf,yein(i+1)
          write(6,"(i4,1P10E15.6)") i-1,x1(i-1),rhoin(i-1)/rho_gf,tempin(i-1),epsin(i-1)/eps_gf,yein(i-1)
          write(6,*) "keyerr: ",keyerr
          call flush(6)
          write(6,*) "Cell center values around cell in which error occured:"
          write(6,*) "rho,temp,eps,ye"
          do j=i-2,i+2
             write(6,"(i4,1P10E15.6)") j,rho(j)/rho_gf,temp(j),eps(j)/eps_gf,ye(j)
          enddo
          write(6,*) "Reconstructed values:"
          write(6,*) "rho,temp,eps,ye"
          j=i
          write(6,"(i4,1P10E15.6)") j,rhoin(j)/rho_gf,tempin(j),epsin(j)/eps_gf,yein(j)
          call flush(6)
          if(.not.fake_neutrinos.and..not.do_leak_ros) then
          ! let's pump in a tiny bit of energy < 0.01*epsin_orig
             itc = 0
             epsin0 = epsin(j)
             do while(keyerr.ne.0.and.itc.lt.25) 
                itc = itc + 1
                epsin(j) = epsin(j) + epsin0 * 1.0001d0
                call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
                     eosdummy(19), &
                     cs2in(i), & 
                     eosdummy(2),&
                     eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                     eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                     eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                     keytemp,keyerr,eoskey,eos_rf_prec)
             enddo
             write(6,*) itc,keyerr
             write(6,*) "############################################"
             if(keyerr.ne.0) then
                stop "problem in reconstruct: eos, could not be fixed."
             endif
          else if(GR.and.temp(i).gt.93.0d0.and.minval(alp(ghosts1+1:n1)).lt.0.45d0) then
             !try fixing the temp and finding the corresponding eps
             keytemp = 1
             keyerr = 0
             call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
                  eosdummy(19), &
                  cs2in(i), & 
                  eosdummy(2),&
                  eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                  eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                  eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                  keytemp,keyerr,eoskey,eos_rf_prec)
             if (keyerr.ne.0) then
                stop "problem in reconstruct: eos, could not be fixed 2."
             endif
             write(*,*) "Success in EOS"
          else if(GR.and.rhoin(i)/rho_gf.gt.1.0d14.and.minval(alp(ghosts1+1:n1)).lt.0.7d0) then
             itc = 0
             write(*,*) "Reducing accuracy"
             !let's reduce the root finding precision by a factor of 10,
             !but don't go below 1.0d-3
             rfeps = eos_rf_prec
             rfeps = min(eos_rf_prec,1.0d-3)
             do while(keyerr.ne.0.and.itc.lt.10) 
                rfeps = min(rfeps*10.0d0,1.0d-3)
                itc = itc+1
                call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
                     eosdummy(19), &
                     cs2in(i), & 
                     eosdummy(2),&
                     eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                     eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                     eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                     keytemp,keyerr,eoskey,rfeps)
             enddo
             if(keyerr.ne.0) then
                stop "problem in reconstruct: eos, could not be fixed."
             endif
          else if(GR.and.do_leak_ros.and.ye(i+2)-ye(i-2).gt.0.2d0) then
             itc = 0
             write(*,*) "Reducing accuracy"
             !let's reduce the root finding precision by a factor of 10,
             !but don't go below 1.0d-3
             rfeps = eos_rf_prec
             rfeps = min(eos_rf_prec,1.0d-3)
             do while(keyerr.ne.0.and.itc.lt.10) 
                rfeps = min(rfeps*10.0d0,1.0d-3)
                itc = itc+1
                call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
                     eosdummy(19), &
                     cs2in(i), & 
                     eosdummy(2),&
                     eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                     eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                     eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                     keytemp,keyerr,eoskey,rfeps)
             enddo
             if(keyerr.ne.0) then
                stop "problem in reconstruct: eos, could not be fixed."
             endif
          else if(GR.and.do_leak_ros.and..false.) then
             itc = 0
             write(*,*) "Reducing accuracy"
             !let's reduce the root finding precision by a factor of 10,
             !but don't go below 1.0d-3
             rfeps = eos_rf_prec
             rfeps = min(eos_rf_prec,1.0d-3)
             do while(keyerr.ne.0.and.itc.lt.10) 
                rfeps = min(rfeps*10.0d0,1.0d-3)
                itc = itc+1
                call eos_full(i,rhoin(i),tempin(i),yein(i),epsin(i),pressin(i),eosdummy(20), & 
                     eosdummy(19), &
                     cs2in(i), & 
                     eosdummy(2),&
                     eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                     eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                     eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                     keytemp,keyerr,eoskey,rfeps)
             enddo
             if(keyerr.ne.0) then
                stop "problem in reconstruct: eos, could not be fixed."
             endif
          else
             stop "problem in reconstruct: eos"
          endif
       endif
  enddo
  
end subroutine reconstruction_eos_call
