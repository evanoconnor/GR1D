!-*-f90-* 
!purpose is to update q_M1(:,:,:,3), the Eddington tensor, and to
!computer higher order and other second moments to store in
!q_M1_extra(:,:,:,:).  We do this for the zone center and the
!reconstructed values

subroutine M1_closure

  use GR1D_module
  implicit none
  
  !local variables
  real*8 :: tol,err
  real*8 :: oneM1en,oneM1flux,oneM1eddy_guess
  real*8 :: alp2,invalp,invalp2,X2,invX,invX2,W2,v2,oneW,onev,onealp,oneX
  real*8 :: udown(2)
  real*8 :: littlehupdown(2,2),Tupmunu(2,2)
  real*8 :: localJ,H2,Hup(2),Kthick(2,2),Kthin(2,2),Kup(2,2)
  real*8 :: ff2,ff3,ff4,chi,oldprrguess
  real*8 :: Lthin,Lthick,Luprrr,Ldownfupfr
  real*8 :: Wuprrr,Wdownfupfr

  integer :: h,k,i,j,ii,jj
  integer :: count
  real*8 :: justshyofone
  
  tol = 1.0d-8 !tolerance for convergence
  justshyofone = 0.9999999999d0

  if (M1_testcase_number.eq.3) then
     q_M1m(:,:,:,3,1) = justshyofone !P_{rr}/E
     q_M1_extram(:,:,:,1,1) = justshyofone
     q_M1(:,:,:,3) = justshyofone !P_{rr}/E
     q_M1_extra(:,:,:,1) =0.0d0 !P_{phi}^{phi}/E,P_{\theta}^{\theta}/E
     q_M1_extra(:,:,:,4) = justshyofone
     q_M1_extra(:,:,:,2) = justshyofone*q_M1(:,:,:,2)
     q_M1_extra(:,:,:,3) = 0.0d0
     q_M1p(:,:,:,3,1) = justshyofone !P_{rr}/E
     q_M1_extrap(:,:,:,1,1) = justshyofone
     return
  endif

  !$OMP PARALLEL DO PRIVATE(alp2,onealp,invalp,invalp2,X2,oneX,invX,invX2,W2, &
  !$OMP oneW,v2,onev,udown,littlehupdown,h,i,j,oneM1en,oneM1flux,oneM1eddy_guess, &
  !$OMP err,count,Tupmunu,localJ,Hup,Kup,ii,jj,H2,ff2,ff3,ff4,chi,Kthin,Kthick,oldprrguess, &
  !$OMP Lthin,Lthick,Luprrr,Ldownfupfr,Wuprrr,Wdownfupfr)
  do k=ghosts1+1,M1_imaxradii

     !constants to help that only depend on spatial zone
     if (GR) then
        alp2 = alp(k)*alp(k)
        onealp = alp(k)
        invalp = 1.0d0/alp(k)
        invalp2 = 1.0d0/alp2
        X2 = X(k)*X(k)
        oneX = X(k)
        invX = 1.0d0/X(k)
        invX2 = 1.0d0/X2
     else
        alp2 = 1.0d0
        onealp = 1.0d0
        invalp = 1.0d0
        invalp2 = 1.0d0
        X2 = 1.0d0
        oneX = 1.0d0
        invX = 1.0d0
        invX2 = 1.0d0
     endif

     if (v_order.eq.-1) then
        if (GR) then
           W2 = W(k)**2
           oneW = W(k)
           v2 = v(k)**2
           onev = v(k)    
        else
           if (v1(k).gt.0.5d0) stop "very relativistic and your are assume newtonian in tranport"
           W2 = 1.0d0/(1.0d0-v1(k)**2)
           oneW = sqrt(W2)
           v2 = v1(k)**2
           onev = v1(k)
        endif

     else if (v_order.eq.0) then
        W2 = 1.0d0
        oneW = 1.0d0
        v2 = 0.0d0
        onev = 0.0d0
     else
        stop "add in vorder"
     endif

     !h = g+uu
     udown(1) = -oneW*onealp
     udown(2) = oneW*onev*oneX
     littlehupdown(1,1) = -v2*W2 !1.0d0-udown(1)*udown(1)*invalp2
     littlehupdown(2,1) = udown(2)*invX2*udown(1)
     littlehupdown(1,2) = -udown(1)*invalp2*udown(2)
     littlehupdown(2,2) = W2 !1.0d0+udown(2)*udown(2)*invX2

     do h=1,3 !minus,center,plus
        do i=1,number_species_to_evolve
           do j=1,number_groups

              if (h.eq.1) then !minus state
                 oneM1en = q_M1m(k,i,j,1,1)/q_M1m(k,i,j,1,1)
                 oneM1flux = q_M1m(k,i,j,2,1)/q_M1m(k,i,j,1,1)
                 oneM1eddy_guess = q_M1(k,i,j,3)
              else if (h.eq.2) then !middle state
                 oneM1en = q_M1(k,i,j,1)/q_M1(k,i,j,1)
                 oneM1flux = q_M1(k,i,j,2)/q_M1(k,i,j,1)
                 oneM1eddy_guess = q_M1(k,i,j,3)
              else if (h.eq.3) then !plus state
                 oneM1en = q_M1p(k,i,j,1,1)/q_M1p(k,i,j,1,1)
                 oneM1flux = q_M1p(k,i,j,2,1)/q_M1p(k,i,j,1,1)
                 oneM1eddy_guess = q_M1(k,i,j,3)
              endif
              
              err = 1.0d0
          
              !begin iteration, in future do this with a NR
              count = 0
              Tupmunu(1,1) = oneM1en*invalp2
              Tupmunu(1,2) = oneM1flux*invX2*invalp
              Tupmunu(2,1) = Tupmunu(1,2)

              do while (err>tol)
                 count = count + 1
                 Tupmunu(2,2) = oneM1eddy_guess*oneM1en*invX2**2

                 localJ = 0.0d0
                 Hup = 0.0d0
                 Kup = 0.0d0
                 do ii=1,2
                    do jj=1,2
                       localJ = localJ + udown(ii)*udown(jj)*Tupmunu(ii,jj)
                       Hup(1) = Hup(1) - udown(ii)*littlehupdown(1,jj)*Tupmunu(ii,jj)  
                       Hup(2) = Hup(2) - udown(ii)*littlehupdown(2,jj)*Tupmunu(ii,jj)  
                    enddo
                 enddo

                 H2 = 0.0d0
                 do ii=1,2
                    do jj=1,2
                       if (ii.eq.1) then
                          H2 = H2 - alp2*littlehupdown(ii,jj)*Hup(ii)*Hup(jj)
                       else
                          H2 = H2 + X2*littlehupdown(ii,jj)*Hup(ii)*Hup(jj)
                       endif
                    enddo
                 enddo

                 ff2 = ((onealp*onev*oneW*Hup(1))**2 + (oneX*oneW*Hup(2))**2 - &
                      2.0d0*onealp*W2*oneX*onev*Hup(2)*Hup(1))/localJ**2*0.999999999999d0

                 if (ff2.ne.ff2) then
                    if (localJ**2.eq.0.0d0) then
                       ff2 = 0.0d0
                    else
                       write(*,*) oneM1en,oneM1flux,nt,k,i,j
                       stop "ff2 is NaNing"
                    endif
                    
                 endif

                 if (M1closure.eq.'ME') then
                    ff3 = ff2*sqrt(ff2)
                    ff4 = ff2*ff2
                    chi = 1.0d0/3.0d0+(3.0d0*ff2-ff3+3.0d0*ff4)*0.1333333333333333333d0
                 else if (M1closure.eq.'LP') then
                    chi = (3.0d0+4.0d0*ff2)/(5.0d0+2.0d0*sqrt(4.0d0-3.0d0*ff2))
                 else
                    stop "define closure"
                 endif
          
                 !have a(J,H^2) = J/H^2 *(3*chi-1)/2
                 !for Kthin, take chi=1
                 if (Hup(1).eq.0.0d0) then
                    Kthin(2,2) = localJ*invX2
                 else if (Hup(2).eq.0.0d0) then
                    Kthin(2,2) = 0.0d0
                 else if (H2.eq.0.0d0) then
                    Kthin(2,2) = 0.0d0
                 else
                    Kthin(2,2) = Hup(2)*Hup(2)*localJ/H2
                 endif

                 if (Kthin(2,2).ne.Kthin(2,2)) then
                    write(*,*) k,i,j
                    write(*,*) Kthin(2,2), "kthin", oneM1en,oneM1flux,H2,Hup,oneM1eddy_guess
                    stop "Kthin is NaN"
                 endif
                 
                 !for kthick, take chi=1/3
                 Kthick(2,2) = localJ*onethird*littlehupdown(2,2)*invX2

                 Kup(2,2) = (3.0d0*chi-1.0d0)/2.0d0*Kthin(2,2) + 3.0d0*(1.0d0-chi)/2.0d0*Kthick(2,2)

                 Tupmunu(2,2) = ((oneW*onev*oneX)**2*localJ +  &
                      2.0d0*oneW*X2*onev*oneX*Hup(2) + X2*X2*Kup(2,2))*invX2**2
                 
                 oldprrguess = oneM1eddy_guess
                 oneM1eddy_guess = Tupmunu(2,2)*X2*X2/oneM1en

                 if (oneM1eddy_guess.ne.oneM1eddy_guess) then
                    write(*,*) "eddy is Nan",Tupmunu(2,2),localJ,Hup(2),Kup(2,2),Kthin(2,2),Kthick(2,2),chi
                    stop
                 endif

                 err = abs(oneM1eddy_guess-oldprrguess)/oneM1eddy_guess

                 !error checking & fixing
                 if (count.gt.900) then
                    !check is issue is because oneM1flux is -1
                    if (oneM1flux.lt.-0.99d0) then
                       oneM1flux = 0.7d0
                       write(*,*) "closure is failing and flux is -1, reset to 0.7d0"
                    endif
                 endif

                 if (count.gt.1000) then
                    if (err.lt.1.0d-5) then
                       write(*,*) "accepting lower error", err,k,i,j,h,onev,oneM1eddy_guess,Hup(2)
                       err = 1.0d-9
                    else
                       if (chi.gt.0.9d0) then
                          write(*,*) "forcing p_rr to X^2, this should only happen near the shock at bounce",k,chi
                          write(*,*) k,i,j,h,oneM1en,oneM1flux,oneM1flux/oneX,q_M1(k,i,j,3),X2
                          oneM1eddy_guess = X2
                          err = 1.0d-9
                       else
                          if (count.gt.1000) then
                             write(*,*) "Error at",count,":", err,k,i,j,h,chi
                             if (count.gt.1020) then
                                write(*,*) "If this is happening close to BH formation (or late stages when", &
                                     " shock has receeded a lot, consider using higher order explicit flux", &
                                     " calculation"
                                stop "Closure is not converging after 1000 iterations"
                             endif
                          endif
                       endif
                    endif
                 endif

              enddo
              
              if (h.eq.1) then !minus state
                 q_M1m(k,i,j,3,1) = oneM1eddy_guess !P_{rr}/E
                 q_M1_extram(k,i,j,1,1) = chi
              else if (h.eq.2) then !middle state
                 if (Hup(2).eq.0.0d0) then
                    Lthin = 0.0d0 !Cardall 108
                 else
                    Lthin = Hup(2)*invX2 !Cardall 108
                 endif
                 Lthick = 0.6d0*Hup(2)*littlehupdown(2,2)*invX2 !Cardall 108
                 Luprrr = (3.0d0*chi-1.0d0)*0.5d0*Lthin + 3.0d0*(1.0d0-chi)*0.5d0*Lthick
                 Ldownfupfr = 3.0d0*(1.0d0-chi)*0.5d0*0.2d0*Hup(2) !Cardall 108

                 Wuprrr = Luprrr + 3.0d0*oneW*onev*invX*Kup(2,2) +  &
                      3.0d0*W2*v2*invX2*Hup(2) + W2*oneW*v2*onev*invX2*invX*localJ !cardall B18
                 Wdownfupfr = oneW*onev*invX*3.0d0*(1.0d0-chi)*0.5d0*localJ*onethird + & 
                      Ldownfupfr !cardall B18

                 q_M1(k,i,j,3) = oneM1eddy_guess !P_{rr}/E
                 q_M1_extra(k,i,j,1) = 3.0d0*(1.0d0-chi)*0.5d0* &
                      localJ*onethird/oneM1en !P_{phi}^{phi}/E,P_{\theta}^{\theta}/E,
                                              !note oneM1en is 1
                 q_M1_extra(k,i,j,4) = chi
                 q_M1_extra(k,i,j,2) = Wuprrr*q_M1(k,i,j,1)
                 q_M1_extra(k,i,j,3) = Wdownfupfr*q_M1(k,i,j,1)
              else if (h.eq.3) then !plus state
                 q_M1p(k,i,j,3,1) = oneM1eddy_guess !P_{rr}/E
                 q_M1_extrap(k,i,j,1,1) = chi
              endif
           enddo
        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO! end do

end subroutine M1_closure

