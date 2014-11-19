!-*-f90-*-
! simple 1D PPM. Inspired by CoCoNuT.
! Colella & Woodward, JcP 54, 174 (1984)

module ppm
  
  implicit none

  real*8,allocatable :: delta_r(:)
  real*8,allocatable :: c_1_r(:),c_2_r(:)
  real*8,allocatable :: c_3_r(:),c_4_r(:),c_5_r(:)
    
  !from Dimmelmeier
  real*8 :: eta1 = 20.0d0
  real*8 :: eta2 = 0.05d0
  real*8 :: epsilon1 = 0.01d0
  real*8 :: epsilon2 = 0.33d0
  real*8 :: K0 = 0.1d0 ! does this value make sense?
  real*8 :: omega1 = 0.75d0
  real*8 :: omega2 = 10.0d0

contains

  subroutine deallocate_ppm
    implicit none
    
    deallocate(delta_r)
    deallocate(c_1_r,c_2_r,c_3_r,c_4_r,c_5_r)
    
  end subroutine deallocate_ppm

  subroutine ppm_coefficients
    ! using some global module vars
    use GR1D_module, only: x1i,x1,ghosts1,n1
    implicit none
    
    integer i
    real*8, allocatable ::  tmp_1(:)
    real*8, allocatable ::  tmp_2(:)
    real*8, allocatable ::  tmp_3(:)
    real*8 tmp

    allocate(c_1_r(n1))
    allocate(c_2_r(n1))
    allocate(c_3_r(n1))
    allocate(c_4_r(n1))
    allocate(c_5_r(n1))

    c_1_r(:) = 0.0d0
    c_2_r(:) = 0.0d0
    c_3_r(:) = 0.0d0
    c_4_r(:) = 0.0d0
    c_5_r(:) = 0.0d0

    allocate(tmp_1(n1))
    allocate(tmp_2(n1))
    allocate(tmp_3(n1))
    allocate(delta_r(n1))

    tmp_1(:) = 0.0d0
    tmp_2(:) = 0.0d0
    tmp_3(:) = 0.0d0
    delta_r(:) = 0.0d0

    do i=ghosts1+1,n1-1
       delta_r(i) = x1i(i+1) - x1i(i)
    enddo
    delta_r(n1) = x1i(n1) + 2.0d0*(x1(n1) - x1i(n1)) !never is used, wrong I think
    delta_r(ghosts1) = delta_r(ghosts1+1)
    delta_r(ghosts1-1) = delta_r(ghosts1+2)
    delta_r(ghosts1-2) = delta_r(ghosts1+3)
    delta_r(ghosts1-3) = delta_r(ghosts1+4)

    do i=ghosts1-2,n1-1 
       tmp_1(i) = delta_r(i-1) + delta_r(i)
       tmp_2(i) = tmp_1(i) + delta_r(i)
       tmp_3(i) = tmp_1(i) + delta_r(i-1)
    enddo

    ! coefficients for delta a (1.7)
    do i=ghosts1-1,n1-2
       tmp = delta_r(i) / &
            (delta_r(i - 1) + delta_r(i) + delta_r(i + 1))
       c_1_r(i) = tmp * tmp_3(i) / tmp_1(i + 1)
       c_2_r(i) = tmp * tmp_2(i + 1) / tmp_1(i)
    enddo


    ! coefficients for delta a_i+1, delta a_i, and
    ! a_i+1 - a_i (1.6)
    do i=ghosts1-1,n1-3
       tmp = 1.0d0 / (tmp_1(i) + tmp_1(i+2))
       c_3_r(i) = -tmp * delta_r(i) * tmp_1(i) / tmp_3(i+1)
       c_4_r(i) = tmp * delta_r(i+1) * tmp_1(i+2) / tmp_2(i+1)
       c_5_r(i) = (delta_r(i) - 2.0d0 * (delta_r(i+1) * c_3_r(i) &
            + delta_r(i) * c_4_r(i))) / tmp_1(i+1)
    enddo
    
    deallocate(tmp_1)
    deallocate(tmp_2)
    deallocate(tmp_3)

  end subroutine ppm_coefficients

  subroutine ppm_interpolate(q,qp,qm)
    
    use GR1D_module, only: n1,ghosts1,x1,x1i,ppm_origin_TVD
    implicit none
    integer i
    real*8 q(n1),qp(n1),qm(n1)
    real*8 delta_q(n1) ! this is put on the stack
    real*8 ppm_tmp(n1) ! this is put on the stack
    real*8 delta_q_sign
    ! for TVD near the origin:
    real*8 dupw, dloc, delta
    real*8 dupwdx, dlocdx, dxl, dxr, dx    

    ppm_tmp(:) = 0.0d0
    qm(:) = 0.0d0
    qp(:) = 0.0d0
    delta_q(:) = 0.0d0


    ! a_i - a_i-1 (for delta a_i) (1.7)
    ! ppm_tmp(i-1) = a_i - a_i-1
    ! ppm_tmp(i) = a_i+1 - a_i
    do i = ghosts1-2,n1-1
       ppm_tmp(i) = q(i+1) - q(i)
    enddo

    do i = ghosts1-1,n1-2
    ! delta a_i (1.7)
       delta_q(i) = c_1_r(i) * ppm_tmp(i) &
            + c_2_r(i) * ppm_tmp(i-1)
    ! convert delta a_i to delta_m a_i (1.8)
       if( ppm_tmp(i)*ppm_tmp(i-1) .gt. 0.0d0) then
          delta_q_sign = sign(1.0d0,delta_q(i))
          delta_q(i) = &
               min( abs(delta_q(i)), &
                    2.0d0*min( abs(ppm_tmp(i)), &
                               abs(ppm_tmp(i-1)))) &
              * delta_q_sign
       else 
          delta_q(i) = 0.0d0
       endif
    enddo

    ! a_i+1/2 (1.6)
    do i=ghosts1-1,n1-3
       qp(i) = q(i) + &
            c_5_r(i) * ppm_tmp(i) + &
            c_3_r(i) * delta_q(i+1) + &
            c_4_r(i) * delta_q(i)
       qm(i+1) = qp(i)
    enddo

    if(ppm_origin_TVD.gt.0) then
       do i=ghosts1-1,ghosts1+(ppm_origin_TVD+1)
          
          dupwdx = x1(i) - x1(i-1)
          dlocdx = x1(i+1) - x1(i)
          dxr = x1i(i+1) - x1(i)
          dxl = x1(i) - x1i(i)
          
          dupw = (q(i) - q(i-1))/dupwdx
          dloc = (q(i+1) - q(i))/dlocdx
          
          ! MC Lim
          if (dupw*dloc < 0.d0) then
             delta=0.d0
          else 
             delta=sign(min(2.0d0*abs(dupw),2.0d0*abs(dloc),&
                  0.5d0*(abs(dupw)+abs(dloc))),(dupw+dloc))
          end if
          
          qm(i) = q(i) - dxl*delta
          qp(i) = q(i) + dxr*delta
       enddo
    endif

  end subroutine ppm_interpolate

  subroutine ppm_steepen
    !Steepening of the density profile
    !near contact discontiuities.
    use GR1D_module,only: rho,rhop,rhom, &
         n1,ghosts1,x1,x1i,press,cs2,eps, GR

    integer i
    real*8 delta_2_rho_minus,delta_2_rho_plus
    real*8 eta, eta_tilde
    real*8 criterion1,criterion2,criterion3,criterion4
    real*8 tmp, gamma_eos_tmp
    real*8 :: tiny = 1.0d-10
    real*8 delta_q(n1) ! put on the stack
    real*8 delta_a(n1), alphaleftD, alpharightD
    
    ! first set up delta_rho:
    do i=ghosts1-2,n1-2
       delta_q(i) = c_1_r(i)*(rho(i+1)-rho(i)) &
            + c_2_r(i) * (rho(i) - rho(i-1))
    enddo


    do i=ghosts1-1,n1-3
       ! C&W (1.17)
       delta_2_rho_minus = & 
            ( (rho(i) - rho(i-1)) / &
              (delta_r(i-1) + delta_r(i)) &
             - (rho(i-1) - rho(i-2)) / &
              (delta_r(i-2) + delta_r(i-1)) ) / &
              (delta_r(i-2)+delta_r(i-1)+delta_r(i)) 

       delta_2_rho_plus = &
            ( (rho(i+2) - rho(i+1)) / &
              (delta_r(i+1) + delta_r(i+2)) &
             - (rho(i+1) - rho(i)) / &
              (delta_r(i) + delta_r(i+1)) ) / &
              (delta_r(i) + delta_r(i+1)+delta_r(i+2))

       tmp = rho(i+1) - rho(i-1)
       if (abs(tmp).lt.tiny) tmp = tiny

       criterion1 = delta_2_rho_plus * delta_2_rho_minus
       criterion2 = abs(tmp) - epsilon1 * &
            min(abs(rho(i+1)),abs(rho(i-1)))
       if(criterion1.gt.0.0d0.or.criterion2.lt.0.0d0) then
          eta_tilde = 0.0d0
       else
          eta_tilde = &
               (delta_2_rho_minus - delta_2_rho_plus) * &
               ( (x1(i)-x1(i-1))**3 + (x1(i+1)-x1(i))**3 ) / &
               ( (x1(i+1) - x1(i-1)) * tmp )
       endif
!       (1.16)
       eta = max(0.0d0,min(eta1*(eta_tilde-eta2),1.0d0))
       ! Now compute the gamma needed for
       ! for shock detection.
       ! The code goes along the lines of (3.2) in C&W
       if (GR) then
          gamma_eos_tmp = cs2(i)*rho(i)*(1.0d0+press(i)/rho(i)+eps(i))/press(i)
       else
          gamma_eos_tmp = cs2(i)*rho(i)/press(i)
       endif
       criterion3 = gamma_eos_tmp * K0 * &
            abs(rho(i+1)-rho(i-1))/min(rho(i+1),rho(i-1))
       criterion4 = abs(press(i+1)-press(i-1)) / &
            min(press(i+1),press(i-1))
       if(criterion3.lt.criterion4) eta = 0.0d0
       ! Now compute the steepened density profile.
       ! Note that this also includes the modification of
       ! (1.14) in C&W.
       rhom(i) = rhom(i) * (1.0d0 - eta) + &
            (rho(i-1) + 0.5d0*delta_q(i-1))*eta

       rhop(i) = rhop(i) * (1.0d0 - eta) + &
            (rho(i+1) - 0.5d0 * delta_q(i+1))*eta

    enddo
    
  end subroutine ppm_steepen

  subroutine ppm_flatten
    ! Flatten profile near strong shocks.
    ! This is somewhat similar to adding
    ! artificial dissipation. See Appendix and
    ! section 4 of C&W.

    ! using a lot of module vars:
    use GR1D_module, only: press,v1,v,      &
         v1p,v1m,vp,vm,eps,epsp,epsm,ye,   &
         yem,yep,rho,rhop,rhom,n1,ghosts1, &
         x1, geometry

    implicit none
    integer i
    real*8 criterion,w,tmp,f
    real*8 delta_v,delta_p_2
    real*8 delta_p(n1)  ! on the stack...
    real*8 f_tilde(n1)  ! on the stack...
    real*8 :: tiny = 1.0d-40

    delta_p(:) = 0.0d0
    f_tilde(:) = 0.0d0

    do i = ghosts1-1,n1-2
       delta_p(i) = press(i+1) - press(i-1)
       if (geometry.eq.2) then
          delta_v = v1(i+1)*x1(i+1)**2 - v1(i-1)*x1(i-1)**2
       else if (geometry.eq.1) then
          delta_v = v1(i+1) - v1(i-1)          
       else
          stop "Add in velocity condition in ppm"
       endif
       ! slightly rewritten A.1
       criterion = epsilon2 * min(press(i+1),press(i-1)) - &
            abs(delta_p(i))
       if(criterion.lt.0.0d0.and.delta_v.lt.0.0d0) then
       ! set w_i = 1
          w = 1.0d0
       else
          w = 0.0d0
       endif

       ! now calculate A.2
       delta_p_2 = press(i+2) - press(i-2)
       ! exception handling (not in paper)
       if(abs(delta_p_2).lt.tiny) then
          if(abs(delta_p(i)).lt.tiny) then
             tmp = - omega1
          else
             tmp = 1.0d0 - omega1
          endif
       else
          tmp = delta_p(i) / delta_p_2 - omega1
       endif
       ! A.2 after some algebra
       f_tilde(i) = min(1.0d0,w*max(0.0d0,tmp*omega2))
    enddo

    ! now we apply the flattening, but first modify
    ! f_tilde according to the shock travel direction
    do i=ghosts1,n1-3
       if(delta_p(i).lt.0.0d0) then
          f = max(f_tilde(i),f_tilde(i+1))
       else
          f = max(f_tilde(i),f_tilde(i-1))
       endif

    ! The flattening is done according to
    ! (4.1) in the C&W paper. NO additional dissipative
    ! flux is added to the velocities.
       
       rhom(i) = rho(i)* f + rhom(i)*(1.0d0 - f)
       rhop(i) = rho(i)* f + rhop(i)*(1.0d0 - f)

       v1m(i)  = v1(i)*f + v1m(i)*(1.0d0 - f)
       v1p(i)  = v1(i)*f + v1p(i)*(1.0d0 - f)

       vm(i)  = v(i)*f + vm(i)*(1.0d0 - f)
       vp(i)  = v(i)*f + vp(i)*(1.0d0 - f)

       epsm(i) = eps(i)*f + epsm(i)*(1.0d0 - f)
       epsp(i) = eps(i)*f + epsp(i)*(1.0d0 - f)

!       yem(i) = ye(i)* f + yem(i)*(1.0d0 - f)
!       yep(i) = ye(i)* f + yep(i)*(1.0d0 - f)

    enddo

  end subroutine ppm_flatten

  subroutine ppm_monotonize(q,qp,qm)
    
    use GR1D_module, only: n1, ghosts1
    implicit none
    integer i
    real*8 q(n1),qp(n1),qm(n1)

! this implements (1.10)
    do i=ghosts1-1,n1-2
       !case 1
       if( (qp(i) - q(i))*(qm(i)-q(i)) .ge. 0.0d0) then
          qp(i) = q(i)
          qm(i) = q(i)
       else 
       !case 2
          if ( (qp(i)-qm(i))*(q(i)-0.5d0*(qm(i)+qp(i))).gt. &
               (qp(i)-qm(i))**2/6.0d0) then
             qm(i) = 3.0d0*q(i) - 2.0d0*qp(i)
          endif
       !case 3
          if ( (qp(i)-qm(i))*(q(i)-0.5d0*(qm(i)+qp(i))).lt. &
               -(qp(i)-qm(i))**2/6.0d0) then
             qp(i) = 3.0d0*q(i) - 2.0d0*qm(i)
          endif
       endif
    enddo
    
  end subroutine ppm_monotonize

end module ppm

