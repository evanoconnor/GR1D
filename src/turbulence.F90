!-*-f90-*-
subroutine turb_diff_terms

   use GR1D_module, only: x1, eps, ye, v_turb, alpha_turb, &
     alpha_turb_K, alpha_turb_e, alpha_turb_ye, dphidr, rhop, &
     pressp, v_turbp, qp, diff_term_eps, diff_term_ye, &
     diff_term_K, n1, ghosts1
   use Grad_module
   implicit none

   integer:: i
   real*8 :: Lambda_mixp
   real*8 :: D_turb_eps, D_turb_ye, D_turb_K
   real*8 :: eps_grad(n1), ye_grad(n1), v2_turb_grad(n1)

   eps_grad = Gradient_int(eps,x1)
   ye_grad = Gradient_int(ye,x1)
   v2_turb_grad = Gradient_int(v_turb**2,x1)
   
   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      Lambda_mixp = alpha_turb * pressp(i) / (rhop(i) * dphidr(i))
      Lambda_mixp = min(Lambda_mixp,x1(i))
      
      ! Calculate diffusion coefficients as given in Couch et al. 2019 eqs. 29-31
      D_turb_eps = alpha_turb_e * v_turbp(i) * Lambda_mixp
      D_turb_ye = alpha_turb_ye * v_turbp(i) * Lambda_mixp
      D_turb_K = alpha_turb_K * v_turbp(i) * Lambda_mixp
   
      diff_term_eps(i) = qp(i,1) * D_turb_eps * eps_grad(i)
      diff_term_ye(i) = qp(i,1) * D_turb_ye * ye_grad(i)
      diff_term_K(i) = qp(i,1) * D_turb_K * v2_turb_grad(i)
   enddo 

end subroutine turb_diff_terms

subroutine turbulence_sources

   use GR1D_module, only: x1, GR, v_turb, omega2_BV, v1, &
     alpha_turb, dphidr, press, rho, n1, ghosts1, turb_source, &
     alp, X, shear, diss, buoy, sqrt_gamma, ishock
   use Grad_module
   implicit none

   integer i
   real*8 Lambda_mix
   real*8 v_grad(n1)
   
   v_grad = Gradient_5pts(v1,x1)
   
   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      Lambda_mix = alpha_turb * press(i)/(rho(i)*dphidr(i))
      Lambda_mix = min(Lambda_mix, x1(i))
 
      shear(i) = - v_turb(i)**2 * v_grad(i)
      diss(i) = v_turb(i)**3 / Lambda_mix
      buoy(i) = v_turb(i) * omega2_BV(i) * Lambda_mix

      turb_source(i,3) = rho(i) * diss(i)
      turb_source(i,6) = rho(i)* (shear(i) + buoy(i) - diss(i))

      if (GR) then
          turb_source(i,3) = alp(i)*X(i)*turb_source(i,3)
          turb_source(i,6) = alp(i)*X(i)*turb_source(i,6)
      else
          turb_source(i,3) = sqrt_gamma(i)*turb_source(i,3)
          turb_source(i,6) = sqrt_gamma(i)*turb_source(i,6)
      endif
   enddo
 
end subroutine turbulence_sources

subroutine Brunt_Vaisala(dts)

   use GR1D_module, only: x1, v_turb, omega2_BV, GR, &
     alpha_turb, Lambda_MLT, v1, v, dphidr, press, rho, eps, cs2, &
     n1, ghosts1, alp, X, ishock, length_gf
   use Grad_module
   implicit none
   
   !local
   integer i 
   real*8 dts
   real*8 Lambda_mix, H_P, v_turb_seed, h
   real*8 rho_grad(n1), press_grad(n1), v_grad(n1)
   real*8 lnrho_grad(n1), lnP_grad(n1)
   
   if (GR) then
       rho_grad = Gradient_3pts(rho(:)*(1.0d0 + eps(:)),x1)
       press_grad = Gradient_3pts(press,x1)
       v_grad = Gradient_5pts(v,x1)
   else
       rho_grad = Gradient_3pts(rho,x1)
       press_grad = Gradient_3pts(press,x1)   
       v_grad = Gradient_5pts(v1,x1)
       
       !lnrho_grad = Gradient_3pts(log(rho),x1)
       !lnP_grad = Gradient_3pts(log(press),x1)
   endif
   
   do i=ghosts1+1,n1-ghosts1
       
      if (GR) then
          h = 1.0d0 + eps(i) + press(i)/rho(i)
          
          omega2_BV(i) = alp(i)**2/(rho(i)*h*X(i)**2)* &
                (dphidr(i) - v1(i)* v_grad(i)) * &
                (rho_grad(i) - press_grad(i)/cs2(i))
      else
          omega2_BV(i) = (dphidr(i) - v(i)* v_grad(i))/rho(i)* &
                (rho_grad(i) - press_grad(i)/cs2(i))
      endif
      
      H_P = press(i)/(rho(i)*dphidr(i))
      Lambda_mix = alpha_turb * H_P
      Lambda_mix = min(Lambda_mix,x1(i))

      Lambda_MLT(i) = Lambda_mix !store for ouptut

      if (omega2_BV(i) .gt. 0.0d0) then
          v_turb_seed = dts*omega2_BV(i)*Lambda_mix
          v_turb(i) = max(v_turb(i),v_turb_seed)
      endif

      if ( x1(i) / length_gf .lt. 5.0d5 ) v_turb(i) = 0.0d0

      ! force covective velocity to be zero in the 5 points stencil around the shock,
      ! i.e. do not convect through the shock
      if ( ishock(1)-2 <= i-ghosts1 .and. i-ghosts1 <= ishock(1)+2 ) v_turb(i) = 0.0d0

   enddo

end subroutine Brunt_Vaisala
