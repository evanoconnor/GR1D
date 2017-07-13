!-*-f90-*-
subroutine mass_interior
        
! calculate baryonic mass interior to cell centers
! for innermost zone: use that zone inner radius and assume that its
! density is constant interior to the inner cell radius
  
  use GR1D_module
  implicit none
  
  integer i
  real*8 dmass,dx1,dx2,discrim

  real*8 h,rhot,massint,phi_bound
  real*8 include_nus
  
  mass(:) = 0.0d0
  if (GR) then	
     mass1 = X * rho * W * volume
  else
     mass1 = rho * volume
  endif

  if (do_M1) then
     include_nus = 1.0d0
  else
     include_nus = 0.0d0
  endif
  
  if (GR) then
     mass(ghosts1+1) = mass1(ghosts1+1) - 4.0d0/3.0d0*pi*rho(ghosts1+1) * X(ghosts1+1) * &
          W(ghosts1+1)*( x1i(ghosts1+1+1)**3 - x1(ghosts1+1)**3 )
     
     do i=ghosts1+2,n1-1
        mass(i) = mass(i-1) + mass1(i) &
             - 4.0d0/3.0d0*pi*rho(i)*X(i)*W(i) * ( x1i(i+1)**3 - x1(i)**3 )  &
             + 4.0d0/3.0d0*pi*rho(i-1)*X(i-1)*W(i-1) * (x1i(i)**3 - x1(i-1)**3)
     enddo
     
     mass(n1) = mass(n1-1) &
          + 4.0d0/3.0d0*pi*rho(n1-1)*X(n1-1)*W(n1-1)*(x1i(n1)**3 - x1(n1-1)**3) &
          + 4.0d0/3.0d0*pi*rho(n1)*X(n1)*W(n1)*( x1(n1)**3 - x1i(n1)**3 ) 
  else
     if (do_effectivepotential) then
        mass(ghosts1+1) = 4.0d0*pi/3.0d0*x1(ghosts1+1)**3*(rho(ghosts1+1)* &
             (1.0d0 + eps(ghosts1+1)) + (energy_nu(ghosts1+1)+mom_nu(ghosts1+1)*v1(ghosts1+1))*include_nus)
        mass(ghosts1+1) = mass(ghosts1+1)*sqrt(1.0d0-2.0d0*mass(ghosts1+1)/x1(ghosts1+1))

        dphidr(ghosts1+1) = (mass(ghosts1+1) + 4.0d0*pi*x1(ghosts1+1)**3* &
             (press(ghosts1+1)+press_nu(ghosts1+1)))/ &
             (x1(ghosts1+1)**2*(1.0d0+v1(ghosts1+1)**2-2.0d0*mass(ghosts1+1)/x1(ghosts1+1)))* &
             (rho(ghosts1+1)+eps(ghosts1+1)*rho(ghosts1+1)+press(ghosts1+1))/rho(ghosts1+1)

        do i=ghosts1+2,n1-1
           mass(i) = mass(i-1) + &
                4.0d0/3.0d0*pi*(rho(i-1)*(1.0d0+eps(i-1))+(energy_nu(i-1)+mom_nu(i-1)*v1(i-1))*include_nus) &
                * ( x1i(i)**3 - x1(i-1)**3 )*sqrt(1.0d0-2.0d0*mass(i-1)/x1(i-1))
           mass(i) = mass(i) + &
                4.0d0/3.0d0*pi*(rho(i)*(1.0d0+eps(i))+(energy_nu(i)+mom_nu(i)*v1(i))*include_nus) * &
                (x1(i)**3 - x1i(i)**3)*sqrt(1.0d0-2.0d0*mass(i)/x1i(i))
           dphidr(i) = (mass(i) + 4.0d0*pi*x1(i)**3*(press(i)+press_nu(i)))/ &
                (x1(i)**2*(1.0d0+v1(i)**2-2.0d0*mass(i)/x1(i)))* &
                (rho(i)+eps(i)*rho(i)+press(i))/rho(i)
        enddo
     
        mass(n1) = mass(n1-1) + &
             4.0d0/3.0d0*pi*(rho(n1-1)*(1.0d0+eps(n1-1))+(energy_nu(n1-1)+mom_nu(n1-1)*v1(n1-1))*include_nus)* &
             (x1i(n1)**3 - x1(n1-1)**3)*sqrt(1.0d0-2.0d0*mass(n1-1)/x1(n1-1))
        mass(n1) = mass(n1) + &
             4.0d0/3.0d0*pi*(rho(n1)*(1.0d0+eps(n1))+(energy_nu(n1)+mom_nu(n1)*v1(n1))*include_nus)* &
             (x1(n1)**3 - x1i(n1)**3)*sqrt(1.0d0-2.0d0*mass(n1)/x1i(n1))         
        dphidr(n1) = (mass(n1) + 4.0d0*pi*x1(n1)**3*(press(n1)+press_nu(n1)))/ &
             (x1(n1)**2*(1.0d0+v1(n1)**2-2.0d0*mass(n1)/x1(n1)))* &
             (rho(n1)+eps(n1)*rho(n1)+press(n1))/rho(n1)    

        phi(ghosts1+1) = 0.0d0 + x1(ghosts1+1)*(dphidr(ghosts1+1))
        phii(ghosts1+1) = 0.0d0
        do i=ghosts1+2,n1-ghosts1+1
           phi(i) = phi(i-1) + (x1i(i)-x1(i-1))*dphidr(i-1) &
                + (x1(i)-x1i(i))*dphidr(i)
           phii(i) = phi(i) - (x1(i)-x1i(i))*dphidr(i)
        enddo
        phii(n1-ghosts1+2) = phi(n1-ghosts1+1) + &
             (x1i(n1-ghosts1+2)-x1(n1-ghosts1+1))*dphidr(n1-ghosts1+1)
        
        phi_bound = 0.5d0 * &
             log(1.0d0 - 2.0d0 * mass(n1-ghosts1+1)/x1(n1-ghosts1+1))
  
        do i=ghosts1+1,n1-ghosts1+1
           phi(i) = phi(i) + (phi_bound-phi(n1-ghosts1+1))
           phii(i) = phii(i) + (phi_bound-phi(n1-ghosts1+1))
        enddo

        do i=ghosts1+1,n1-ghosts1+1
           alp(i) = exp(phi(i))
           alpm(i) = exp(phii(i))
           alpp(i) = exp(phii(i+1))
        enddo

        dphidr(ghosts1+1) = (phi(ghosts1+2)-phi(ghosts1+1))/(x1(ghosts1+2)-x1(ghosts1+1))
        do i=ghosts1+2,n1-ghosts1
           dphidr(i) = (phi(i+1)-phi(i-1))/(x1(i+1)-x1(i-1))
        enddo

     else
        mass(ghosts1+1) = mass1(ghosts1+1) - 4.0d0/3.0d0*pi*rho(ghosts1+1)* &
             ( x1i(ghosts1+1+1)**3 - x1(ghosts1+1)**3 )
        dphidr(ghosts1+1) =  mass(ghosts1+1)/x1(ghosts1+1)**2

        do i=ghosts1+2,n1-1
           mass(i) = mass(i-1) + mass1(i) &
                - 4.0d0/3.0d0*pi*rho(i) * ( x1i(i+1)**3 - x1(i)**3 )  &
                + 4.0d0/3.0d0*pi*rho(i-1) * (x1i(i)**3 - x1(i-1)**3)
           dphidr(i) =  mass(i)/x1(i)**2
        enddo
     
        mass(n1) = mass(n1-1) &
             + 4.0d0/3.0d0*pi*rho(n1-1)*(x1i(n1)**3 - x1(n1-1)**3) &
             + 4.0d0/3.0d0*pi*rho(n1)*( x1(n1)**3 - x1i(n1)**3 ) 
        dphidr(n1) =  mass(n1)/x1(n1)**2
     endif
  endif

end subroutine mass_interior

