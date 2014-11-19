!-*-f90-*-
subroutine mass_interior
        
! calculate baryonic mass interior to cell centers
! for innermost zone: use that zone inner radius and assume that its
! density is constant interior to the inner cell radius
  
  use GR1D_module
  implicit none
  
  integer i
  real*8 dmass,dx1,dx2,discrim

  real*8 h,rhot,massint

  mass(:) = 0.0d0
  if (GR) then	
     mass1 = X * rho * W * volume
  else
     mass1 = rho * volume
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
     mass(ghosts1+1) = mass1(ghosts1+1) - 4.0d0/3.0d0*pi*rho(ghosts1+1)* &
          ( x1i(ghosts1+1+1)**3 - x1(ghosts1+1)**3 )
     
     do i=ghosts1+2,n1-1
        mass(i) = mass(i-1) + mass1(i) &
             - 4.0d0/3.0d0*pi*rho(i) * ( x1i(i+1)**3 - x1(i)**3 )  &
             + 4.0d0/3.0d0*pi*rho(i-1) * (x1i(i)**3 - x1(i-1)**3)
     enddo
     
     mass(n1) = mass(n1-1) &
          + 4.0d0/3.0d0*pi*rho(n1-1)*(x1i(n1)**3 - x1(n1-1)**3) &
          + 4.0d0/3.0d0*pi*rho(n1)*( x1(n1)**3 - x1i(n1)**3 ) 
  endif

end subroutine mass_interior

