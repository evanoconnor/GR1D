module Grad_module

  use GR1D_module, only: ghosts1, n1
  implicit none
  public

contains

function Gradient(var,r) result(grad)
 
  real*8 var(n1),r(n1)
  integer i
  real*8 grad(n1), dumgrad

  grad(:) = 0.0d0
  ! calculate derivativees with standard gradient formula
  !do i=ghosts1+2,n1-ghosts1-1
  !   grad(i) = (var(i+1) - var(i-1))/(r(i+1) - r(i-1))
  !enddo
  
  ! calculate derivatives with the 3 point formula
  do i=ghosts1+2,n1-ghosts1-1
     grad(i) = var(i+1)*(r(i)-r(i-1))/((r(i+1)-r(i-1))*(r(i+1)-r(i))) &
             + var(i)*(r(i+1)-2.0d0*r(i)+r(i-1))/((r(i+1)-r(i))*(r(i)-r(i-1))) &
             - var(i-1) * (r(i+1)-r(i))/((r(i)-r(i-1))*(r(i+1)-r(i-1))) 
  enddo
  
  ! calculate one-sided derivatives
  grad(ghosts1+1) = (var(ghosts1+2) - var(ghosts1+1))/ &
                  (r(ghosts1+2) - r(ghosts1+1))

  grad(n1-ghosts1) = (var(n1-ghosts1) - var(n1-ghosts1-1))/ &
                     (r(n1-ghosts1) - r(n1-ghosts1-1))

end function Gradient 

function Gradient_3pts(var,r) result(grad)

  real*8 var(n1),r(n1)
  integer i
  real*8 grad(n1)
  real*8 h1,h2

  grad(:) = 0.0d0
  
  ! calculate derivatives with the 3 point formula
  do i=ghosts1+2,n1-ghosts1-1
     
     h1 = r(i)   - r(i-1)
     h2 = r(i+1) - r(i)
    
     grad(i) = - var(i-1)*(h2/(h1*(h1+h2))) &
               + var(i  )*((h2-h1)/(h1*h2)) &
               + var(i+1)*((h1/(h2*(h1+h2))))
  
  enddo
                    
  ! calculate one-sided derivatives
  grad(ghosts1+1) = (var(ghosts1+2) - var(ghosts1+1))/ &
                  (r(ghosts1+2) - r(ghosts1+1))

  grad(n1-ghosts1) = (var(n1-ghosts1) - var(n1-ghosts1-1))/ &
                     (r(n1-ghosts1) - r(n1-ghosts1-1))

end function Gradient_3pts

function Gradient_5pts(var,r) result(grad)

  real*8 var(n1),r(n1)
  integer i
  real*8 grad(n1)
  real*8 h1,h2,h3,h4,Ht1,Ht2

  grad(:) = 0.0d0
  
  ! calculate derivatives with the 5 point formula
  do i=ghosts1+3,n1-ghosts1-2
     
     h1 = r(i-1) - r(i-2)
     h2 = r(i)   - r(i-1)
     h3 = r(i+1) - r(i)
     h4 = r(i+2) - r(i+1)
     Ht1 = h1+h2+h3
     Ht2 = Ht1+h4
    
     grad(i) = var(i-2)*(h2*h3*(h3+h4))/(h1*(h1+h2)*Ht1*Ht2) &
             - var(i-1)*((h1+h2)*h3*(h3+h4))/(h1*h2*(h2+h3)*(Ht2-h1)) &
             + var(i)*((h1+2.0d0*h2)*h3*(h3+h4)-(h1+h2)*h2*(2.0d0*h3+h4)) &
                     /((h1+h2)*h2*h3*(h3+h4)) &
             + var(i+1)*((h2+h1)*h2*(h3+h4))/(Ht1*(h2+h3)*h3*h4) &
             - var(i+2)*(h2*(h1+h2)*h3)/(Ht2*(Ht2-h1)*(h3+h4)*h4)

  enddo
  
  ! simple 2 points far away from the shock should suffice
  grad(ghosts1+2) = (var(ghosts1+3) - var(ghosts1+1))/ &
                  (r(ghosts1+3) - r(ghosts1+1))
  grad(n1-ghosts1-1) = (var(n1-ghosts1) - var(n1-ghosts1-2))/ &
                     (r(n1-ghosts1) - r(n1-ghosts1-2))
                     
  ! calculate one-sided derivatives
  grad(ghosts1+1) = (var(ghosts1+2) - var(ghosts1+1))/ &
                  (r(ghosts1+2) - r(ghosts1+1))

  grad(n1-ghosts1) = (var(n1-ghosts1) - var(n1-ghosts1-1))/ &
                     (r(n1-ghosts1) - r(n1-ghosts1-1))

end function Gradient_5pts
 
function Gradient_int(var,radius) result(grad)

  real*8, intent(in) :: var(n1),radius(n1)
  integer i
  real*8 grad(n1)
  
  grad(:) = 0.0d0
  ! calculate derivatives
  do i=ghosts1+1,n1-ghosts1-1
     grad(i) = (var(i+1) - var(i))/(radius(i+1) - radius(i))
  enddo

  grad(n1-ghosts1) = grad(n1-ghosts1-1)

end function Gradient_int

subroutine QuadraticInterpolation1D(xOld, yOld, nX, &
                                 X_New, Y_new)

  real*8, intent(out)   :: Y_New
  real*8, intent(in)    :: X_New
  real*8, intent(in)    :: xOld(nX)
  real*8, intent(in)    :: yOld(nX)
  integer, intent(in)     :: nX

  integer  :: i
  real*8 :: L2_0, L2_1, L2_2

  if (X_New <= xOld(1)) then

    !Extrapolate DOwn
    Y_New = yOld(1) - (xOld(1) - X_New) * (yOld(2) - yOld(1)) &
                                          / (xOld(2) - xOld(1))

    return

  else if ( X_New <= xOld(2) ) then

    Y_New = yOld(1) + (X_New    - xOld(1)) &
                     * (yOld(2) - yOld(1)) &
                     / (xOld(2) - xOld(1))

  else if (X_New >= xOld(nX)) then

    !Exptrapolate up
    Y_New = yOld(nX) + (X_New - xOld(nX)) * (yOld(nX) - yOld(nX-1)) &
                                            / (xOld(nX) - xOld(nX-1))
    return

  else if ( X_New >= xOld(nX-1) ) then

    Y_New = yOld(nX-1) + (X_New     - xOld(nX-1)) &
                        * (yOld(nX) - yOld(nX-1)) &
                        / (xOld(nX) - xOld(nX-1))

  else

    do i = 1,nX

      if( xOld(i) >= X_New ) then

        L2_0 = ( (X_New      - xOld(i)) * (X_New      - xOld(i+1)) ) / &
               ( (xOld(i-1) - xOld(i)) * (xOld(i-1) - xOld(i+1)) )
        L2_1 = ( (X_New    - xOld(i-1)) * (X_New    - xOld(i+1)) ) / &
               ( (xOld(i) - xOld(i-1)) * (xOld(i) - xOld(i+1)) )

        L2_2 = ( (X_New      - xOld(i-1)) * (X_New      - xOld(i)) ) / &
               ( (xOld(i+1) - xOld(i-1)) * (xOld(i+1) - xOld(i)) )

        Y_New = yOld(i-1)*L2_0 + yOld(i)*L2_1 + yOld(i+1)*L2_2

        return

      end if

    end do

  end if

end subroutine QuadraticInterpolation1D

end module Grad_module
