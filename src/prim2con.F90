!-*-f90-*-
subroutine prim2con

  use GR1D_module
  implicit none
  
  call prim2con_1
  
end subroutine prim2con


subroutine prim2con_1
  
  use GR1D_module
  implicit none
  
  integer i
  real*8 h
  
  if (GR) then 
     !Romero's formulation
     do i=1,n1-1
        h = (1.0d0 + press(i)/rho(i) + eps(i) )
        q(i,1) = X(i)*W(i)*rho(i)
        q(i,2) = rho(i)*h*W(i)**2*v(i)
        q(i,3) = rho(i)*h*W(i)**2 - press(i) - q(i,1) 
        q(i,4) = X(i)*W(i)*rho(i)*ye(i)

        if(do_rotation) then
           q(i,5) = rho(i)*h*W(i)**2*vphi(i)*x1(i)
        endif
     enddo
  else 
     !Newtonian
     do i=1,n1
        q(i,1) = rho(i)
        q(i,2) = rho(i)*v1(i)
        q(i,3) = rho(i)*eps(i) + 0.5d0*rho(i)*v1(i)**2 
        if(do_rotation) then
           q(i,3) = q(i,3) + 0.5d0*rho(i)*twothirds*vphi1(i)**2
        endif
        q(i,4) = rho(i)*ye(i)
     enddo
     if(do_rotation) then
        do i=1,n1
           q(i,5) = rho(i)*vphi1(i)*x1(i)
        enddo
     endif

  endif
  
end subroutine prim2con_1

subroutine prim2con_if
  
  use GR1D_module
  
  implicit none
  integer i
  real*8 hp,hm
  
  if (GR) then 
     do i=ghosts1-1,n1-ghosts1+1

        hp = (1.0d0 + pressp(i)/rhop(i)+ epsp(i))
        hm = (1.0d0 + pressm(i)/rhom(i)+ epsm(i))

        qp(i,1) = Xp(i) * Wp(i) * rhop(i)
        qp(i,2) = rhop(i) * hp * Wp(i)**2 * vp(i)
        qp(i,3) = rhop(i) * hp * Wp(i)**2 - pressp(i) - qp(i,1) 
        qp(i,4) = Xp(i) * Wp(i) * rhop(i) * yep(i)
        
        qm(i,1) = Xm(i) * Wm(i) * rhom(i)
        qm(i,2) = rhom(i) * hm * Wm(i)**2 * vm(i)
        qm(i,3) = rhom(i) * hm * Wm(i)**2 - pressm(i) - qm(i,1) 
        qm(i,4) = Xm(i) * Wm(i) * rhom(i) * yem(i)

        if(do_rotation) then
           qp(i,5) = rhop(i) * hp * Wp(i)**2 * vphip(i) * x1i(i+1)
           qm(i,5) = rhom(i) * hm * Wm(i)**2 * vphim(i) * x1i(i)
        endif

     enddo
  else 
     !Newtonian
     do i=ghosts1-1,n1-ghosts1+1
        qp(i,1) = rhop(i)
        qp(i,2) = rhop(i)*v1p(i)
        qp(i,3) = rhop(i)*epsp(i) + 0.5d0*rhop(i)*v1p(i)**2  
        qp(i,4) = rhop(i)*yep(i)
       
        qm(i,1) = rhom(i)
        qm(i,2) = rhom(i)*v1m(i)
        qm(i,3) = rhom(i)*epsm(i) + 0.5d0*rhom(i)*v1m(i)**2 
        qm(i,4) = rhom(i)*yem(i)

        if(do_rotation) then
           qp(i,5) = rhop(i)*vphi1p(i)*x1i(i+1)
           qm(i,5) = rhom(i)*vphi1m(i)*x1i(i)
           qp(i,3) = qp(i,3) + 0.5d0*rhop(i)*twothirds*vphi1p(i)**2  
           qm(i,3) = qm(i,3) + 0.5d0*rhom(i)*twothirds*vphi1m(i)**2 
        endif

     enddo
  endif
  
end subroutine prim2con_if



