!-*-f90-*-
subroutine tvd_reconstruction(nx, ghosts1a, orig, bextp, bextm,tvd_lim)

  use GR1D_module, only : x1,x1i
  implicit none

  integer nx, ghosts1a
  integer i

  character*(*) tvd_lim

  real*8, dimension(nx) :: orig, bextp, bextm
  real*8 dupw, dloc, delta, ratio, hdelta
  real*8 dupwdx, dlocdx, dxl,dxr,dx 
  real*8 minmod

  bextp = 0.0d0
  bextm = 0.0d0

  if(tvd_lim.eq."MC") then
     do i=2, nx-ghosts1a+1
        dupwdx = x1(i) - x1(i-1)
        dlocdx = x1(i+1) - x1(i)
        dxr = x1i(i+1) - x1(i)
        dxl = x1(i) - x1i(i)
        
        dupw = (orig(i) - orig(i-1))/dupwdx
        dloc = (orig(i+1) - orig(i))/dlocdx
        
        ! MC Lim
        if (dupw*dloc < 0.d0) then
           delta=0.d0
        else 
           delta=sign(min(2.0d0*abs(dupw),2.0d0*abs(dloc),&
                0.5d0*(abs(dupw)+abs(dloc))),(dupw+dloc))
        end if
        
        bextm(i) = orig(i) - dxl*delta
        bextp(i) = orig(i) + dxr*delta
     enddo
  else if(tvd_lim.eq."minmod") then
     do i=ghosts1a-1, nx-ghosts1a+1
        dupwdx = x1(i) - x1(i-1)
        dlocdx = x1(i+1) - x1(i)
        dxl = x1i(i+1) - x1(i)
        dxr = x1(i) - x1i(i)
        
        dupw = (orig(i) - orig(i-1))/dupwdx
        dloc = (orig(i+1) - orig(i))/dlocdx
        
        ! MINMOD
        delta = minmod(dupw,dloc)
        
        bextm(i) = orig(i) - dxl*delta
        bextp(i) = orig(i) + dxr*delta
     enddo

  else 
     write(6,*) tvd_lim
     flush(6)
     stop "Limiter not implemented :-/"
  endif
  
end subroutine tvd_reconstruction

real*8 function minmod(a,b)

  implicit none
  real*8 a,b
  
  if( a*b < 0.0d0 ) then
     minmod = 0.0d0
  else if( abs(a) < abs(b) ) then
     minmod = a
  else
     minmod = b
  endif
      
end function minmod
