!-*-f90-*-
subroutine sedov
  use GR1D_module
  use ideal_eos_module
  implicit none

  integer i

  real*8 rho1,rho2
  real*8 rchange
  real*8 rblast, eblast, gammablast, rhoblast, chi
  real*8 buffer, buffer2
  real*8 xmin,xmax,mindx

  ! initial timestep etc, forcing variables regardless of parameter file
  dt = 1.0d-10*time_gf
  eoskey = 4
  geometry = 2
  gravity_active = .false.
  xmin = 0.0d0

  if (GR) then

     stop "Does not perform well"

     xmax = (9.4605284d17)*length_gf !one light year
     mindx = xmax/10.0d0/real(radial_zones)
     if (gridtype.ne.'unigrid') stop "Work this gridtype in sedov"
     call grid(xmin,xmax,mindx)

     eblast = 1.0d51*energy_gf
     gammablast = 30.0d0
     rhoblast = amu_cgs*rho_gf
     rblast = (17.0d0*eblast/(8.0d0*pi*gammablast**2*rhoblast))**(1.0d0/3.0d0)
     
     idealgamma = 4.0d0/3.0d0
     idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
     
     do i=ghosts1+1,n1-ghosts1-1
        
        if (x1(i).lt.rblast) then
           
           chi = 1.0d0+8.0d0*gammablast**2*(1.0d0-x1(i)/rblast)
           W(i) = sqrt(1.0d0+gammablast**2/(2.0d0*chi))
           v(i) = (1.0d0-W(i)**(-2.0d0))**(0.5d0)
           v1(i) = v(i)
           rho(i) = 2.0d0*rhoblast*gammablast**2*chi**(-7.0d0/4.0d0)/W(i)
           eps(i) = rhoblast*gammablast**2 &
                *(sqrt(2.0d0)*gammablast*chi**(-23.0d0/12.0d0)- &
                2.0d0*chi**(-7.0d0/4.0d0))/W(i)/rho(i)
           call ideal_eos(rho(i),eps(i),press(i),buffer,buffer2,cs2(i),0)
        else
           rho(i) = rhoblast
           eps(i) = 1.0d-4
           v(i) = 0.0d0
           v1(i) = v(i)
           W(i) = 1.0d0
           call ideal_eos(rho(i),eps(i),press(i),buffer,buffer2,cs2(i),0)
           
        endif
     enddo
     
  else 
     !parameters from Tasker et al. (2008)
     idealgamma = 5.0d0/3.0d0
     idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
     rchange = 0.0875d0
     xmax = 10.0d0
     mindx = xmax/10.0d0/real(radial_zones)
     if (gridtype.ne.'unigrid') stop "Work this gridtype in sedov"
     call grid(xmin,xmax,mindx)
     
     do i=ghosts1+1,n1
        
        rho(i) = 1.0d0
        eps(i) = 1.0d-5
        press(i) = eps(i)*rho(i)*(idealgamma-1.0d0)
        
        if(x1i(i+1).le.rchange) then
           eps(i) = 1.0d5/(4.0d0/3.0d0*pi*rchange**3)
           press(i) = eps(i)*rho(i)*(idealgamma-1.0d0)
        else if(x1i(i).lt.rchange) then
           eps(i) = 1.0d5/(4.0d0/3.0d0*pi*rchange**3)*(rchange**3-x1i(i)**3)/(x1i(i+1)**3-x1i(i)**3)
           press(i) = eps(i)*rho(i)*(idealgamma-1.0d0)
        endif
        
        v1(i) = 0.0d0
        v(i) = 0.0d0
        ye(i) = 0.0d0
        
     enddo
     
  endif
  
end subroutine sedov
