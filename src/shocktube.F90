!-*-f90-*-
subroutine shocktube
  use GR1D_module
  use ideal_eos_module
  implicit none
  
  integer i
  
  real*8 rho1,rho2,eps1,eps2,press1,press2
  real*8 rchange
  real*8 xmin,xmax,mindx
  logical have_press
  
  ! initial timestep
  dt = 1.0d-7*time_gf
        
  !just to be sure
  eoskey = 4
  geometry = 1      
  gravity_active = .false.

  xmin = 0.0d0
  xmax = 1.0d0
  mindx = xmax/real(radial_zones)
  
  if (gridtype.ne.'unigrid') stop "Work this gridtype in shocktube"
  call grid(xmin,xmax,mindx)

  rchange = 0.5d0
 
  if (shocktube_problem.eq.1) then
     write(*,*) "Shocktube: M&M LLR Problem 1"
     !M&M LLR Problem 1
     have_press = .true.
     rho1 = 10.0d0
     rho2 = 1.0d0
     press1 = 13.33d0
     press2 = 1.0d-6
     idealgamma = 5.0d0/3.0d0

  elseif (shocktube_problem.eq.2) then
     write(*,*) "Shocktube: M&M LLR Problem 2"
     !M&M LLR Problem 2
     have_press = .true.
     rho1 = 1.0d0
     rho2 = 1.0d0
     press1 = 1000.0d0
     press2 = 0.01d0
     idealgamma = 5.0d0/3.0d0

  elseif (shocktube_problem.eq.3) then
     write(*,*) "Shocktube: RST values"
     !RST values
     have_press = .false.
     rho1 = 1.0d0
     rho2 = 0.125d0
     eps1 = 2.5d0
     eps2 = 2.0d0
     idealgamma = 1.4d0

  elseif (shocktube_problem.eq.4) then
     write(*,*) "Shocktube: SCHN values"
     !SCHN values
     have_press = .false.
     rho1=10.0d0
     rho2=1.0d0
     eps1=2.0d0
     eps2=1.0d-6
     idealgamma = 5.0d0/3.0d0

  else
     write(*,*) "Shocktube Problem #", shocktube_problem, " not setup"
     stop
  endif
     
  
  idealK1 = 1.0d0/(press_gf/rho_gf**idealgamma)

  do i=ghosts1,n1
     
     if(x1(i).le.rchange) then
        rho(i) = rho1
        eps(i) = eps1
        press(i) = press1
     else
        rho(i) = rho2
        eps(i) = eps2
        press(i) = press2
     endif
     
     if (have_press) then
        eps(i) = press(i)/(idealgamma-1.0d0)/rho(i)
     else
        press(i) = (idealgamma-1.0d0)*rho(i)*eps(i)
     endif
     v1(i) = 0.0d0
     ye(i) = 0.0d0

  enddo

end subroutine shocktube
