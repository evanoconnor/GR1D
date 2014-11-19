!-*-f90-*-
subroutine neutrino_pressure
  
  use GR1D_module
  implicit none

  !4/3*pi*(1MeV)^4/(hc)^3
  real*8 :: press_nu_constant = 3.52127727d24*press_gf
  real*8 :: pito4 = 97.4090910340024
  real*8 :: pito2 = 9.86960440108936
  real*8 eta,FD  
  real*8 h1,h2,s1,s2
  integer i

  !this calculates the neutrino pressure in each cell and also the pressure gradiant (to second order)
  !From Christian's CoConut Microphysics talk
  
  do i=1,n1
     if (rho(i).gt.2.0d12*rho_gf) then
        eta = nuchem(i)/temp(i)
        FD = 7.0d0*pito4/60.0d0+0.5d0*(eta**2)*(pito2+0.5d0*(eta**2))
        press_nu(i) = press_nu_constant*(temp(i)**4)*FD
        energy_nu(i) = 0.0d00 !because the energy is still in the matter
     else
        press_nu(i) = 0.0d0
        energy_nu(i) = 0.0d00 !because the energy is still in the matter
     endif
  enddo

  dnupdr(:) = 0.0d0

  do i=1+1,n1-1
     h1=x1(i)-x1(i-1)
     h2=x1(i+1)-x1(i)
     s1=(press_nu(i)-press_nu(i-1))/h1
     s2=(press_nu(i+1)-press_nu(i))/h2

     dnupdr(i) = (s1*h2+s2*h1)/(h1+h2)
     if (dnupdr(i).ne.0.0d0) then
     ! to prevent glitches in dnupdr near rho=2.0d12
        if (press_nu(i+1).eq.0.0d0) then
           dnupdr(i) = 0.5d0*dnupdr(i-1)
        else if (press_nu(i).eq.0.0d0) then
           dnupdr(i) = 0.5d0*dnupdr(i-1)
        else if (press_nu(i-1).eq.0.0d0) then
           dnupdr(i) = 0.5d0*dnupdr(i-1)
        endif
     endif
  enddo
  

end subroutine neutrino_pressure


subroutine nu_press_sources

  use GR1D_module
  implicit none

  integer i

  !do nothing in planar geometry
  if(geometry.eq.1) then
     return
  else
     if (GR) then
        do i=2,n1-1
           presssource(i,2) = presssource(i,2) - alp(i)*W(i)*dnupdr(i)
           presssource(i,3) = presssource(i,3) - alp(i)*v(i)*W(i)*dnupdr(i) 
        enddo
     else
        do i=2,n1-1
           presssource(i,2) = presssource(i,2) - dnupdr(i)*sqrt_gamma(i)
           presssource(i,3) = presssource(i,3) - v1(i)*dnupdr(i)*sqrt_gamma(i)
        enddo
     endif
  endif

end subroutine nu_press_sources
