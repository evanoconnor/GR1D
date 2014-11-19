!-*-f90-*-
module atmos

  real*8 atmo_rho_thr, atmo_rho
  
end module atmos

subroutine atmosphere_init

  use GR1D_module
  use atmos

  if(atmo_rho_rel_min.gt.0.0d0.and.atmo_rho_abs_min.lt.1.0d-13) then
     atmo_rho_thr = maxval(rho)*atmo_rho_rel_min*rho_gf
  else
     atmo_rho_thr = atmo_rho_abs_min*rho_gf
  endif
  atmo_rho = atmo_rho_thr * atmo_fac


end subroutine atmosphere_init

