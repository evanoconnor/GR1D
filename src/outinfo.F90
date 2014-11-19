!-*-f90-*-
subroutine outinfo

  use GR1D_module
  implicit none

  if(eoskey.eq.3) then
     if(GR) then
        write(*,"(i8,1P10E15.6)") nt,time,dt/time_gf,&
             maxval(rho/rho_gf),minval(ye),alp(ghosts1+1)
     else
        write(*,"(i8,1P10E15.6)") nt,time,dt/time_gf,&
             maxval(rho/rho_gf),minval(ye)
     endif
  else
     if(GR) then
        write(*,"(i8,1P10E15.6)") nt,time,dt/time_gf,&
             rho(ghosts1+1)/rho_gf,alp(ghosts1+1)
     else
        write(*,"(i8,1P10E15.6)") nt,time,dt/time_gf,&
             rho(ghosts1+1)/rho_gf
     endif
  endif

end subroutine outinfo
