!-*-f90-*-
subroutine output_control
  
  use GR1D_module
  implicit none

  real*8 t_pb, minalp
  integer i

  ! switch back to normal output after this
  ! amount of time
  t_pb = 10.0d-3

  if(initial_data.eq."Collapse") then
     
     if(nt.eq.0) then
        switch1 = .false.
        switch2 = .false.
        switch3 = .false.
        bounce = .false.
     endif

     if(maxval(rho).lt.5.0d12*rho_gf) then

        if(.not.switch1) then
           switch1 = .true.
           ntinfo = ntinfo0
           
           ntout = ntout0
           ntout_scalar = ntout_scalar0
           dtout = dtout0
           dtout_scalar = dtout_scalar0
  
       endif

     else if(maxval(rho).lt.1.0d13*rho_gf) then
        if(.not.switch2) then
           switch2 = .true.
           ntinfo = 10
           ntout = -1
           ntout_scalar = -1
           dtout = .05d-3
           dtout_scalar = 0.01d-3
           ! output now!
           tdump = time
           tdump_scalar = time
        endif
        
     else
        if(.not.bounce) then
           if (eoskey.eq.3) then
              !use entropy definition for bounce
              i = ghosts1+1
              do while(x1(i).lt.3.0d6*length_gf.or.i.gt.n1)
                 if (entropy(i).gt.3.0d0) then
                    bounce = .true.
                    t_bounce = time
                    write(*,"(A20,1P10E15.6)") "bounce! NUCEOS, entropy > 3 in inner core", t_bounce
                    open(unit=666,file=trim(adjustl(outdir))//"/tbounce.dat", &
                         status='unknown',form='formatted',position='rewind')
                    write(666,"(E27.18)") t_bounce
                    close(666)
                    i = n1-ghosts1
                 else
                    i = i+1
                 endif
              enddo
           else
              !use 2e14 as bounce definition
              if(maxval(rho).gt.2.0d14*rho_gf) then
                 bounce = .true.
                 t_bounce = time
                 write(*,"(A20,1P10E15.6)") "bounce! Hybrid, rho > 2.0d14", t_bounce
                 open(unit=666,file=trim(adjustl(outdir))//"/tbounce.dat", &
                      status='unknown',form='formatted',position='rewind')
                 write(666,"(E27.18)") t_bounce
                 close(666)
              endif
           endif
        else
           if(time-t_bounce .gt. t_pb) then
              if(.not.switch3) then
                 switch3 = .true.
                 write(*,*) "Switching to postbounce dtdump"
                 ntinfo = ntinfo0
                 ntout = ntout0
                 ntout_scalar = ntout_scalar0
                 dtout = dtout0
                 dtout_scalar = dtout_scalar0
              endif
           endif
        endif
     endif
  endif

  !min lapse to increase output frequency near BH formation
  if (alp(ghosts1+1).lt.0.4) then
     ntout = 10
     ntout_scalar = 10
     if (alp(ghosts1+1).lt.0.3) then
        ntinfo = 1
        ntout = 1
        ntout_scalar = 1
        if (initial_data.eq."OSC") then
           ntinfo = 10
           ntout = 10
           ntout_scalar = 10
        endif
     endif
  endif

end subroutine output_control
