!-*-f90-*-
module timers

  real*8 :: t_c2p, t_fdhlle

contains
  subroutine start_timers
    implicit none

    t_c2p = 0.0d0
    t_fdhlle = 0.0d0

  end subroutine start_timers

  subroutine output_timers
    use GR1D_module
    implicit none
    open(666,file=trim(adjustl(outdir))//"/timers.dat",status="unknown",position="append")
    write(666,"(i8,1P10E15.6)") nt, time, t_c2p, t_fdhlle
    close(666)

  end subroutine output_timers

end module timers
