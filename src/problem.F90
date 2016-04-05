! -*-f90-*-
subroutine problem
     
  use GR1D_module
  use hybrid_eos_module
  use poly_eos_module
  use ideal_eos_module
  implicit none
  
  integer i
  
  real*8 xmin,xmax,mindx
  character(len=100) outfilename
  
  ! which eos do we want to use; eos scalars
  write(*,*) "Setting up EOS: ", eoskey
  if(eoskey.eq.1) then
     gamma1 = hybridgamma1
     gamma2 = hybridgamma2
     gammath = hybridgamma_th
     K1 = 1.2435d15 * (0.5d0**(4.d0/3.d0))
     rhonuc = 2.0d14
     call init_hybrid_eos
  else if (eoskey.eq.2) then
     polyK = 1.2435d15 * (0.5d0**(4.d0/3.d0))
     polygamma = 5.0d0/3.0d0
     !polyK = 1.455e5 !for PNS migration test
     !polygamma = 2.0d0 !for PNS migration test
  else if (eoskey.eq.4) then
     idealgamma = 1.4d0
     idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
  else if(eoskey.eq.3) then
     write(6,*) "Using the NUC_EOS interface."
  else
     stop "choice of eos not implemented"
  endif
  
  if(initial_data.eq."Collapse") then

     write(*,*) "Collapsing a star!"
     call collapse
     
  else if (initial_data .eq. "Shocktube") then
     
     write(*,*) "Shock Tube!"
     call shocktube
     
  else if(initial_data.eq."Sedov") then
     
     write(*,*) "Sedov Blast Wave!"
     call sedov
     
  else if(initial_data.eq."OSC") then

     write(*,*) "OSCing!"
     call OSC

  else if(initial_data.eq."M1test") then

     !setup everything in for the M1test
     write(*,*) "Setting up M1 test case #",M1_testcase_number
     call M1test(M1_testcase_number)

  else 
     
     stop "Please specify initial data!"
     
  endif
  
  call boundaries(0,0)
     
  outfilename = trim(adjustl(outdir))//"/volume.xg"
  call output_single(volume,outfilename)

  ! set up "metric"
  if(geometry.eq.2) then
     sqrt_gamma(:) = X(:)*x1(:)**2
  else
     sqrt_gamma(:) = 1.0d0
  endif
  
  
end subroutine problem
