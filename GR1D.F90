!-*-f90-*-
program GR1D	
    	  
  use GR1D_module
  implicit none

  !Welcome to GR1D
  write(*,*) "#################################################"
  write(*,*) "#################################################"
  write(*,*) "########### GR1D SPHERICAL HYDRO v2 #############"
  write(*,*) "######### Now with Neutrino Transport ###########"
  write(*,*) "################# Nov ??, 2014 ##################"
  write(*,*) "#################################################"

  ! Call problem setup and allocate/initialize variables 
  call start
  write(*,*) "Done with initial data :-)"

  write(*,*) "Begin time integration loop:"
  IntegrationLoop: do 

     call SetTimeStep

     call handle_output

!!   Integrate
     call Step(dt)

     call postStep_analysis
     
  enddo IntegrationLoop
      
  write(*,*) "Shutting down!"
  write(*,*) " "

end program GR1D
