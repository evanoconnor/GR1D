!-*-f90-*-
subroutine M1test

  use GR1D_module
  implicit none
  real*8 xmin,xmax,mindx
  integer i,j

  real*8 kappa,b,rsurface
  real*8, dimension(18), parameter :: myenergies = (/1.0d0, 3.0d0, &
       5.23824d0, 8.00974d0, 11.4415d0, 15.6909d0, 20.9527d0, 27.4681d0, &
       35.5357d0, 45.5254d0, 57.8951d0, 73.2117d0, 92.1775d0, 115.662d0, &
       144.741d0, 180.748d0, 225.334d0, 280.542d0/)

  real*8 :: top(18),bottom(18)
  real*8 :: v_cut,v_width,v_max
  real*8 :: M1_testcase_travelling_pulse,M1_testcase_diffusion_wave
  real*8 :: constant_rho_TC4

  if (M1_testcase_number.eq.1) then

     stop "not a valid test case number"

  else if (M1_testcase_number.eq.2) then
     !dense radiating sphere, abdikamalov test
     !doublecheck on settings
     if (GR) stop "This test is Newtonian!"
     if (n1.ne.100+2*ghosts1) stop "You need 100 zones"

     M1_control_flag = .false.
     
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     include_energycoupling_imp = .false.
     M1_do_backwardfix = 1

     ! initial timestep
     dt = 1.0d-8*time_gf

     M1_maxradii = 50.0d5
     
     xmin = 0.0d0
     xmax = 50.0d5*length_gf
     
     do_hydro = .false.

     mindx = xmax/real(n1-2*ghosts1)
     gridtype = 'unigrid'
     
     call grid(xmin,xmax,mindx)

     kappa = 250.0d-6
     b = 10.0d0
     rsurface = 10.0d5

     if (x1(1)/length_gf.gt.rsurface) then
        stop "you don't have a surface!"
     endif

     do i=1,n1
        if (x1(i)/length_gf.lt.rsurface) then
           !emissivities and opacities
           eas(i,:,:,1) = kappa*nulib_opacity_gf*b*nulib_emissivity_gf
           eas(i,:,:,2) = kappa*nulib_opacity_gf
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
              
           !initial conditions
           q_M1(i,:,:,1) = eas(i,1,1,1)/eas(i,1,1,2)
           q_M1(i,:,:,2) = 1.0d-10*eas(i,1,1,1)/eas(i,1,1,2)
        else
           !emissivities and opacities
           eas(i,:,:,1) = 0.0d0
           eas(i,:,:,2) = 0.0d0
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
              
           !initial conditions
           q_M1(i,:,:,1) = eas(1,1,1,1)/eas(1,1,1,2)*(rsurface*length_gf/x1(i))**2 !falls as 1/r^2 from the surface
           q_M1(i,:,:,2) = 0.1*eas(1,1,1,1)/eas(1,1,1,2)*(rsurface*length_gf/x1(i))**2 !half the energy density
           
        endif
     enddo


  else if (M1_testcase_number.eq.3) then
     !travelling pulses
     !doublecheck on settings
     if (GR) stop "This test is Newtonian!"
     if (n1.ne.200+2*ghosts1) stop "You need 200 zones"

     M1_control_flag = .false.
     
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     include_energycoupling_imp = .false.
     M1_do_backwardfix = 1

     ! initial timestep
     dt = 1.0d-8*time_gf
     dtout = 0.5d0/time_gf
     dtout_scalar = 0.5d0/time_gf
     tend = 6.0d0/time_gf

     M1_maxradii = 10.2/length_gf
     
     xmin = 0.2
     xmax = 10.2
     
     do_hydro = .false.

     mindx = (xmax-xmin)/200.0
     gridtype = 'unigrid'
     
     call grid(xmin,xmax,mindx)

     do i=1,n1
        !emissivities and opacities
        eas(i,:,:,1) = 0.0d0
        eas(i,:,:,2) = 1.0d-40*nulib_opacity_gf
        eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
        
        !initial conditions, availabe from M1_testcase_travelling_pulse(radius,time,0/1) routine
        q_M1(i,:,:,1) = M1_testcase_travelling_pulse(x1(i),0.0d0,0)
        q_M1(i,:,:,2) = M1_testcase_travelling_pulse(x1(i),0.0d0,1)
        q_M1(i,:,:,3) = 1.0d0

     enddo

  else if (M1_testcase_number.eq.4) then
     !Smit test
     !not so dense radiating sphere
     if (GR) stop "This test is for NW only!"
     if (n1.ne.800+2*ghosts1) stop "You need 800 zones"

     M1_control_flag = .false.
     
     M1_do_backwardfix = 1
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     include_energycoupling_imp = .false.

     ! initial timestep
     dt = 1.0d-7*time_gf

     M1_maxradii = 3.0d6
     xmin = 0.0d0
     xmax = 3.0d6*length_gf
     
     do_hydro = .false.
     
     mindx = xmax/real(n1-2*ghosts1)
     gridtype = 'unigrid'
     
     call grid(xmin,xmax,mindx)

     kappa = 4.0e-6
     b = 0.8d0
     rsurface = 1.0d6

     if (x1(1)/length_gf.gt.rsurface) then
        stop "you don't have a surface!"
     endif

     bottom(1) = 0.0d0
     top(1) = myenergies(1)*2.0d0
     do j=2,18
        bottom(j) = top(j-1)
        top(j) = (myenergies(j)-bottom(j))*2.0d0+bottom(j)
     enddo

     do i=1,n1
        if (x1(i)/length_gf.lt.rsurface) then
           !emissivities and opacities
           eas(i,:,:,1) = kappa*nulib_opacity_gf*b*nulib_emissivity_gf
           eas(i,:,:,2) = kappa*nulib_opacity_gf
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
              
           !initial conditions
           q_M1(i,:,:,1) = eas(i,1,1,1)/eas(i,1,1,2)
           q_M1(i,:,:,2) = 1.0d-10*eas(i,1,1,1)/eas(i,1,1,2)
        else
           !emissivities and opacities
           eas(i,:,:,1) = 0.0d0
           eas(i,:,:,2) = 0.0d0
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
              
           !initial conditions
           q_M1(i,:,:,1) = eas(1,1,1,1)/eas(1,1,1,2)*(rsurface*length_gf/x1(i))**2 !falls as 1/r^2 from the surface
           q_M1(i,:,:,2) = 0.1*eas(1,1,1,1)/eas(1,1,1,2)*(rsurface*length_gf/x1(i))**2 !half the energy density
           
        endif
     enddo

  else if (M1_testcase_number.ge.5.and.M1_testcase_number.le.7) then
     !Mueller et al. (2010)
     !Three tests, all static hydro
     !Test #5: background velocity field with shock
     !Test #6: background gravitational field
     !Test #7: combined #5 and #6

     M1_control_flag = .false.
     
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     v_order = -1
     include_energycoupling_imp = .true.
     number_species_to_evolve = 1
     M1_do_backwardfix = 0

     if (M1_testcase_number.eq.5) then
        constant_rho_TC4 = 1.0d8
        v_cut = 135.0d5
        v_width = 15.0d5 
        v_max = -0.2d0
     elseif (M1_testcase_number.eq.6) then
        constant_rho_TC4 = 9.0d14
        v_cut = 135.0d5
        v_width = 15.0d5
        v_max = 0.0d0        
     else
        constant_rho_TC4 = 9.0d14
        v_cut = 135.0d5
        v_width = 15.0d5
        v_max = -0.2d0
     endif

     if (.not.GR) stop "This test is for GR only!"
     if (n1.ne.300+2*ghosts1) stop "You need 300 zones"

     ! initial timestep
     dt = 1.0d-8*time_gf

     M1_maxradii = 1.0d9
     xmin = 0.0d0
     xmax = 1.0d9*length_gf
     
     kappa = 4.0d-4
     b = 0.8d0
     rsurface = 10.0d5

     do_hydro = .false.

     mindx = 1.0d4*length_gf
     gridtype = 'log'

     call grid(xmin,xmax,mindx)

     if (x1(1)/length_gf.gt.rsurface) then
        stop "you don't have a surface!"
     endif

     bottom(1) = 0.0d0
     top(1) = myenergies(1)*2.0d0
     do j=2,18
        bottom(j) = top(j-1)
        top(j) = (myenergies(j)-bottom(j))*2.0d0+bottom(j)
     enddo

     do i=1,n1-1
        if (x1(i)/length_gf.lt.v_cut) then
           v(i) = 0.0d0
        elseif(x1(i)/length_gf.lt.v_cut+v_width) then
           v(i) = (x1(i)/length_gf-v_cut)/v_width*v_max
        else 
           v(i) = v_max*((v_cut+v_width)/(x1(i)/length_gf))**2
        endif
        if (x1i(i)/length_gf.lt.v_cut) then
           vm(i) = 0.0d0
        elseif(x1i(i)/length_gf.lt.v_cut+v_width) then
           vm(i) = (x1i(i)/length_gf-v_cut)/v_width*v_max
        else 
           vm(i) = v_max*((v_cut+v_width)/(x1i(i)/length_gf))**2
        endif
        if (x1i(i+1)/length_gf.lt.v_cut) then
           vp(i) = 0.0d0
        elseif(x1i(i+1)/length_gf.lt.v_cut+v_width) then
           vp(i) = (x1i(i+1)/length_gf-v_cut)/v_width*v_max
        else 
           vp(i) = v_max*((v_cut+v_width)/(x1i(i+1)/length_gf))**2
        endif
        v_prev(i) = v(i)
        W(i) = 1.0d0/(1.0d0-v(i)**2)**0.5d0
        Wm(i) = 1.0d0/(1.0d0-vm(i)**2)**0.5d0
        Wp(i) = 1.0d0/(1.0d0-vp(i)**2)**0.5d0
        
        ye(i) = 0.5d0
        temp(i) = 1.0d0

        if (x1(i)/length_gf.lt.rsurface) then
           rho(i) = constant_rho_TC4*rho_gf
           do j=1,18
              eas(i,:,j,1) = kappa*nulib_opacity_gf*b*nulib_emissivity_gf* &
                   (myenergies(j)/5.0d0)**3/(exp(myenergies(j)/5.0d0)+1.0d0)*(top(j)-bottom(j))
              eas(i,:,j,2) = kappa*nulib_opacity_gf
              
              !initial conditions
              q_M1(i,:,j,1) = eas(i,:,j,1)/eas(i,:,j,2)
              q_M1(i,:,j,2) = q_M1(i,:,j,1)/1.0d2
           enddo
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf              
        else
           !emissivities and opacities
           rho(i) = 1.0d-10*rho_gf
           eas(i,:,:,1) = 0.0d0
           eas(i,:,:,2) = 0.0d0
           eas(i,:,:,3) = 1.0d-40*nulib_opacity_gf
              
           q_M1(i,:,:,1) = eas(1,:,:,1)/eas(1,:,:,2)* &
                (rsurface*length_gf/x1(i))**2 !falls as 1/r^2 from the surface
           q_M1(i,:,:,2) = 0.7d0*eas(1,:,:,1)/eas(1,:,:,2)* &
                (rsurface*length_gf/x1(i))**2 !half the energy density
        endif
     enddo
     
     call prim2con
     call con2GR
     call GR_alp

  else if (M1_testcase_number.eq.8) then
     !weak diffusion wave
     if (GR) stop "This test is not for GR!"
     if (n1.ne.100+2*ghosts1) stop "You need 100 zones"

     M1_control_flag = .false.
     
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     include_energycoupling_imp = .false.
     number_species_to_evolve = 1

     ! initial timestep
     dt = 1.0d-12*time_gf

     M1_maxradii = 1.0d0/length_gf
     xmin = 0.0d0
     xmax = 1.0d0
     
     dtout = 1.d0/time_gf
     dtout_scalar = 1.0d0/time_gf
     tend = 5.0d0/time_gf

     time = 0.0d0

     kappa = 100.0d0/nulib_opacity_gf
     b = 0.0d0

     do_hydro = .false.

     mindx = xmax/(n1-2*ghosts1)
     gridtype = 'unigrid'
     
     call grid(xmin,xmax,mindx)

     do i=1,n1
        do j=1,number_groups
           eas(i,:,j,1) = 1.0d-80*nulib_opacity_gf
           eas(i,:,j,2) = 1.0d-80*nulib_opacity_gf
           eas(i,:,j,3) = kappa*nulib_opacity_gf

           !initial conditions
           q_M1(i,1,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),0,M1_testcase_number)
           q_M1(i,1,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),1,M1_testcase_number)
           q_M1(i,1,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),0,M1_testcase_number)/3.0d0
           q_M1(i,2,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),0,M1_testcase_number)
           q_M1(i,2,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),1,M1_testcase_number)
           q_M1(i,2,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),0,M1_testcase_number)/3.0d0
           q_M1(i,3,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),0,M1_testcase_number)
           q_M1(i,3,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),1,M1_testcase_number)
           q_M1(i,3,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),0,M1_testcase_number)/3.0d0
        enddo
     enddo
     q_M1_fluid = q_M1

  else if (M1_testcase_number.eq.9) then
     !heavy diffusion wave
     if (GR) stop "This test is not for GR!"
     if (n1.ne.100+2*ghosts1) stop "You need 100 zones"

     M1_control_flag = .false.
     
     include_epannihil_kernels = .false.
     include_Ielectron_exp = .false.
     include_Ielectron_imp = .false.
     include_Ielectron = .false.
     include_energycoupling_exp = .false.
     include_energycoupling_imp = .false.
     number_species_to_evolve = 1

     ! initial timestep
     dt = 1.0d-12*time_gf

     M1_maxradii = 1.0d0/length_gf
     xmin = 0.0d0
     xmax = 1.0d0
     
     dtout = 20.0d0/time_gf
     dtout_scalar = 20.0d0/time_gf
     tend = 200.0d0/time_gf

     time = 0.0d0

     kappa = 100000.0d0/nulib_opacity_gf
     b = 0.0d0

     do_hydro = .false.

     mindx = xmax/(n1-2*ghosts1)
     gridtype = 'unigrid'
     
     call grid(xmin,xmax,mindx)

     do i=1,n1
        do j=1,number_groups
           eas(i,:,j,1) = 1.0d-80*nulib_opacity_gf
           eas(i,:,j,2) = 1.0d-80*nulib_opacity_gf
           eas(i,:,j,3) = kappa*nulib_opacity_gf

           !initial conditions
           q_M1(i,1,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),0,M1_testcase_number)
           q_M1(i,1,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),1,M1_testcase_number)
           q_M1(i,1,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,1,j,3),0,M1_testcase_number)/3.0d0
           q_M1(i,2,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),0,M1_testcase_number)
           q_M1(i,2,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),1,M1_testcase_number)
           q_M1(i,2,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,2,j,3),0,M1_testcase_number)/3.0d0
           q_M1(i,3,j,1) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),0,M1_testcase_number)
           q_M1(i,3,j,2) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),1,M1_testcase_number)
           q_M1(i,3,j,3) = M1_testcase_diffusion_wave(x1(i),0.0d0,eas(i,3,j,3),0,M1_testcase_number)/3.0d0
        enddo
     enddo
     q_M1_fluid = q_M1

  else
     stop "Add in this test case"
  endif

end subroutine M1test

function M1_testcase_travelling_pulse(radius,time,moment) result(moment_value)

  implicit none

  real*8 radius,time
  integer moment
  real*8 moment_value

  if (moment.eq.0) then
     moment_value = exp(-(radius-time)**2)/radius**2 + exp(-((radius+time)-9.0d0)**2)/radius**2
  else if (moment.eq.1) then
     moment_value = 0.99999999d0*(exp(-(radius-time)**2)/radius**2 - exp(-((radius+time)-9.0d0)**2)/radius**2)
  else
     stop "request correct moment"
  endif

end function M1_testcase_travelling_pulse

function M1_testcase_diffusion_wave(radius,time,kappa,moment,TCnumber) result(moment_value)

  implicit none

  real*8 radius,time,kappa
  integer moment,TCnumber
  real*8 moment_value
  real*8 timeshift

  if (TCnumber.eq.8) then
     timeshift = 1.0d0
  else if (TCnumber.eq.9) then
     timeshift = 200.0d0
  else
     stop "add in TCnumber"
  endif

  if (moment.eq.0) then
     moment_value = (kappa/(time+timeshift))**(3.0d0/2.0d0)* &
          exp(-3.0d0*kappa*radius**2/(4.0d0*(time+timeshift)))
  else if (moment.eq.1) then
     moment_value = (kappa/(time+timeshift))**(3.0d0/2.0d0)* &
          exp(-3.0d0*kappa*radius**2/(4.0d0*(time+timeshift)))*radius/2.0d0/(time+timeshift)
  else
     stop "request correct moment"
  endif

end function M1_testcase_diffusion_wave
