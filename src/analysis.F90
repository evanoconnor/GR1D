!-*-f90-*-
subroutine findaccretionradii

  use GR1D_module
  implicit none

  integer i
  real*8 interiormass(n1)

  !have 11 spots for accretion radii
  !shock,10,15,20,25,50,100,150,200,250,300 km

  accretion_radii(:) = ghosts1+1

  if (GR) then
     interiormass = mgravi
  else
     interiormass = mass
  endif

  accretion_radii(1) = max(ishock(1),ghosts1+1)
  accreted_mass(1) = interiormass(ishock(1))

  accretion_analysis = .true.
  
  if (x1i(n1-ghosts1)/length_gf.lt.300.0d5) then
     write(*,*) "Change accretion radii if accretion rates required."
     write(*,*) "turning accretion analysis off..."
     accretion_analysis = .false.
  endif

  do i=ghosts1+1,n1-ghosts1
     if (accretion_analysis) then
        if (x1i(i).lt.10.0d5*length_gf) then
           accretion_radii(2) = i
           accreted_mass(2) = interiormass(i)
        elseif (x1i(i).lt.15.0d5*length_gf) then
           accretion_radii(3) = i
           accreted_mass(3) = interiormass(i)
        elseif (x1i(i).lt.20.0d5*length_gf) then
           accretion_radii(4) = i
           accreted_mass(4) = interiormass(i)
        elseif (x1i(i).lt.25.0d5*length_gf) then
           accretion_radii(5) = i
           accreted_mass(5) = interiormass(i)
        elseif (x1i(i).lt.50.0d5*length_gf) then
           accretion_radii(6) = i
           accreted_mass(6) = interiormass(i)
        elseif (x1i(i).lt.100.0d5*length_gf) then
           accretion_radii(7) = i
           accreted_mass(7) = interiormass(i)
        elseif (x1i(i).lt.150.0d5*length_gf) then
           accretion_radii(8) = i
           accreted_mass(8) = interiormass(i)
        elseif (x1i(i).lt.200.0d5*length_gf) then
           accretion_radii(9) = i
           accreted_mass(9) = interiormass(i)
        elseif (x1i(i).lt.250.0d5*length_gf) then
           accretion_radii(10) = i
           accreted_mass(10) = interiormass(i)
        elseif (x1i(i).lt.300.0d5*length_gf) then
           accretion_radii(11) = i     
           accreted_mass(11) = interiormass(i)   
        endif
     endif
  enddo

  accreted_mass(1:11) = 0.0d0

  if (accretion_radii(11).eq.0) then
     stop "Missing outer accretion radii, decrease shell radius"
  endif

  !Write out headers
  open(unit=666,file=trim(adjustl(outdir))//trim(adjustl("/accretion_rates.dat")),status= &
       "unknown", form='formatted',position="append")
  write(666,"('#Radii: ',1P20E18.9)") x1(accretion_radii(1))/length_gf, &
        x1(accretion_radii(2))/length_gf, x1(accretion_radii(3))/length_gf, &
        x1(accretion_radii(4))/length_gf, x1(accretion_radii(5))/length_gf, &
        x1(accretion_radii(6))/length_gf, x1(accretion_radii(7))/length_gf, &
        x1(accretion_radii(8))/length_gf, x1(accretion_radii(9))/length_gf, &
        x1(accretion_radii(10))/length_gf, x1(accretion_radii(11))/length_gf
  close(666)
  open(unit=666,file=trim(adjustl(outdir))//trim(adjustl("/accreted_mass.dat")),status= &
       "unknown", form='formatted',position="append")
  write(666,"('#Radii: ',1P20E18.9)") x1(accretion_radii(1))/length_gf, &
        x1(accretion_radii(2))/length_gf, x1(accretion_radii(3))/length_gf, &
        x1(accretion_radii(4))/length_gf, x1(accretion_radii(5))/length_gf, &
        x1(accretion_radii(6))/length_gf, x1(accretion_radii(7))/length_gf, &
        x1(accretion_radii(8))/length_gf, x1(accretion_radii(9))/length_gf, &
        x1(accretion_radii(10))/length_gf, x1(accretion_radii(11))/length_gf
   close(666)


end subroutine findaccretionradii

subroutine mass_analysis

  use GR1D_module
  implicit none

  integer i
  real*8 maxX
  integer maxXi,mass12i
  logical contflag
  integer inner_core_i

  totalmass = 0.0d0

  if (GR) then
     totalmass = mgravi(n1-ghosts1) 
     !accretion rates in solar masses/s
     if (gravity_active.and.accretion_analysis) then
        accretion_rates(1) = -time_gf*4.0d0*pi*shock_radius**2 &
             * rho(ishock(1))*v(ishock(1))
        accreted_mass(1) = mgravi(ishock(1))

        do i=2,11 
           if(accretion_radii(i).gt.n1-ghosts1.or. &
                accretion_radii(i).lt.ghosts1+1) then
              accretion_rates(i) = 0.0d0
              cycle
           endif
           accretion_rates(i) = -time_gf*4.0d0*pi * &
                x1(accretion_radii(i))**2 * &
                rho(accretion_radii(i))*v(accretion_radii(i))
           accreted_mass(i) = mgravi(accretion_radii(i))
        enddo
     endif
  else 
     totalmass = mass(n1-ghosts1)
     if (gravity_active.and.accretion_analysis) then
        accretion_rates(1) = -time_gf*4.0d0*pi*shock_radius**2 &
             * rho(ishock(1))*v1(ishock(1))
        accreted_mass(1) = mass(ishock(1))

        do i=2,11
           if(accretion_radii(i).gt.n1-ghosts1.or. &
                accretion_radii(i).lt.ghosts1+1) then
              accretion_rates(i) = 0.0d0
              cycle
           endif 
           accretion_rates(i) = -time_gf*4.0d0*pi * &
                x1(accretion_radii(i))**2 * &
                rho(accretion_radii(i))*v1(accretion_radii(i))
           accreted_mass(i) = mass(accretion_radii(i))
        enddo
     endif
  endif

  !also want various masses
  maxX = 1.0d0
  maxXi = ghosts1+1 
  mass12i = 0
  contflag = .true.
  do i=ghosts1+1,n1
     if (X(i).gt.maxX) then
        maxX = X(i)
        maxXi = i
     endif
     
     if (rho(i)/rho_gf.gt.1.0d12.and.contflag) then
        mass12i = i
     else
        contflag = .false.
     endif

  enddo

  if(maxXi.gt.ghosts1.and.maxXi.le.n1-ghosts1) then
     mgravX = mgrav(maxXi)
     mbaryX = mass(maxXi)
     rXmax = x1(maxXi)
  else
     mgravX = 0.0d0
     mbaryX = 0.0d0
     rXmax = 0.0d0
  endif

  if(mass12i.gt.ghosts1.and.mass12i.le.n1-ghosts1) then
     mgrav12 = mgrav(mass12i)
     mbary12 = mass(mass12i)
     r12max = x1(mass12i)
  else
     mgrav12 = 0.0d0
     mbary12 = 0.0d0
     r12max = 0.0d0
  endif

  !now find inner core size
  !first define core to be r where v>cs.  If not the case then must be early on in evolution, just take r=r(min(v1))

  inner_core_i = ghosts1+1
  contflag = .true.
  do i=ghosts1+1,n1-ghosts1-1
     if (contflag) then
        if (abs(v1(i)).lt.sqrt(cs2(i))) then
           inner_core_i = i
        else
           contflag = .false.
        endif
     endif
  enddo

  if (inner_core_i.gt.n1-ghosts1-2) then
     inner_core_i = ishock(1)
  endif

  mass_inner_core = mass(inner_core_i)

end subroutine mass_analysis

subroutine get_shock_radius

  use GR1D_module, only : v1,ishock,shock_radius,x1,bounce, &
       n1,ghosts1
  use atmos

  implicit none

  integer i

  if (.not.bounce) then
     ishock(1) = ghosts1+1
     shock_radius = 0.0d0
     return
  endif

  ishock = minloc(v1( (ghosts1+1) : (n1-ghosts1-1) ))
  shock_radius = x1(ishock(1))

end subroutine get_shock_radius

subroutine rotation_analysis
  
  use GR1D_module
  implicit none 

  integer i

  if(GR) then
     omega(:) = vphi(:)/x1(:)
     angular_momentum = 0.0d0
     rotational_energy = 0.0d0
     internal_energy = 0.0d0
     do i=ghosts1+1,n1-ghosts1
        angular_momentum = angular_momentum + X(i)*W(i)**2*rho(i)*volume(i)*vphi(i)* &
             x1(i)*(1.0d0+eps(i)+press(i)/rho(i))*twothirds
        rotational_energy = rotational_energy + 0.5d0*X(i)*W(i)**2*rho(i)* &
             volume(i)*vphi(i)**2*(1.0d0+eps(i)+press(i)/rho(i))*twothirds
        internal_energy = internal_energy + X(i)*rho(i)*eps(i)*W(i)*volume(i)
        ToverW(i) = rotational_energy / abs(mgrav(i)- &
       mass(i)-internal_energy-rotational_energy)
     enddo
  else
     omega(:) = vphi1(:)/x1(:)
     angular_momentum = 0.0d0
     rotational_energy = 0.0d0
     internal_energy = 0.0d0
     do i=ghosts1+1,n1-ghosts1
        angular_momentum = angular_momentum + rho(i)*volume(i)*vphi1(i)*x1(i)*twothirds
        rotational_energy = rotational_energy + 0.5d0*rho(i)*volume(i)*vphi1(i)**2*twothirds
        !internal energy is playing the role of gravitational energy here
        internal_energy = internal_energy - rho(i)*volume(i)*mass(i)/x1(i)
        ToverW(i) = rotational_energy/abs(internal_energy)
     enddo
  endif

end subroutine rotation_analysis

subroutine get_binding_energy

  use GR1D_module
  implicit none

  integer i

  binding_energy_total = 0.0d0
  !rough binding energy
  do i=ghosts1+1,n1-ghosts1
     if (x1(i).gt.shock_radius) then
        binding_energy(i) = mgravi(i)*rho(i)*volume(i)/x1(i)
        binding_energy_total = binding_energy_total + binding_energy(i)
     else
        binding_energy(i) = 0.0d0
     endif
  enddo
  
  binding_energy_total = binding_energy_total + binding_energy_envelope

end subroutine get_binding_energy

subroutine map_envelope_binding_energy(lprofile_name)
  
  !this is the binding energy of the envelope not included in our grid

  use GR1D_module, only : mgrav, n1,ghosts1,grid_rmax,binding_energy_envelope, &
       energy_gf,mass_gf,rho_cut,rho_gf,rho,length_gf
  implicit none
  
  character*(*) lprofile_name
  integer profile_zones

  real*8 buffer, dmass, dx
  integer i,ibuffer
  real*8 radius_cut
  integer rho_cut_i
  
  integer keytemp,keyerr,eosflag
  real*8 binding_energy_total
  real*8 binding_energy_shell

  real*8, allocatable :: pradius(:), &
       pmass(:),prho(:),ptemp(:), &
       ppress(:),peps(:),pvel(:),&
       pye(:),pomega(:)


  real*8 :: kboltz_cgs = 1.380662d-16
  real*8 :: mgrav_interior

! read profile      
  open(666,file=trim(lprofile_name),status='unknown', & 
       form='formatted',action='read')
  read(666,*) profile_zones
  
  allocate(pradius(profile_zones),pmass(profile_zones))
  allocate(prho(profile_zones),ppress(profile_zones))
  allocate(ptemp(profile_zones))
  allocate(peps(profile_zones))
  allocate(pvel(profile_zones))
  allocate(pye(profile_zones))
  allocate(pomega(profile_zones))
  
  do i=1,profile_zones
     read(666,*) ibuffer,pmass(i),pradius(i),&
          ptemp(i),prho(i),pvel(i),pye(i), &
          buffer
  enddo
  
  close(666)
  
  do i=1,profile_zones
     if (prho(i).gt.rho(n1-ghosts1-1)/rho_gf) then
        rho_cut_i = i
     endif
  enddo
      
  mgrav_interior = mgrav(n1-ghosts1-1)
  binding_energy_total = 0.0d0

  !convert to solar mass units 
  if (pmass(2).gt.1.0d10) then
     pmass(:) = pmass(:)*mass_gf
  endif

  do i=rho_cut_i,profile_zones
     binding_energy_shell = (pmass(i)-pmass(i-1))*mgrav_interior/ &
          (pradius(i)*length_gf)
     mgrav_interior = mgrav_interior + binding_energy_shell + & 
          (pmass(i)-pmass(i-1))
     binding_energy_total = binding_energy_total + binding_energy_shell
  enddo

  binding_energy_envelope = binding_energy_total

  write(*,*) "Binding energy of envelope (ergs):", binding_energy_envelope/energy_gf

  deallocate(pradius,pmass)
  deallocate(prho,ppress)
  deallocate(ptemp)
  deallocate(peps)
  deallocate(pvel)
  deallocate(pye)
  deallocate(pomega)

end subroutine map_envelope_binding_energy
  
subroutine analytic_OSC_alpha(time,radius_o,M,alpha,rho,vel,X,maxr)

  use GR1D_module, only : pi,x1,n1,ghosts1
  implicit none
  
  real*8 radius_guess, radius_o
  real*8 time, time_guess
  real*8 M
  real*8 alpha(n1),rho(n1),W(n1),e(n1),vel(n1),X(n1)
  real*8 maxr

  real*8 chi_s
  real*8 chi
  real*8 eta_c
  real*8 eta_c_star ! maximum
  real*8 oldeta_low,oldchi_low
  real*8 oldeta_high,oldchi_high
  real*8 tol
  real*8 eta_guess,chi_guess
  real*8 radius_t
  
  integer i,gi,counter
  integer index
  logical cont

  tol = 1.0d-8

  !this routine finds the lapse for a given r,t,M,R(t=0)
  chi_s = asin(sqrt(2.0d0*M/radius_o))

  ! first eta_c
  ! note eta_c increases linear with time therefore if eta_c guess
  ! gives time too low, choose higher time.
  eta_c_star = -2.0d0*acos((1.0d0-2.0d0*M/radius_o)**(0.75d0))

  cont = .true.
  i = 0
  oldeta_low = -pi
  oldeta_high = eta_c_star
 
  if (time.eq.0.0d0) then
     cont=.false.
     eta_c = -pi
  endif

  do while(cont)
     i = i+1
     eta_guess = (oldeta_low+oldeta_high)/2.0d0
     call get_OSC_ana_time(eta_guess,M,chi_s,time_guess)
     if (time_guess.gt.time) then
        oldeta_high = eta_guess
     else
        oldeta_low = eta_guess
     endif
     
     if (abs(time-time_guess)/time.lt.tol) then
        cont = .false.
        eta_c = eta_guess
     endif
  enddo

  radius_t = radius_o*(1.0d0-(cos(eta_c/2.0d0))**2/cos(chi_s))
  maxr = radius_t

  do counter=ghosts1+1,n1-ghosts1
     if (x1(counter).lt.radius_t) then
        !have eta_c, now get chi
        cont = .true.
        i = 0
        oldchi_low = 0.0d0
        oldchi_high = chi_s
        do while(cont)
           i = i+1
           chi_guess = (oldchi_low+oldchi_high)/2.0d0
           call get_OSC_ana_r(chi_guess,radius_o,eta_c,chi_s,radius_guess)
           if (radius_guess.gt.x1(counter)) then
              oldchi_high = chi_guess
           else
              oldchi_low = chi_guess
           endif
           
           if (abs(x1(counter)-radius_guess)/x1(counter).lt.tol) then
              cont = .false.
              chi = chi_guess
           endif
        enddo
        
        alpha(counter) = (cos(chi)-(cos(eta_c/2.0d0))**2)*((cos(chi_s))**3 - & 
             (cos(eta_c/2.0d0))**2)/(((cos(chi))**3- &
             (cos(eta_c/2.0d0))**2)**(0.5d0)*(cos(chi_s)- &
             (cos(eta_c/2.0d0))**2)**(1.5d0))
        W(counter) = cos(chi)*((cos(chi)-(cos(eta_c/2.0d0))**2)/ &
             ((cos(chi))**3-(cos(eta_c/2.0d0))**2))**(0.5d0)
        X(counter) = sqrt((cos(chi)-(cos(eta_c/2.0d0))**2)/ &
             ((cos(chi))**3-(cos(eta_c/2.0d0))**2))
        e(counter) = 3.0d0*(sin(chi_s))**6/(32.0d0*pi*M**2)*(cos(chi))**3/ &
             (cos(chi)-(cos(eta_c/2.0d0))**2)**3
        rho(counter) = e(counter)
        vel(counter) = -cos(eta_c/2.0d0)*tan(chi)/sqrt(cos(chi)-(cos(eta_c/2.0d0))**2)
        

     else
        alpha(counter) = sqrt(1.0d0-2.0d0*M/x1(counter))
        X(counter) = 1.0d0/sqrt(1.0d0-2.0d0*M/x1(counter))
        rho(counter) = 0.0d0
        vel(counter) = 0.0d0
     endif

  enddo

  ! boundaries
  gi = ghosts1
  do counter=ghosts1,1,-1
     gi = gi + 1
     alpha(counter) = alpha(gi)
     W(counter) = W(gi)
     X(counter) = X(gi)
     rho(counter) = rho(gi)
     vel(counter) = vel(gi)
  enddo

  gi = n1-ghosts1
  do counter=n1-ghosts1+1,n1
     alpha(counter) = alpha(gi)
     W(counter) = W(gi)
     X(counter) = X(gi)
     rho(counter) = rho(gi)
     vel(counter) = vel(gi)
  enddo

end subroutine analytic_OSC_alpha

subroutine get_OSC_ana_time(eta_guess,M,chi_s,time_guess)

  use GR1D_module, only : pi
  implicit none

  real*8 eta_guess
  real*8 M
  real*8 chi_s
  real*8 time_guess
  
  real*8 eta_s

  eta_s = -2.0d0*acos(cos(eta_guess/2.0d0)/sqrt(cos(chi_s)))
  
  time_guess = 2.0d0*M*((eta_s + pi + (eta_s+pi- &
       sin(eta_s))/(2.0d0*(sin(chi_s))**2))/tan(chi_s)+ &
       log((tan(eta_s/2.0d0)-tan(chi_s))/(tan(eta_s/2.0d0)+ &
       tan(chi_s))))

end subroutine get_OSC_ana_time

subroutine get_OSC_ana_r(chi_guess,radius_o,eta_c,chi_s,radius_guess)

  use GR1D_module, only : pi
  implicit none

  real*8 chi_guess
  real*8 radius_o
  real*8 eta_c,chi_s
  real*8 radius_guess
  
  radius_guess = radius_o*(1.0d0-(cos(eta_c/2.0d0))**2/cos(chi_guess))* &
       sin(chi_guess)/sin(chi_s)

end subroutine get_OSC_ana_r
