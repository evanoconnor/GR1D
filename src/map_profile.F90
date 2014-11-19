! -*-f90-*-
subroutine map_profile(lprofile_name)

  use GR1D_module
  implicit none

  character*(*) lprofile_name
  integer profile_zones
  
  real*8 buffer, dmass, dx
  integer i,ibuffer
  
  integer keytemp,keyerr,eosflag
  
  real*8, allocatable :: pradius(:), &
       pmass(:),prho(:),ptemp(:), &
       ppress(:),peps(:),pvel(:),&
       pye(:),pomega(:)

  real*8,allocatable,save :: pradius_new(:)
  
  real*8 :: kboltz_cgs = 1.380662d-16
  
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

! new short format:  
  if(do_rotation) then
     do i=1,profile_zones
        read(666,*) ibuffer,pmass(i),pradius(i),&
             ptemp(i),prho(i),pvel(i),pye(i), &
             pomega(i)
     enddo
  else
     do i=1,profile_zones
        read(666,*) ibuffer,pmass(i),pradius(i),&
             ptemp(i),prho(i),pvel(i),pye(i), &
             buffer
     enddo
  endif

  ! go to c=G=Msun=1
  if (pmass(profile_zones).lt.1.0d10) then
  else
     pmass = pmass * mass_gf
  endif
  
  !is pradius the cell center or cell interface?  depends on input
  !models, GR1D assumes pradius is cell centers (along with the
  !density, etc.) unless you set WHW02profile to 1, then the code
  !below executes.  GR1D does not use the mass coordinate for creating
  !the initial profile

  if (WHW02) then
     allocate(pradius_new(profile_zones))
     pradius_new(1) = pradius(1)/2.0d0
     do i=2,profile_zones
        pradius_new(i) = (pradius(i)+pradius(i-1))/2.0d0
     enddo
     pradius(:) = pradius_new(:)
     !end of suggested code
  endif

  pradius = pradius * length_gf
  prho = prho * rho_gf
  pvel = pvel / clite
  pomega = pomega / time_gf
  
  if(eoskey.eq.3) then
     ptemp = ptemp * kboltz_cgs / mev_to_erg
  endif
  
  do i=1,n1
     call map_map(rho(i),x1(i),prho,pradius,profile_zones)
     call map_map(eps(i),x1(i),peps,pradius,profile_zones)
     call map_map(press(i),x1(i),ppress,pradius,profile_zones)
     call map_map(temp(i),x1(i),ptemp,pradius,profile_zones)
     call map_map(v1(i),x1(i),pvel,pradius,profile_zones)
     call map_map(ye(i),x1(i),pye,pradius,profile_zones)
     if (do_rotation) then
        call map_map(omega(i),x1(i),pomega,pradius,profile_zones)
     endif
  enddo

  if (do_rotation.and.set_omega) then
     !omega_c and omega_A already in code units
     omega(:) = omega_c/(1.0d0+(x1(:)/omega_A)**2)
  endif


  if(do_rotation) then
     write(*,*) "Have Rotation"
  ! set up vphi
     if(GR) then
        vphi(:) = omega(:)*x1(:)
     else
        vphi1(:) = omega(:)*x1(:)
     endif
  endif

  ! set up atmosphere
  call atmosphere_init
  
  do i=1,n1
         
! call eos
! first: reset eps
     keytemp = 1 ! coming in with temperature
     eosflag = 4 ! we want eps to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),eps(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        write(*,*) keyerr
        write(*,"(i3,i5,1P10E15.6)") eoskey,i,rho(i),temp(i),ye(i),eps(i)
        stop "problem in initial data: eos: eps"
     endif
     
     
! now get new pressures
     keytemp = 0 ! coming in with eps
     eosflag = 1 ! we want press to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),press(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in initial data: eos: press"
     endif

     if(eoskey.eq.3) then
! now get the entropy for analysis purposes
! only do this if we are using the hot nuclear eos!
        keytemp = 1 ! coming in with temperature
        eosflag = 8 ! we want entropy to be reset
        keyerr = 0
        call eos(i,rho(i),temp(i),ye(i),eps(i),entropy(i), &
             keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           stop "problem in initial data: eos: entropy"
        endif
     endif
! speed of sound
     keytemp = 1 ! coming in with temperature
     eosflag = 6 ! we want cs2 to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in initial data: eos: entropy"
     endif
  enddo

end subroutine map_profile

! **************************************************************
subroutine map_linterp(x1,x2,y1,y2,x,y)

! perform linear interpolation      
  implicit none

  real*8 slope,x1,x2,y1,y2,x,y

  if (x2.lt.x1) then
     stop "Error in linterp!"
  endif

  slope = (y2 - y1) / (x2 - x1)

  y = slope*(x-x1) + y1
 
end subroutine  map_linterp

! ***************************************************************

subroutine map_find_index(zones,array,goal,upper_index,lower_index)
  
! bisection search
  implicit none
  
  integer zones,i
  real*8 array(*)
  real*8 goal
  integer middle_index,upper_index,lower_index

  lower_index = 1
  upper_index = zones
  
  do while ( (upper_index - lower_index) .gt. 1 )
     middle_index = (lower_index + upper_index) * 0.5
     if ( (goal .ge. array(lower_index)) &
          .and. (goal .le. array(middle_index)) ) then
        upper_index = middle_index
     else
        if ( (goal .ge. array(middle_index)) &
             .and. (goal .le. array(upper_index)) ) then
           lower_index = middle_index
        endif
     endif
  enddo
      
end subroutine map_find_index

! ******************************************************************

subroutine map_map(point_value,point_radius0,parray,pradius,zones)

  implicit none
  
  real*8 point_value, point_radius, point_radius0
  real*8 pradius(*), parray(*)
  integer zones
  integer upper_index, lower_index

  point_radius = abs(point_radius0)
  
  if (point_radius .ge. pradius(1) .and. & 
       point_radius .lt. pradius(zones) )  then
     
     call map_find_index(zones,pradius,point_radius, &
          upper_index,lower_index)
     
     call map_linterp( pradius(lower_index),pradius(upper_index), &
          parray(lower_index), parray(upper_index),  & 
          point_radius, point_value )

  else if (point_radius .lt. pradius(1)) then
     ! linear extrapolation
     call map_linterp(pradius(1),pradius(2), & 
          parray(1),parray(2),point_radius,point_value)

  else if (point_radius .gt. pradius(zones)) then
     ! linear extrapolation
     call map_linterp(pradius(zones-1),pradius(zones), & 
          parray(zones-1),parray(zones),point_radius,point_value)
  endif
  
  
end subroutine map_map

subroutine map_limit(lprofile_name)
  use GR1D_module, only: grid_rmax,rho_cut
  implicit none
      
  character*(*) lprofile_name
  integer profile_zones
  
  real*8 buffer, dmass, dx
  integer i,ibuffer
  
  integer keytemp,keyerr,eosflag
  
  real*8, allocatable :: pradius(:), &
       pmass(:),prho(:),ptemp(:), &
       ppress(:),peps(:),pvel(:),&
       pye(:),pomega(:)
  
  
  real*8 :: kboltz_cgs = 1.380662d-16

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
  
  do i=1,profile_zones
     if (prho(i).gt.rho_cut) then
        grid_rmax = pradius(i)
     endif
  enddo
  
  write(*,*) "Maximum Radius (at rho=",rho_cut,") set to: ", grid_rmax
  
  close(666)
  
  deallocate(pradius,pmass)
  deallocate(prho,ppress)
  deallocate(ptemp)
  deallocate(peps)
  deallocate(pvel)
  deallocate(pye)
  deallocate(pomega)
  

end subroutine map_limit
