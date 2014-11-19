!-*-f90-*-
subroutine grid(xmin,xmax,mindx)
  
  use GR1D_module
  
  implicit none
  
  real*8 :: xmin,xmax,mindx
  real*8 :: dxfactor
  real*8 :: dx
  integer :: i,j,gi
  integer :: nconstant, nlog
  
  ! setup innermost updated zone
  real*8 :: inner_grid
  real*8 :: increment
  real*8 :: cellwidth
  
  if(gridtype .eq. "unigrid") then
     dxfactor = 1.0d0
     dx = xmax/(n1-ghosts1*2)
     x1(ghosts1+1) = dx/2.0d0
     x1i(ghosts1+2) = dx
  else if(gridtype .eq. "log") then
     call series2(n1-ghosts1*2,xmin,xmax,mindx,dxfactor)
     dx = mindx
     x1(ghosts1+1) = mindx/2.0d0
     x1i(ghosts1+2) = mindx
  else if (gridtype.eq."custom") then
     nconstant = nint(grid_custom_rad1/grid_custom_dx1)
     if (nconstant.eq.0) then
        stop "grid_custom parameters wrong, or use log"
     endif
     nlog=n1-ghosts1*2-nconstant
     xmin = real(nconstant)*grid_custom_dx1*length_gf
     call series2(nlog,xmin,xmax,grid_custom_dx1*length_gf,dxfactor)
     x1(ghosts1+1) = grid_custom_dx1/2.0d0*length_gf
     x1i(ghosts1+2) = grid_custom_dx1*length_gf
  else if (gridtype.eq."custom2") then       

!!!!!!!!!!!!!!!custom2 setup!!!!!!!!!!!!!!!!
! In custom2 there are 3 zones: Zone 1 begins large and logly gets
! smaller over a number of cells, this is to minimize the effect of
! the 1/r^2 origin issue (or so we think).  Zone 2 is a constant cell
! size zone, this is to maintain resolution around the PNS surface.
! Zone 3 is a log zone that gradually increases as you get farther
! away from the center
!
! Inputs needed:
!    grid_custom_dx1: smallest cell size ~100m
!    grid_custom_inner: size of innermost cell ~1000m
!    grid_custom_number: number of cells in Zone 1 ~10
!    grid_custom_rad1: radius where the constant zoning stops ~20km
!    grid_rmax: maximum radius in simulation, does not have to be exact
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !Part 1: inner zone where size of cells decrease logly
     if (grid_custom_number.eq.0) then
        stop "grid_custom_number cannot be 0, use custom instead"
     endif
     dxfactor = (grid_custom_inner/grid_custom_dx1)** &
          (1.0d0/(real(grid_custom_number)))
     if (dxfactor.gt.1.4d0) then
        stop "Inner grid cells change quickly"
     endif
     x1i(ghosts1+1) = 0.0d0
     gi = grid_custom_number
     do i=ghosts1+2,ghosts1+grid_custom_number+1
        cellwidth = (grid_custom_dx1*length_gf)*dxfactor**(gi)
        x1i(i) = x1i(i-1) + cellwidth
        gi = gi-1
     enddo

     !Part 2: middle zone where width is constant
     nconstant = int((grid_custom_rad1- &
          x1i(ghosts1+grid_custom_number+1)/length_gf)/grid_custom_dx1)
     
     do i=ghosts1+grid_custom_number+2,ghosts1+grid_custom_number+1+nconstant
        cellwidth = grid_custom_dx1*length_gf
        x1i(i) = x1i(i-1)+cellwidth
     enddo

     !Part 3: log increasing grid to the edge of the sim
     !find dx for remaining cells
     nlog = n1-ghosts1*2-nconstant-grid_custom_number
     xmax = grid_rmax*length_gf
     xmin = x1i(ghosts1+grid_custom_number+1+nconstant)
     
     call series2(nlog,xmin,xmax,grid_custom_dx1*length_gf,dxfactor)

     dx = grid_custom_dx1*length_gf
     do i=ghosts1+grid_custom_number+nconstant+2,n1
        x1i(i) = x1i(i-1)+dx
        dx = dx*dxfactor
     enddo

     !Now find position of center of cell
     do i=ghosts1+1,n1-1
        x1(i) = 0.5d0*(x1i(i+1)+x1i(i))
     enddo
     x1(n1) = 2.0d0*x1i(n1) - x1(n1-1)
     
  else
     stop "gridtype not implemented"
  endif

  x1i(ghosts1+1) = 0.0d0

  ! setup grid and outer ghost zones
  if (gridtype.eq."custom2") then
     !everything is done above....
  else if (gridtype.eq."custom") then
     nconstant = nint(grid_custom_rad1/grid_custom_dx1)
     do i=ghosts1+3,ghosts1+nconstant
        dx = grid_custom_dx1*length_gf
        x1i(i) = x1i(i-1) + dx
     enddo
     
     do i=ghosts1+nconstant+1,n1
        dx = dx*dxfactor
        x1i(i) = x1i(i-1)+dx
     enddo
     
     do i=ghosts1+2,n1-1
        x1(i) = 0.5d0*(x1i(i+1)+x1i(i))
     enddo
     
  else
     do i=ghosts1+3,n1
        dx = dx*dxfactor
        x1i(i) = x1i(i-1) + dx
     enddo
     
     do i=ghosts1+2,n1-1
        x1(i) = 0.5d0 * (x1i(i+1)+x1i(i))
     enddo
  endif
  
  x1(n1) = 2.0d0*x1i(n1) - x1(n1-1)
  
  ! setup inner ghost zones
  j=ghosts1+1
  do i=ghosts1,1,-1
     x1(i) = -x1(j)
     x1i(i) = -x1i(j+1)
     j=j+1
  enddo
  
  ! calculate cell coordinate volume
  do i = ghosts1+1,n1-1
     volume(i) = 4.0d0/3.0d0 * pi * ( x1i(i+1)**3 - x1i(i)**3 )
  enddo
  volume(n1) = 4.0d0/3.0d0 * pi * ( &
       (x1i(n1)+(x1i(n1)-x1i(n1-1)))**3 - x1i(n1)**3)
  j=ghosts1+1
  do i=ghosts1,1,-1
     volume(i) = volume(j)
     j=j+1
  enddo
  
end subroutine grid


subroutine series2(nzones,xmin,xmax,mindx,dxfac)
! This routine is motivated by the "series2" subroutine of
! Cala resp. Prometheus.
!
! It solves for a factor dxfac by which each dx is slightly
! larger than the preceding dx.

  implicit none
  
  real*8 dxfac
  integer nzones
  real*8 xmin,xmax,mindx
  
  ! internal vars
  real*8 tol
  real*8 al,aold,ferror,sum,F,dsum,dFda,anew
  integer k,i,itermax
  
  tol = 1.0d-6
  itermax = 100
  
  ! solve for dxfac
  dxfac=0.0d0

  ! estimate
  al = log( (xmax-xmin)/mindx )/ ( real(nzones-2) )
  aold = exp(al)
  k = 1
  ferror = 1.0d0

  !-------------------------------------------------
  ! Solve: F = (xmax-xmin)/mindx - (Sum[ a^j],j=0,N)
  ! let x = a, y(x) = F
  !-------------------------------------------------

  ! evaluate F
  do while( (ferror.gt.tol).and.(k.lt.itermax))
     sum = 0.0d0
     do i=1,nzones-1
        sum = sum + aold**(i-1)
     enddo
     F = ( (xmax-xmin)/mindx) - sum
     
     ! evaluate dFDa
     dsum = 1.0d0
     do i=4,nzones
        dsum = dsum + (i-2)*(aold**(i-3))
     enddo
     dFda = -1.0d0*dsum
     ! next root
     anew = aold - F/dFda
     ferror = abs(anew-aold)/aold
     k = k + 1
     aold = anew
  enddo
  dxfac = anew
  
end subroutine  series2
     
