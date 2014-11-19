!-*-f90-*-
      SUBROUTINE intp3d_many_mod( x, y, z, f, ft, nx, ny, nz, nvars, xt, yt, zt)
!
      implicit none
!                                                          
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  8-point Lagrange linear interpolation formula          
!
!     x        input vector of first  variable
!     y        input vector of second variable
!     z        input vector of third  variable
!
!     f        output vector of interpolated function values
!
!     ft       3d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     nz       z-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!     zt       vector of z-coordinates of table
!
!---------------------------------------------------------------------

      integer nx,ny,nz,iv,nvars
      real*8 :: ft(nx,ny,nz,nvars)

      real*8 x,y,z,f(nvars)
      real*8 xt(nx),yt(ny),zt(nz)
      real*8 d1,d2,d3
!
!
      real*8  fh(8,nvars), delx, dely, delz, &
           a1(nvars), a2(nvars), a3(nvars), a4(nvars), &
           a5(nvars), a6(nvars), a7(nvars), a8(nvars)

      real*8 dx,dy,dz,dxi,dyi,dzi,dxyi,dxzi,dyzi,dxyzi
      integer n,ix,iy,iz

!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / dble(nx-1)
      dy    = (yt(ny) - yt(1)) / dble(ny-1)
      dz    = (zt(nz) - zt(1)) / dble(nz-1)
!
      dxi   = 1.0d0 / dx
      dyi   = 1.0d0 / dy
      dzi   = 1.0d0 / dz
!
      dxyi  = dxi * dyi
      dxzi  = dxi * dzi
      dyzi  = dyi * dzi
!
      dxyzi = dxi * dyi * dzi

!
!------- determine location in (equidistant!!!) table 
!                                                                  
      ix = 2 + int( (x - xt(1) - 1.0d-10) * dxi )
      iy = 2 + int( (y - yt(1) - 1.0d-10) * dyi )
      iz = 2 + int( (z - zt(1) - 1.0d-10) * dzi )
!                                                     
      ix = MAX( 2, MIN( ix, nx ) )
      iy = MAX( 2, MIN( iy, ny ) )
      iz = MAX( 2, MIN( iz, nz ) )

!
!------- set-up auxiliary arrays for Lagrange interpolation
!                                                                 
      delx = xt(ix) - x
      dely = yt(iy) - y
      delz = zt(iz) - z
!      
      do iv = 1, nvars
         fh(1,iv) = ft(ix  , iy  , iz, iv  )                             
         fh(2,iv) = ft(ix-1, iy  , iz, iv  )                             
         fh(3,iv) = ft(ix  , iy-1, iz, iv  )                             
         fh(4,iv) = ft(ix  , iy  , iz-1, iv)                             
         fh(5,iv) = ft(ix-1, iy-1, iz, iv  )                             
         fh(6,iv) = ft(ix-1, iy  , iz-1, iv)                             
         fh(7,iv) = ft(ix  , iy-1, iz-1, iv)                             
         fh(8,iv) = ft(ix-1, iy-1, iz-1, iv)                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
            !                                                    
         a1(iv) = fh(1,iv)                             
         a2(iv) = dxi   * ( fh(2,iv) - fh(1,iv) )       
         a3(iv) = dyi   * ( fh(3,iv) - fh(1,iv) )       
         a4(iv) = dzi   * ( fh(4,iv) - fh(1,iv) )       
         a5(iv) = dxyi  * ( fh(5,iv) - fh(2,iv) - fh(3,iv) + fh(1,iv) )
         a6(iv) = dxzi  * ( fh(6,iv) - fh(2,iv) - fh(4,iv) + fh(1,iv) )
         a7(iv) = dyzi  * ( fh(7,iv) - fh(3,iv) - fh(4,iv) + fh(1,iv) )
         a8(iv) = dxyzi * ( fh(8,iv) - fh(1,iv) + fh(2,iv) + fh(3,iv) + &
              fh(4,iv) - fh(5,iv) - fh(6,iv) - fh(7,iv) )
!
         f(iv)  = a1(iv) +  a2(iv) * delx                         &
              +  a3(iv) * dely                         &
              +  a4(iv) * delz                         &
              +  a5(iv) * delx * dely               &
              +  a6(iv) * delx * delz               &
              +  a7(iv) * dely * delz               &
              +  a8(iv) * delx * dely * delz     
!
      enddo
      
    end SUBROUTINE intp3d_many_mod

    SUBROUTINE intp2d_many_mod( x, y, f, ft, nx, ny, nvars, xt, yt)
!
      implicit none
!                                                          
!---------------------------------------------------------------------
!
!     purpose: interpolation of a function of three variables in an
!              equidistant(!!!) table.
!
!     method:  4-point Lagrange linear interpolation formula          
!
!     x        input vector of first  variable
!     y        input vector of second variable
!
!     f        output vector of interpolated function values
!
!     ft       2d array of tabulated function values
!     nx       x-dimension of table
!     ny       y-dimension of table
!     xt       vector of x-coordinates of table
!     yt       vector of y-coordinates of table
!
!---------------------------------------------------------------------

      integer nx,ny,iv,nvars
      real*8 :: ft(nx,ny,nvars)

      real*8 x,y,f(nvars)
      real*8 xt(nx),yt(ny)
      real*8 d1,d2
!
!
      real*8  fh(4,nvars), delx, dely, &
           a1(nvars), a2(nvars), a3(nvars), a4(nvars)
      real*8 dx,dy,dxi,dyi,dxyi
      integer n,ix,iy

!
!------  determine spacing parameters of (equidistant!!!) table
!
      dx    = (xt(nx) - xt(1)) / dble(nx-1)
      dy    = (yt(ny) - yt(1)) / dble(ny-1)
!
      dxi   = 1.0d0 / dx
      dyi   = 1.0d0 / dy
!
      dxyi  = dxi * dyi

!
!------- determine location in (equidistant!!!) table 
!                                                                  
      ix = 2 + int( (x - xt(1) - 1.0d-10) * dxi )
      iy = 2 + int( (y - yt(1) - 1.0d-10) * dyi )
!                                                     
      ix = MAX( 2, MIN( ix, nx ) )
      iy = MAX( 2, MIN( iy, ny ) )

!
!------- set-up auxiliary arrays for Lagrange interpolation
!                                                                 
      delx = xt(ix) - x
      dely = yt(iy) - y
!      
      do iv = 1, nvars
         fh(1,iv) = ft(ix  , iy  , iv  )                             
         fh(2,iv) = ft(ix-1, iy  , iv  )                             
         fh(3,iv) = ft(ix  , iy-1, iv  )                             
         fh(4,iv) = ft(ix-1, iy-1, iv  )                             
!              
!------ set up coefficients of the interpolation polynomial and 
!       evaluate function values 
         a1(iv) = fh(1,iv)                             
         a2(iv) = dxi   * ( fh(2,iv) - fh(1,iv) )       
         a3(iv) = dyi   * ( fh(3,iv) - fh(1,iv) )       
         a4(iv) = dxyi  * ( fh(4,iv) - fh(2,iv) - fh(3,iv) + fh(1,iv) )

         f(iv)  = a1(iv) +  a2(iv) * delx &
              +  a3(iv) * dely &
              +  a4(iv) * delx * dely

      enddo
      
    end SUBROUTINE intp2d_many_mod

