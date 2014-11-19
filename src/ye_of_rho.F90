!-*-f90-*-
module ye_of_rho

  implicit none

  real*8,allocatable,save :: yeinlogrho(:)
  real*8,allocatable,save :: rhoinlogrho(:)
  real*8,allocatable,save :: logrhoinlogrho(:)
  real*8,allocatable,save :: slopeinlogrho(:)
  real*8,save :: minrho
  real*8,save :: maxrho
  integer,save :: maxindex

  real*8,save :: yeofrho_logrho1
  real*8,save :: yeofrho_logrho2
  real*8,save :: yeofrho_ye1
  real*8,save :: yeofrho_ye2
  real*8,save :: yeofrho_yec

contains
  subroutine adjust_ye
    
    use GR1D_module
    
    implicit none
  
    real*8 out_ye, delta_ye, delta_S, out_munu
    integer eosflag, keytemp, keyerr, i
    real*8 eosdummy(15)
    
    do i=1,n1
       nuchem(i) = 0.0d0
       if (rho(i).gt.1.0d6*rho_gf) then
          if((.not.bounce).or.(i.gt.ishock(1))) then
             if (do_yeofrhofit) then
                call fit_ye(rho(i)/rho_gf,out_ye)
             else
                call interpolate_ye(rho(i)/rho_gf,out_ye)
             endif
             delta_ye = min(0.0d0, out_ye-ye(i))
             if (eoskey.eq.3) then
                eosflag = 9 !get mu_nu
                keytemp = 0
                call eos(i,rho(i),temp(i),ye(i),eps(i),out_munu,keytemp,keyerr, &
                     eosflag,eoskey,eos_rf_prec)
                nuchem(i) = out_munu
                
                if (out_munu.lt.10.0d0) then
                   delta_S = 0.0d0
                elseif (rho(i).lt.2.0d12*rho_gf) then
                   delta_S = -delta_ye*(out_munu-10.0d0)/temp(i)
                else
                   delta_S = 0.0d0
                endif

                !this updates temp using constant entropy, then finds new energy 
                !associated with that new temp, the one call does both.                
                keytemp = 2
                entropy(i) =  entropy(i) + delta_S
                ye(i) = ye(i) + delta_ye
                keyerr = 0
                call eos_full(i,rho(i),temp(i),ye(i),eps(i),press(i),pressth(i), & 
                     entropy(i),cs2(i),eosdummy(2),&
                     eosdummy(3),eosdummy(4),eosdummy(5),eosdummy(6), &
                     eosdummy(7),eosdummy(8),eosdummy(9),eosdummy(10), &
                     eosdummy(11),eosdummy(12),eosdummy(13),eosdummy(14), &
                     keytemp,keyerr,eoskey,eos_rf_prec)
                if (keyerr.ne.0) then
                   stop "problem in ye_of_rho: eos"
                endif
             else
                ye(i) = ye(i) + delta_ye
             endif
          endif
       endif
    enddo

  end subroutine adjust_ye

  subroutine read_yeprofile

    !in CGS
    use GR1D_module
    implicit none
    
    character(len=100) filename
    integer plen, i, j
    logical proceed
    real*8,allocatable :: profileye(:)
    real*8,allocatable :: profilerho(:)
    real*8 tempvar, tempvar2, tempvar3
    real*8 slope

    filename = trim(adjustl(yeprofile_name))
    open(unit=666,file=trim(adjustl(filename)))    
    
    read(666,"(I6)") plen
    allocate(profileye(plen))
    allocate(profilerho(plen))

    read(666,"(E21.5,E20.5,E12.5)") profilerho(plen), profileye(plen), tempvar
    maxrho = profilerho(plen)
    minrho = profilerho(plen)
    do i=1,plen-1
       read(666,"(E21.5,E20.5,E12.5)") profilerho(plen-i), profileye(plen-i), tempvar
       if (profilerho(plen-i).gt.maxrho) maxrho=profilerho(plen-i)
       if (profilerho(plen-i).lt.minrho) minrho=profilerho(plen-i)
    enddo

    !limits of lookup table
    if (minrho.le.1.0d0) stop "rho very small, adjust method"
    maxindex = (int(log10(maxrho)*100.0d0)+1)
    !allocate lookup table
    allocate(yeinlogrho(maxindex))
    allocate(rhoinlogrho(maxindex))
    allocate(logrhoinlogrho(maxindex))
    allocate(slopeinlogrho(maxindex))

    !density as a function of index is (10.0d0)**(real(index)/100.0)
    do i=1,maxindex
       rhoinlogrho(i) = (10.0d0)**(real(i)/100.0d0)
       logrhoinlogrho(i) = real(i)/100.0d0
    enddo

    !monotonize the ye
    do i=2,plen
       if (profileye(i).gt.profileye(i-1)) then
          profileye(i) = profileye(i-1)
       endif
    enddo

    !pick first lookup index, interate through profile till we straddle the density we want
    yeinlogrho(1) = profileye(1)
    write(*,*) maxindex,maxrho,rhoinlogrho(1333),rhoinlogrho(1332),profilerho(80)
    do i=2,maxindex
       j=1

       proceed = .true.
       if (rhoinlogrho(i).gt.maxrho) then
          proceed = .false.
          yeinlogrho(i) = yeinlogrho(i-1)
       endif
       do while ((j.lt.(plen+1)).and.proceed)
          if (profilerho(j).gt.rhoinlogrho(i)) then
             proceed = .false.
             !set ye it interpolated value between next profile rho and last one
             if (j.eq.1) then
                !special case, can't interpolate, just assume ye does not change
                yeinlogrho(i) = yeinlogrho(i-1)
             else 
                !here we interpolate
                slope = (profileye(j)-profileye(j-1))/(&
                     (log10(profilerho(j))-log10(profilerho(j-1))))
                yeinlogrho(i) = slope*(logrhoinlogrho(i)-log10(profilerho(j)))+profileye(j)
             endif
          endif
          j=j+1
       enddo
       if (j.eq.(plen+1).and.proceed) then
          write(*,*) "Error",i
          stop
       endif
    enddo
    
    do i=1,maxindex-1
       slopeinlogrho(i) = (yeinlogrho(i+1)-yeinlogrho(i))/(logrhoinlogrho(i+1)-logrhoinlogrho(i))
    enddo

  end subroutine read_yeprofile

  subroutine interpolate_ye(incomingrho, outgoingye)

    ! in CGS
    use GR1D_module
    implicit none
    
    real*8 incomingrho,outgoingye
    integer index
    real*8 logincomingrho

    if (incomingrho.lt.1.03d0) then
       stop "Figure something else out"
    end if

    logincomingrho = log10(incomingrho)
    index = int(logincomingrho*100.0d0)


    if (index.gt.maxindex-1) then
       outgoingye = yeinlogrho(maxindex)
    else
       outgoingye = slopeinlogrho(index)*(logincomingrho-logrhoinlogrho(index))+yeinlogrho(index)
    endif

  end subroutine interpolate_ye

  subroutine fit_ye(incomingrho, outgoingye)

    ! in CGS
    implicit none
    
    real*8 incomingrho,outgoingye
    real*8 logincomingrho
    
    real*8 x,absx

    if (incomingrho.lt.1.03d0) then
       stop "Figure something else out"
    end if

    logincomingrho = log10(incomingrho)

    x = max(-1.0d0,min(1.0d0,(2.0d0*logincomingrho - yeofrho_logrho2 - yeofrho_logrho1) &
         /(yeofrho_logrho2-yeofrho_logrho1)))

    absx = abs(x)
    
    outgoingye = 0.5d0*(yeofrho_ye2+yeofrho_ye1) + x/2.0d0*(yeofrho_ye2-yeofrho_ye1) &
         + yeofrho_yec*(1.0d0-absx+4.0d0*absx*(absx-0.5d0)*(absx-1.0d0))

  end subroutine fit_ye

end module ye_of_rho
    
