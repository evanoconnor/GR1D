!-*-f90-*-
subroutine eos_full(i,xrhoi,    &
     xtemp,                   &
     xye,                     &
     xenri,                   &
     xprs,                    &
     xprs_th,                 & 
     xent,                    &
     xcs2,                    &
     xdedt,                   &
     xdpderho,                &
     xdpdrhoe,                & 
     xxa,xxh,xxn,xxp,         & 
     xabar,xzbar,             &
     xmu_e,xmu_n,xmu_p,xmunu,& 
     keytemp,keyerr,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer,intent(in)    :: i
  real*8, intent(inout) :: xrhoi  ! inout because of atmosphere
  real*8, intent(in)    :: xye
  real*8, intent(inout) :: xtemp,xenri,xent
  real*8, intent(out)   :: xprs,xprs_th,xcs2,xdedt
  real*8, intent(out)   :: xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp
  real*8, intent(out)   :: xabar,xzbar,xmu_e,xmu_n,xmu_p,xmunu
  real*8, intent(in)    :: rfeps
  integer, intent(in)   :: keytemp
  integer, intent(in)   :: eoskey
  integer, intent(out)  :: keyerr

  real*8 :: xrho, xenr
  

  ! convert to CGS:
  xrho = xrhoi/rho_gf
  xenr = xenri/eps_gf

!  if(keyerr.eq.0) then
     if(xrhoi.le.atmo_rho_thr) then
        call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
        return
     endif
!  endif

  if(eoskey.eq.1) then
     ! hybrid EOS
     call hybrid_eos(xrho,xenr,xprs,xprs_th,&
          xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs_th*press_gf
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.4) then
     ! Gamma-Law EOS
     call ideal_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.2) then
     ! Poly EOS
     call poly_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)
     xprs = xprs*press_gf
     xprs_th = 0.0d0
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = xmunu
     if(keytemp.eq.1.or.keytemp.eq.2) then
        xenri = xenr*eps_gf
     endif
     if(keyerr.eq.667) then
        if(i.ne.n1) then
           if(atmo(i+1).ne.0) then
              atmo(i) = 1
              call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
              keyerr = 0
              stop "eh"
           endif
        endif
     endif
#else
     stop "eoskey=3 impossible, since NUC_EOS not present."
#endif
  else
     write(6,*) "eoskey ",eoskey," not implemented with eos_full"
     stop "This is bad! Fix me please!"
  endif


end subroutine eos_full

subroutine eos(i,ri,tio,y,eio,xx,keytemp,keyerr,eosflag,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer, intent(in) :: i
  real*8, intent(in) :: y,rfeps
  real*8, intent(inout) :: ri,eio, tio
  real*8, intent(out)    :: xx

  real*8 r,e,p_th,tp
  real*8 px,sx,cs2x,gammax,dedt,ex
  integer eosflag,eoskey,keytemp,keyerr,keyerrt
  !internal
  real*8 prs,soundsqr
  real*8 dpdrho,dpde

#if HAVE_NUC_EOS
  real*8 xmunu,xent,xdedt
#endif

  xx = 0.0d0

  ! convert to CGS:
  r = ri/rho_gf
  e = eio/eps_gf

! eosflags:                                                                     
! 1 --> pressure                                                                
! 2 --> dpdrhoe                                                                 
! 3 --> dpderho                                                                 
! 4 --> eps                                                                     
! 5 --> temp                                                                    
! 6 --> cs2                                                                     
! 7 --> gamma                                                                   
! 8 --> entropy                                                                 
! 9 --> munu             

! eoskey:
! 1 --> Hybrid EOS
! 2 --> Poly EOS
! 3 --> finite-T EOS
! 4 --> Gamma-Law EOS

     if(ri.le.atmo_rho_thr) then
        call atmos_eos(i,ri,prs,eio,soundsqr)
        
        select case (eosflag)
        case(1)
           !pressure
           xx=prs
        case(2)
           stop "eh 2"

        case(3)
           stop "eh 3"

        case(4)
           xx=eio

        case(5)
           stop "ehbah"

        case(6)
           xx=soundsqr
 
        case(7)
           xx=soundsqr*ri/prs

        case(8)
           xx=xx

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

        end select
        
        return
     endif

  if(eoskey.eq.1) then

     select case (eosflag)

        case(1)
           !pressure
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf
           tio=p_th
!           write(*,"(1P10E15.6)") t
        case(2)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for hybrid eos"

        case(6)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=soundsqr*r/prs

        case(8)
           stop "eosflag=8 not implemented for hybrid eos"

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

     end select

  else if(eoskey.eq.2) then

     select case (eosflag)

        case(1)
           !pressure
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for polytropic eos"

        case(6)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for polytropic eos"

        case(8)
           stop "eosflag=8 not implemented for polytropic eos"

        case(9)
           stop "eosflag=9 not implemented for polytropic eos"

        case default
           stop "eosflag not implemented for polytropic eos"

     end select

  else if(eoskey.eq.4) then

     select case (eosflag)

        case(1)
           !pressure
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf


        case(3)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for ideal fluid"

        case(6)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for ideal fluid"

        case(8)
           stop "eosflag=8 not implemented for ideal fluid"

        case(9)
           stop "eosflag=9 not implemented for ideal fluid"

        case default
           stop "eosflag not implemented for ideal fluid"

     end select

  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     select case (eosflag)
        case(1)
           !pressure 
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=prs*press_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(2)
           !dpdrhoe
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=dpdrho*press_gf/rho_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(3)
           !dpde
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=dpde*press_gf/eps_gf
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif

        case(4)
           !eps
           keytemp=1
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for eoskey = 3"

        case(6)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           xx = soundsqr*r/prs

        case(8)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if(keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           xx = xent

        case(9)
           call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
                dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
           if (keytemp.eq.1.or.keytemp.eq.2) then
              eio = e*eps_gf
           endif
           xx = xmunu

        case default
           stop "eosflag not implemented for eoskey = 3"

        end select

        if(keyerr.eq.667) then
           if(i.ne.n1) then
              if(atmo(i+1).ne.0) then
                 atmo(i) = 1
                 keyerr = 0
                 call atmos_eos(i,ri,prs,eio,soundsqr)
                 select case (eosflag)
                 case(1)
                    !pressure
                    xx=prs
                 case(2)
                    stop "eh 4"
                    
                 case(3)
                    stop "eh 5"

                 case(4)
                    xx=eio
                    
                 case(5)
                    stop "ehbah"
                    
                 case(6)
                    xx=soundsqr
                    
                 case(7)
                    xx=soundsqr*ri/prs
                    
                 case(8)
                    xx=xx
                    
                 case(9)
                    stop "eosflag=9 not implemented for hybrid eos"
                    
                 case default
                    stop "eosflag not implemented for hybrid eos"
                    
                 end select
              endif
           endif
        endif

#endif

  else
     write(6,*) "eoskey: ", eoskey
     call flush(6)
     stop "eos choice not implemented, sorry..."

  endif

end subroutine eos

module hybrid_eos_module

  real*8 gamma1
  real*8 gamma2
  real*8 gammath
  real*8 K1,K2,E1,E2,E3
  real*8 rhonuc

end module hybrid_eos_module

subroutine hybrid_eos(rho,enr,prs,pth,dpdrho_o,dpde_o,soundsqr,keytemp)

  use hybrid_eos_module
  implicit none

! Input
  real*8 rho,enr,prs,soundsqr
  integer keytemp
! Output
  real*8 dpdrho_o,dpde_o


! Local 
!
  real*8 pco,pth,eth,dpth_drho,dpth_denr
  real*8 Gx,Ex,Kx,Ex3,dpco_drho,dpco_denr
  real*8 up,dp_drho,dpde

  if(keytemp.eq.1) then
!	energy wanted
     prs = K1*rho**gamma1 
     enr=prs/rho/(gamma1-1.d0)
     soundsqr=gamma1*prs/rho
  endif

  if(rho .lt. rhonuc) then
     Kx=K1
     Ex=E1
     Gx=gamma1
     Ex3=0.d0
  else
     Kx=K2
     Ex=E2
     Gx=gamma2
     Ex3=E3
  endif

!                Thermal
    up=Ex*rho**Gx+Ex3*rho
    eth = enr*rho - up
    pth = (gammath - 1)*eth
    dpth_drho=(gammath - 1)*(enr-Ex*Gx*rho**(Gx-1.d0)-Ex3)
    dpth_denr=(gammath - 1)*rho

    if (pth.lt.0.0d0) then
       pth = 0.0d0
       dpth_drho = 0.0d0
       dpth_denr = 0.0d0
    endif


!                 Cold 
    pco = Kx*rho**Gx
    dpco_drho=Gx*pco/rho
    dpco_denr=0.d0
!
    prs=pco+pth
    dpde=dpco_denr+dpth_denr
    dp_drho=dpco_drho+dpth_drho
    soundsqr= dp_drho+dpde*prs/rho**2

    dpdrho_o = dp_drho
    dpde_o = dpde

end subroutine hybrid_eos


subroutine init_hybrid_eos

  use hybrid_eos_module

  implicit none

  E1 = K1/(gamma1-1.d0)
  E2 = (gamma1 - 1.d0)/(gamma2-1.d0)*E1*rhonuc**(gamma1-gamma2)
  K2 = (gamma2 - 1.d0)*E2
  E3 = (gamma2 - gamma1)/(gamma2-1.d0)*E1*rhonuc**(gamma1-1.d0)

end subroutine init_hybrid_eos

module poly_eos_module

  real*8 polygamma
  real*8 polyK

end module poly_eos_module

subroutine poly_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use poly_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  integer keytemp
  real*8 dpdrho,dpde

! Local 
!
  if(keytemp.eq.1) then
!	energy wanted
     enr=prs/rho/(polygamma-1.d0)
     soundsqr=polygamma*prs/rho
  endif

  prs = polyK * rho**polygamma
  soundsqr = polygamma*polyK*rho**(polygamma-1.d0)

  dpdrho = polyK * rho**(polygamma - 1) * polygamma
  dpde = 0.0d0

end subroutine poly_eos

module ideal_eos_module

  real*8 idealgamma
  real*8 idealK1

end module ideal_eos_module

subroutine ideal_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use ideal_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  real*8 dpdrho,dpde
  integer keytemp


  if(keytemp.eq.1) then
!	energy wanted
     prs=idealK1*rho**(idealgamma)
     enr=prs/rho/(idealgamma-1.d0)
     soundsqr=idealgamma*prs/rho

  endif
     
  prs = (idealgamma - 1.0d0) *rho*enr

  dpde = (idealgamma - 1.0d0 ) * rho
  dpdrho = (idealgamma - 1.0d0 ) * enr

  soundsqr= dpdrho+dpde*prs/rho**2


end subroutine ideal_eos

subroutine atmos_eos(i,xrho,xprs,xenr,xcs2)

  use atmos
  use GR1D_module, only: rho_gf,press_gf, &
       v1,v,W,atmo
  implicit none

  integer, intent(in) :: i
  real*8 :: xrho,xprs, xenr, xcs2
  real*8,parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
  real*8,parameter :: idealgamma = 1.66666666666d0

  xrho = atmo_rho
  xprs=idealK1*((atmo_rho/rho_gf)**(idealgamma))*press_gf
  xenr=xprs/atmo_rho/(idealgamma-1.d0)
  xcs2=idealgamma*xprs/atmo_rho
  v(i) = 0.0d0
  v1(i) = 0.0d0
  W(i) = 1.0d0
  atmo(i) = 1

end subroutine atmos_eos
