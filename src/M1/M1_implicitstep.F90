!-*-f90-*- 

!this routine returns the new q_M1 after doing the implicit solve. It
!needs q_M1_old; B_M1, C_M1, D_M1 which are constants in terms of the
!implicit step.  They are the old moments, the explicit spatial flux,
!explicit energy flux, and the explicit inelastic scattering terms,
!respectively. The final radiation moments are stored in q_M1. If
!implicit electron scattering, pair-processes of energy coupling is
!turned on, it is done here.

subroutine M1_implicitstep(dts,implicit_factor)

  use GR1D_module
  use nulibtable, only : nulibtable_inv_energies,nulibtable_ewidths, &
       nulibtable_energies,nulibtable_etop,nulibtable_ebottom, &
       nulibtable_logenergies,nulibtable_logetop
  implicit none
  
  real*8 :: dts,implicit_factor
  
  real*8 :: Sr,Stnalpha,Stnum,Stzone,Srzone,oneM1en,oneM1flux,oneeddy
  integer :: sign_one
  real*8 :: RF(2*number_groups),NL_jacobian(2*number_groups,2*number_groups)
  real*8 :: new_RF(2*number_groups),new_NL_jacobian(2*number_groups,2*number_groups)
  real*8 :: old_jacobian(2*number_groups,2*number_groups),old_RF(2*number_groups)
  real*8 :: inverse(2*number_groups,2*number_groups),det
  real*8 :: NLsolve_x(2*number_groups),oldx(2*number_groups)
  real*8 :: pivot(2*number_groups)
  real*8 :: nucubed,nucubedprime,R0out,R0in,R1out,R1in,ies_temp,species_factor,ispecies_factor
  real*8 :: R0pro,R0ann,R1pro,R1ann,epannihil_temp
  real*8 :: ies_sourceterms(2*number_groups)
  real*8 :: epannihil_sourceterms(2*number_groups)
  real*8 :: local_M(2,2),local_J(number_groups),local_H(number_groups,2)
  real*8 :: local_Mbar(2,2),local_Jbar(number_groups),local_Hbar(number_groups,2)
  real*8 :: Ebar(number_groups),Fbar(number_groups),eddybar(number_groups)
  real*8 :: local_L(number_groups,2,2),local_Hdown(number_groups,2)
  real*8 :: local_Ltilde(number_groups,2,2)
  real*8 :: local_dMdE(2,2),local_dJdE(number_groups),local_dHdE(number_groups,2)
  real*8 :: local_dLdE(number_groups,2,2),local_dHdowndE(number_groups,2)
  real*8 :: local_dLtildedE(number_groups,2,2)
  real*8 :: local_dMdF(2,2),local_dJdF(number_groups),local_dHdF(number_groups,2)
  real*8 :: local_dLdF(number_groups,2,2),local_dHdowndF(number_groups,2)
  real*8 :: local_dLtildedF(number_groups,2,2)
  real*8 :: JoverE(number_groups)
  real*8 :: JoverF(number_groups),HoverE(number_groups,2),HoverF(number_groups,2)
  real*8 :: LoverE(number_groups,2,2),LoverF(number_groups,2,2)

  real*8 :: Z(6),Yupr(6),Xuprr(6),Xupff(6),heatterm_NL(6),heattermff_NL(6)
  real*8 :: velocity_coeffs(6,2)
  real*8 :: littlefactors(number_groups,6)
  real*8 :: velocity_center(number_groups,2)
  real*8 :: dmdr,dmdt,dWvuprdr,dWvuprdt,dXdr,dXdt,Kdownrr,dWdt,dWdr
  real*8 :: dvdt(n1),div_v(n1)
  real*8 :: heatterm(number_groups),heattermff(number_groups)
  real*8 :: logdistro(number_groups)
  real*8 :: loginterface_distroj(number_groups)
  real*8 :: interface_distroj(number_groups)
  real*8 :: NLenergyfluxterms(number_groups*2)
  real*8 :: dNLenergyfluxtermsdx(2*number_groups,3)
  real*8 :: dFRdx(2*number_groups),dFLdx(2*number_groups)
  real*8 :: FL(2*number_groups),FR(2*number_groups)
  real*8 :: calculate_enext,xi(number_groups)
  
  real*8 :: M1en_Exp_term(number_groups)
  real*8 :: M1flux_Exp_term(number_groups)
  real*8 :: eddy(number_groups),eddytt(number_groups),eddyff(number_groups),chi(number_groups)

  real*8 :: sourceS(2*number_groups,3)
  real*8 :: sourceG(2*number_groups,3)

  real*8 :: h,invalp,invalp2,invX,invX2,X2,alp2,W2,v2,invr,invr2,onev,oneW,onealp,oneX
  real*8 :: local_u(2),local_littleh(2,2),local_uup(2),local_littlehdowndown(2,2)
  real*8 :: local_littlehupup(2,2)

  integer :: i,j,k,ii,jj,count,j_prime
  integer :: location(1)
  integer :: info

  logical :: nothappenyet1,nothappenyet2,stillneedconvergence
  logical :: problem_fixing,trouble_brewing
  integer :: problem_zone
  integer :: myloc(1)
  real*8 :: maxRF

  problem_fixing = .false.
  problem_zone = 0
  trouble_brewing = .true.

  press_nu = 0.0d0
  energy_nu = 0.0d0
  ynu = 0.0d0

  nothappenyet1 = .true.
  nothappenyet2 = .true.
  !$OMP PARALLEL DO PRIVATE(i,j,div_v,dvdt,M1en_Exp_term,M1flux_Exp_term,eddy,eddytt, &
  !$OMP eddyff,chi,heatterm,heattermff,NLsolve_x,sourceS,sourceG,h,invalp,invalp2, &
  !$OMP alp2,onealp,invX,invX2,X2,oneX,W2,oneW,v2,onev,invr,invr2,local_u,local_uup, &
  !$OMP local_littleh,local_littlehdowndown,local_littlehupup,stillneedconvergence, &
  !$OMP count,ies_sourceterms,epannihil_sourceterms,RF, &
  !$OMP NL_jacobian,local_M,local_J,local_H,local_Mbar,local_Jbar,local_Hbar,local_L, &
  !$OMP local_dMdE,local_dJdE,local_dHdE,local_dLdE,local_dMdF,local_dJdF,local_dHdF, &
  !$OMP local_dLdF,Ebar,Fbar,eddybar,ii,jj,local_Hdown,local_dHdowndE,local_dHdowndF, &
  !$OMP local_Ltilde,local_dLtildedE,local_dLtildedF,j_prime,nucubed, &
  !$OMP nucubedprime,R0pro,R0ann,R1pro,R1ann,epannihil_temp,R0out,R1out,R0in,R1in, &
  !$OMP ies_temp,NLenergyfluxterms,dNLenergyfluxtermsdx,dmdr,dmdt,dXdr,dXdt,Kdownrr, &
  !$OMP dWdt,dWdr,dWvuprdt,dWvuprdr,Z,Yupr,Xuprr,Xupff,heatterm_NL,heattermff_NL, &
  !$OMP velocity_coeffs,littlefactors,velocity_center,logdistro,loginterface_distroj, &
  !$OMP interface_distroj,xi,FL,dFLdx,calculate_enext,FR,dFRdx,inverse,det,old_RF, &
  !$OMP new_NL_jacobian,new_RF,old_jacobian,oldx,myloc,problem_fixing,problem_zone, &
  !$OMP maxRF,Sr,Stnalpha,Stnum,oneM1en,oneM1flux,oneeddy,Stzone,Srzone,sign_one,pivot, &
  !$OMP info,trouble_brewing,species_factor,ispecies_factor)
  do k=ghosts1+1,M1_imaxradii
     do i=1,number_species_to_evolve
        
        if (GR) then
           if (k.eq.ghosts1+1) then
              div_v(k) = v(k)/x1(k)
           else
              div_v(k) = (v(k+1)-v(k-1))/(x1(k+1)-x1(k-1))
           endif
           
           if (nt.eq.0) then
              dvdt(k) = 0.0d0
           else
              dvdt(k) = (v(k)-v_prev(k))/dts
           endif
        else
           if (k.eq.ghosts1+1) then
              div_v(k) = v1(k)/x1(k)
           else
              div_v(k) = (v1(k+1)-v1(k-1))/(x1(k+1)-x1(k-1))
           endif
           
           if (nt.eq.0) then
              dvdt(k) = 0.0d0
           else
              dvdt(k) = (v1(k)-v_prev(k))/dts
           endif
        endif

        !get E and F for this species and grid point
        M1en_Exp_term = B_M1(k,i,:,1) + C_M1(k,i,:,1) + D_M1(k,i,:,1) !on the RHS
        M1flux_Exp_term = B_M1(k,i,:,2) + C_M1(k,i,:,2) + D_M1(k,i,:,2) !on the RHS
        eddy(:) = q_M1(k,i,:,3)
        eddytt(:) = q_M1_extra(k,i,:,1)
        eddyff(:) = q_M1_extra(k,i,:,1)
        chi(:) = q_M1_extra(k,i,:,4)
        heatterm(:) = q_M1_extra(k,i,:,2)
        heattermff(:) = q_M1_extra(k,i,:,3)
        
        !initialize vector of variables to solve
        NLsolve_x(1:number_groups) = q_M1(k,i,:,1)
        NLsolve_x(number_groups+1:2*number_groups) = q_M1(k,i,:,2)
                
        !we'll compute the coefficients of the source terms first,
        !these do not need to be calculated each iteration
        sourceS = 0.0d0
        sourceG = 0.0d0           
        h = (1.0d0+eps(k)+press(k)/rho(k))

        if (GR) then
           invalp = 1.0d0/alp(k)
           invalp2 = invalp**2
           alp2 = alp(k)**2
           onealp = alp(k)
           invX = 1.0d0/X(k)
           invX2 = invX**2
           X2 = X(k)**2
           oneX = X(k)
        else
           if (do_effectivepotential) then
              invalp = 1.0d0/alp(k)
              invalp2 = invalp*invalp
              alp2 = alp(k)**2
              onealp = alp(k)
           else
              invalp = 1.0d0
              invalp2 = 1.0d0
              alp2 = 1.0d0
              onealp = 1.0d0
           endif
           invX = 1.0d0
           invX2 = 1.0d0
           X2 = 1.0d0
           oneX = 1.0d0
        endif

        if (v_order.eq.-1) then
           if (GR) then
              W2 = W(k)**2
              oneW = W(k)
              v2 = v(k)**2
              onev = v(k)
           else
              W2 = 1.0d0/(1.0d0-v1(k)**2)
              oneW = sqrt(W2)
              v2 = v1(k)**2
              onev = v1(k)
           endif
        else if (v_order.eq.0) then
           W2 = 1.0d0
           oneW = 1.0d0
           v2 = 0.0d0
           onev = 0.0d0
           div_v(k) = 0.0d0
           dvdt(k) = 0.0d0
        else
           stop "add in v_order"
        endif
        
        invr = 1.0d0/x1(k)
        invr2 = invr*invr
           
        !h^a_b
        local_u(1) = -oneW*onealp
        local_u(2) = oneW*onev*oneX
        local_uup(1) = -local_u(1)*invalp2
        local_uup(2) = local_u(2)*invX2
           
        local_littleh(1,1) = -v2*W2 !1.0d0-local_u(1)*local_u(1)/alp(k)**2
        local_littleh(2,1) = local_u(2)*invX2*local_u(1)
        local_littleh(1,2) = -local_u(1)*invalp2*local_u(2)
        local_littleh(2,2) = W2 !1.0d0+local_u(2)*local_u(2)/X(k)**2

        local_littlehdowndown(1,1) = -alp2 + local_u(1)*local_u(1)
        local_littlehdowndown(1,2) = local_u(1)*local_u(2)
        local_littlehdowndown(2,1) = local_u(2)*local_u(1)
        local_littlehdowndown(2,2) = X2 + local_u(2)*local_u(2)

        local_littlehupup(1,1) = -invalp2 + local_uup(1)*local_uup(1)
        local_littlehupup(1,2) = local_uup(1)*local_uup(2)
        local_littlehupup(2,1) = local_uup(2)*local_uup(1)
        local_littlehupup(2,2) = invX2 + local_uup(2)*local_uup(2)

        do j=1,number_groups
           
           !these are terms from the RHS of the energy equation
           !that go as the LHS variable (i.e. dEg/dt + kEg ~= 0,
           !\bar{kappa}_E in thesis)
           sourceS(j,1) = implicit_factor*dts*onealp*W2*( &
                eas(k,i,j,2)*oneW+eddy(j)*eas(k,i,j,2)*oneW*v2*invX2 - &
                (eas(k,i,j,2)+eas(k,i,j,3))*v2*oneW*(1.0d0+eddy(j)*invX2))
           if (GR) then
              sourceG(j,1) = 1.0d0 - implicit_factor*dts*onealp*W2*( &
                   4.0d0*pi*x1(k)*rho(k)*h*onev*oneX*(1.0d0+eddy(j)*invX2))
           else
              sourceG(j,1) = 1.0d0
           endif

           !these are terms from the RHS of the energy equation
           !that are independent of Eg and Fg, (i.e. dEg/dt ~=
           !\eta, \bar{\eta}_E in thesis) this also has the
           !spatial flux term (via B) in it and the old energy
           sourceS(j,3) = implicit_factor*dts*onealp*eas(k,i,j,1)*oneW
           sourceG(j,3) = M1en_Exp_term(j) + q_M1_old(k,i,j,1)

           !these are terms from the RHS of the energy equation
           !that go as the flux of the same group (i.e. dEg/dt +
           !(k_a+k_s)Fg ~= 0, -\bar{\psi}_E in thesis, careful
           !of signs here, we are setting up the matrix)
           sourceS(j,2) = -implicit_factor*dts*onealp*W2*(- &
                (eas(k,i,j,2)+eas(k,i,j,3))*oneW*onev*(1.0d0+v2)*invX + &
                2.0d0*eas(k,i,j,2)*oneW*onev*invX)
           if (GR) then
              sourceG(j,2) = implicit_factor*dts*onealp*W2* &
                   4.0d0*pi*x1(k)*rho(k)*h*(1.0d0+v2)
           else
              sourceG(j,2) = 0.0d0
           endif

           !these are terms from  the RHS of the flux equation
          !that go as the LHS variable (i.e. dfg/dt + (k_s+k_a)Fg ~= 0,
           !\bar{kappa}_F in thesis)
           sourceS(j+number_groups,2) = implicit_factor*dts*(onealp*eas(k,i,j,3)* &
                W2*oneW*(1.0d0+v2) + onealp*eas(k,i,j,2)*oneW)
           if (GR) then
              sourceG(j+number_groups,2) = 1.0d0 - implicit_factor*dts* &
                   X(k)*4.0d0*pi*x1(k)*onealp*rho(k)*h*W(k)**2*onev
           else
              sourceG(j+number_groups,2) = 1.0d0
           endif

           !these are terms from the RHS of the flux equation
           !that are independent of Eg and Fg, (i.e. dFg/dt ~=
           !v*\eta, \bar{\eta}_F in thesis) this also has the
           !spatial flux term (via B) in it and the old Flux
           sourceS(j+number_groups,3) = implicit_factor*dts*onealp* &
                oneX*eas(k,i,j,1)*oneW*onev
           sourceG(j+number_groups,3) = M1flux_Exp_term(j) + q_M1_old(k,i,j,2)
           
           !these are terms from the RHS of the flux equation
           !that go as the energy of the same group (i.e. dFg/dt
           !+ -k_a*v*Eg ~= 0, -\bar{\psi}_F in thesis, careful
           !of signs here, we are setting up the matrix)
           sourceS(j+number_groups,1) = -implicit_factor*dts*onealp*(-&
                oneX*eas(k,i,j,2)*W2*oneW*onev - &
                eas(k,i,j,2)*eddy(j)*W2*oneW*v2*onev*invX + &
                oneX*(eas(k,i,j,2)+eas(k,i,j,3))*W2*oneW*onev* &
                (1.0d0+eddy(j)*invX2))

           if (GR) then
              sourceG(j+number_groups,1) = implicit_factor*dts*onealp*( &
                   X2*(mgrav(k)*invr2 + 4.0d0*pi*x1(k)* &
                   (press(k)+rho(k)*h*W2*v2)) - &
                   (eddyff(j)+eddytt(j))*invr)
           else
              sourceG(j+number_groups,1) = implicit_factor*dts*onealp*( -&
                   (eddyff(j)+eddytt(j))*invr)
           endif

        enddo

        !start implicit loop
        stillneedconvergence = .true.
        count = 0
        do while (stillneedconvergence)
           !reset calculated variables
           ies_sourceterms = 0.0d0
           epannihil_sourceterms = 0.0d0
           RF = 0.0d0
           NL_jacobian = 0.0d0       

           !ep annihilation or Ielectron scattering determine fluid terms
           if ((include_epannihil_kernels.and.i.eq.3).or.include_Ielectron_imp) then
              local_M = 0.0d0
              local_J = 0.0d0
              local_H = 0.0d0
              local_Mbar = 0.0d0
              local_Jbar = 0.0d0
              local_Hbar = 0.0d0
              local_L = 0.0d0
              local_dMdE = 0.0d0
              local_dJdE = 0.0d0
              local_dHdE = 0.0d0
              local_dLdE = 0.0d0
              local_dMdF = 0.0d0
              local_dJdF = 0.0d0
              local_dHdF = 0.0d0
              local_dLdF = 0.0d0

              if (i.eq.3.and.number_species.eq.3) then
                 species_factor = 4.0d0
                 ispecies_factor = 0.25
              else if (i.eq.3) then
                 stop "add in proper i>3 species support for thermal process"
              else
                 species_factor = 1.0d0
                 ispecies_factor = 1.0d0
              endif

              if (i.eq.1) then
                 !implicit neutrino pair production/annihilation isn't tested for i.eq.1 or i.eq.2
                 Ebar(:) = q_M1_old(k,2,:,1) !antinue at time (n), because order of species loop
                 Fbar(:) = q_M1_old(k,2,:,2) !antinue at time (n), because order of species loop
                 eddybar(:) = q_M1_old(k,2,:,3) !antinue at time (n), because order of species loop
              else if (i.eq.2) then
                 !implicit neutrino pair production/annihilation isn't tested for i.eq.1 or i.eq.2
                 Ebar(:) = q_M1_old(k,1,:,1) !nue at time (n+1), because order of species loop
                 Fbar(:) = q_M1_old(k,1,:,2) !nue at time (n+1), because order of species loop
                 eddybar(:) = q_M1_old(k,1,:,3) !nue at time (n+1), because order of species loop
              else if (i.eq.3.and.number_species.eq.3) then
                 Ebar(:) = NLsolve_x(1:number_groups) !nux at time (n)
                 Fbar(:) =  NLsolve_x(number_groups+1:2*number_groups)!nux at time (n)
                 eddybar(:) = q_M1(k,i,:,3) !nux at time (n)
              else
                 stop "add in i>3 species support for thermal processes"
              endif

              !linearize the block terms
              do j=1,number_groups
                 !M^{ab}
                 local_M(1,1) = NLsolve_x(j)*invalp2*ispecies_factor
                 local_M(1,2) = NLsolve_x(j+number_groups)*invX2*invalp*ispecies_factor
                 local_M(2,1) = local_M(1,2)
                 local_M(2,2) = eddy(j)*NLsolve_x(j)*invX2**2*ispecies_factor

                 local_Mbar(1,1) = Ebar(j)*invalp2*ispecies_factor
                 local_Mbar(1,2) = Fbar(j)*invX2*invalp*ispecies_factor
                 local_Mbar(2,1) = local_Mbar(1,2)
                 local_Mbar(2,2) = eddybar(j)*Ebar(j)*invX2**2*ispecies_factor

                 local_dMdE(1,1) = invalp2
                 local_dMdF(1,2) = invX2*invalp
                 local_dMdF(2,1) = local_dMdF(1,2)
                 local_dMdE(2,2) = eddy(j)*invX2**2
                 

                 do ii=1,2
                    do jj=1,2
                       local_J(j) = local_J(j) + local_M(ii,jj)*local_u(ii)*local_u(jj)
                       local_Jbar(j) = local_Jbar(j) + local_Mbar(ii,jj)*local_u(ii)*local_u(jj)
                       local_dJdE(j) = local_dJdE(j) + local_dMdE(ii,jj)*local_u(ii)*local_u(jj)
                       local_dJdF(j) = local_dJdE(j) + local_dMdF(ii,jj)*local_u(ii)*local_u(jj)

                       !H^{a}
                       local_H(j,1) = local_H(j,1) - local_M(ii,jj)*local_u(ii)*local_littleh(1,jj)
                       local_Hbar(j,1) = local_Hbar(j,1) - local_Mbar(ii,jj)*local_u(ii)*local_littleh(1,jj)
                       local_dHdE(j,1) = local_dHdE(j,1) - local_dMdE(ii,jj)*local_u(ii)*local_littleh(1,jj)
                       local_dHdF(j,1) = local_dHdF(j,1) - local_dMdF(ii,jj)*local_u(ii)*local_littleh(1,jj)
                       local_H(j,2) = local_H(j,2) - local_M(ii,jj)*local_u(ii)*local_littleh(2,jj)
                       local_Hbar(j,2) = local_Hbar(j,2) - local_Mbar(ii,jj)*local_u(ii)*local_littleh(2,jj)
                       local_dHdE(j,2) = local_dHdE(j,2) - local_dMdE(ii,jj)*local_u(ii)*local_littleh(2,jj)
                       local_dHdF(j,2) = local_dHdF(j,2) - local_dMdF(ii,jj)*local_u(ii)*local_littleh(2,jj)
                    enddo
                 enddo

                 do ii=1,2
                    do jj=1,2
                       !L^{ab} -- thin limit and thick limit interpolated with chi!
                       local_L(j,1,1) = local_L(j,1,1) + &
                            (local_M(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,1)*local_J(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdE(j,1,1) = local_dLdE(j,1,1) + &
                            (local_dMdE(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,1)*local_dJdE(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdF(j,1,1) = local_dLdF(j,1,1) + &
                            (local_dMdF(ii,jj)*local_littleh(1,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,1)*local_dJdF(j)*onethird)*(1.5d0-1.5d0*chi(j))

                       local_L(j,1,2) = local_L(j,1,2) + &
                            (local_M(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,2)*local_J(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdE(j,1,2) = local_dLdE(j,1,2) + &
                            (local_dMdE(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,2)*local_dJdE(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdF(j,1,2) = local_dLdF(j,1,2) + &
                            (local_dMdF(ii,jj)*local_littleh(1,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(1,2)*local_dJdF(j)*onethird)*(1.5d0-1.5d0*chi(j))

                       local_L(j,2,1) = local_L(j,2,1) + &
                            (local_M(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,1)*local_J(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdE(j,2,1) = local_dLdE(j,2,1) + &
                            (local_dMdE(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,1)*local_dJdE(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdF(j,2,1) = local_dLdF(j,2,1) + &
                            (local_dMdF(ii,jj)*local_littleh(2,ii)*local_littleh(1,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,1)*local_dJdF(j)*onethird)*(1.5d0-1.5d0*chi(j))

                       local_L(j,2,2) = local_L(j,2,2) + &
                            (local_M(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,2)*local_J(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdE(j,2,2) = local_dLdE(j,2,2) + &
                            (local_dMdE(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,2)*local_dJdE(j)*onethird)*(1.5d0-1.5d0*chi(j))
                       local_dLdF(j,2,2) = local_dLdF(j,2,2) + &
                            (local_dMdF(ii,jj)*local_littleh(2,ii)*local_littleh(2,jj))*(1.5d0*chi(j)-0.5d0) + &
                            (local_littlehupup(2,2)*local_dJdF(j)*onethird)*(1.5d0-1.5d0*chi(j))
                    enddo
                 enddo

                 local_Hdown(j,1) = sum(local_H(j,:)*local_littlehdowndown(1,:))
                 local_dHdowndE(j,1) = sum(local_dHdE(j,:)*local_littlehdowndown(1,:))
                 local_dHdowndF(j,1) = sum(local_dHdF(j,:)*local_littlehdowndown(1,:))

                 local_Hdown(j,2) = sum(local_H(j,:)*local_littlehdowndown(2,:))
                 local_dHdowndE(j,2) = sum(local_dHdE(j,:)*local_littlehdowndown(2,:))
                 local_dHdowndF(j,2) = sum(local_dHdF(j,:)*local_littlehdowndown(2,:))

                 !L^a_b tilde as in Shibata et al. 2011
                 local_Ltilde(j,1,1) = -local_L(j,1,1)*alp2 - local_J(j)*local_littleh(1,1)*onethird
                 local_dLtildedE(j,1,1) = -local_dLdE(j,1,1)*alp2 - local_dJdE(j)*local_littleh(1,1)*onethird
                 local_dLtildedF(j,1,1) = -local_dLdF(j,1,1)*alp2 - local_dJdF(j)*local_littleh(1,1)*onethird

                 local_Ltilde(j,1,2) = local_L(j,1,2)*X2 - local_J(j)*local_littleh(1,2)*onethird
                 local_dLtildedE(j,1,2) = local_dLdE(j,1,2)*X2 - local_dJdE(j)*local_littleh(1,2)*onethird
                 local_dLtildedF(j,1,2) = local_dLdF(j,1,2)*X2 - local_dJdF(j)*local_littleh(1,2)*onethird

                 local_Ltilde(j,2,1) = -local_L(j,2,1)*alp2 - local_J(j)*local_littleh(2,1)*onethird
                 local_dLtildedE(j,2,1) = -local_dLdE(j,2,1)*alp2 - local_dJdE(j)*local_littleh(2,1)*onethird
                 local_dLtildedF(j,2,1) = -local_dLdF(j,2,1)*alp2 - local_dJdF(j)*local_littleh(2,1)*onethird

                 local_Ltilde(j,2,2) = local_L(j,2,2)*X2 - local_J(j)*local_littleh(2,2)*onethird
                 local_dLtildedE(j,2,2) = local_dLdE(j,2,2)*X2 - local_dJdE(j)*local_littleh(2,2)*onethird
                 local_dLtildedF(j,2,2) = local_dLdF(j,2,2)*X2 - local_dJdF(j)*local_littleh(2,2)*onethird

              enddo
              
           endif

           !ep-annihilation for nux only
           if (include_epannihil_kernels.and.i.eq.3) then
              if (i.ne.3) stop "check update_eas to see if all kernels are getting interpolated i.ne.3 is highly experimental"

              do j=1,number_groups !\omega
                 do j_prime=1,number_groups !\omega^\prime, integrate over this

                    !\nu^3 in units of E,F.  If f( or fprime)=1,
                    !nucubed (or nucubedprime) = J (or Jprime). I do
                    !not include the 4\pi with the \nu^3 term seen in
                    !Shibata et al. (2011) as my variables are /srad.
                    !nucubed has the \delta E for the bin width in it
                    !already, so do the E and J, therefore we don't
                    !have to worry about any bin widths, just add up
                    !the individual contributions

                    !J should never be bigger than nu3, if it is we
                    !are over populated, this happens in high
                    !velocity, degenerate conditions.  It never gets
                    !too big, like <1% above 1

                    !we multiply by 4pi because we just integrated
                    !over the other neutrino. We only integrated over
                    !energy, this is 4\pi is for integraing over solid
                    !angle (which is done kinda, but still off by
                    !4\pi)

                    !we include the species_factor so that the
                    !emission is 4x if we are using nux =
                    !numu+anumu+nutau+anutau

                    nucubed = M1_moment_to_distro_inverse(j) 
                    nucubedprime = M1_moment_to_distro_inverse(j_prime)

                    R0pro = 0.5d0*epannihil(k,i,j,j_prime,1)
                    R0ann = 0.5d0*epannihil(k,i,j,j_prime,2)

                    R1pro = 1.5d0*epannihil(k,i,j,j_prime,3)
                    R1ann = 1.5d0*epannihil(k,i,j,j_prime,4)

                    if (R0pro.lt.0.0d0) stop "R0pro should not be less than 0"
                    if (R0ann.lt.0.0d0) stop "R0ann should not be less than 0"

                    !shibata eq. 4.19, alpha=t , evaluate, time by 4*pi*dt*alp^2 and move to LHS for RF
                    epannihil_temp = -species_factor*implicit_factor*dts*alp2*4.0d0*pi*( - &
                         ((local_J(j)-nucubed)*local_uup(1) + local_H(j,1))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,1)*onethird*((nucubed-local_J(j))*R1pro + local_J(j)*R1ann) + &
                         (sum(local_Hdown(j,:)*local_Hbar(j_prime,:))*local_uup(1) + &
                         sum(local_Ltilde(j,1,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_J(j)*local_uup(1)+local_H(j,1))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                    RF(j) = RF(j) + epannihil_temp
                    epannihil_sourceterms(j) = epannihil_sourceterms(j) - epannihil_temp

                    !dStdE
                    NL_jacobian(j,j) = NL_jacobian(j,j) - species_factor*implicit_factor*dts*alp2*4.0d0*pi*( - &
                         (local_dJdE(j)*local_uup(1) + local_dHdE(j,1))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,1)*onethird*(-local_dJdE(j)*R1pro + local_dJdE(j)*R1ann) + &
                         (sum(local_dHdowndE(j,:)*local_Hbar(j_prime,:))*local_uup(1) + &
                         sum(local_dLtildedE(j,1,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_dJdE(j)*local_uup(1)+local_dHdE(j,1))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                    !dStdF
                    NL_jacobian(j,j+number_groups) = NL_jacobian(j,j+number_groups) - &
                         species_factor*implicit_factor*dts*alp2*4.0d0*pi*( - &
                         ((local_dJdF(j))*local_uup(1) + local_dHdF(j,1))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,1)*onethird*(-local_dJdF(j)*R1pro + local_dJdF(j)*R1ann) + &
                         (sum(local_dHdowndF(j,:)*local_Hbar(j_prime,:))*local_uup(1) + &
                         sum(local_dLtildedF(j,1,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_dJdF(j)*local_uup(1)+local_dHdF(j,1))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                    !shibata eq. 4.19, alpha=r , evaluate, time by 4*pi*dt*alp*X2 and move to LHS for RF
                    epannihil_temp = -species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( - &
                         ((local_J(j)-nucubed)*local_uup(2) + local_H(j,2))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,2)*onethird*((nucubed-local_J(j))*R1pro + local_J(j)*R1ann) + &
                         (sum(local_Hdown(j,:)*local_Hbar(j_prime,:))*local_uup(2) + &
                         sum(local_Ltilde(j,2,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_J(j)*local_uup(2)+local_H(j,2))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                    RF(j+number_groups) = RF(j+number_groups) + epannihil_temp

                    epannihil_sourceterms(j+number_groups) = epannihil_sourceterms(j+number_groups) - &
                         epannihil_temp

                    !dSrdE
                    NL_jacobian(j+number_groups,j) = NL_jacobian(j+number_groups,j) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( - &
                         (local_dJdE(j)*local_uup(2) + local_dHdE(j,2))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,2)*onethird*(-local_dJdE(j)*R1pro + local_dJdE(j)*R1ann) + &
                         (sum(local_dHdowndE(j,:)*local_Hbar(j_prime,:))*local_uup(2) + &
                         sum(local_dLtildedE(j,2,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_dJdE(j)*local_uup(2)+local_dHdE(j,2))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                    !dSrdF
                    NL_jacobian(j+number_groups,j+number_groups) = &
                         NL_jacobian(j+number_groups,j+number_groups) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( - &
                         ((local_dJdF(j))*local_uup(2) + local_dHdF(j,2))* &
                         (nucubedprime-local_Jbar(j_prime))*R0pro - &
                         local_Hbar(j_prime,2)*onethird*(-local_dJdF(j)*R1pro + local_dJdF(j)*R1ann) + &
                         (sum(local_dHdowndF(j,:)*local_Hbar(j_prime,:))*local_uup(2) + &
                         sum(local_dLtildedF(j,2,:)*local_Hbar(j_prime,:)))*(R1pro-R1ann) - &
                         (local_dJdF(j)*local_uup(2)+local_dHdF(j,2))*local_Jbar(j_prime)*R0ann)* &
                         nulibtable_inv_energies(j_prime)

                 enddo
                 epannihil_sourceterm(k,i,j,1) = epannihil_sourceterms(j)/(implicit_factor*dts*alp2)
                 epannihil_sourceterm(k,i,j,2) = epannihil_sourceterms(j+number_groups)/(implicit_factor*dts*X2)
              enddo
              
           endif

           !implicit scattering
           if (include_Ielectron_imp) then
              do j=1,number_groups !\omega
                 do j_prime=1,number_groups !\omega^\prime, integrate over this

                    !\nu^3 in units of E,F.  If f(or fprime)=1,
                    !nucubed(or nucubedprime) = J(or Jprime). I do not
                    !include the 4\pi on the \nu^3 seen in Shibata et
                    !al. (2011) as my variables are /srad.  nucubed
                    !has the \delta E for the bin width in it already,
                    !so do the E and J, therefore we don't have to
                    !worry about any bin widths, just add up the
                    !individual contributions

                    !J should never be bigger than nu3, if it is we are
                    !over populated, this happens in high velocity,
                    !degenerate conditions.  It never gets too big.

                    !multiply by 4pi because we just integrated over
                    !incoming neutrino. We only integrated over energy, this
                    !is integraing over solid angle (which is done kinda,
                    !but still off by 4\pi)

                    nucubed = M1_moment_to_distro_inverse(j) 
                    nucubedprime = M1_moment_to_distro_inverse(j_prime)

                    R0out = 0.5d0*ies(k,i,j,j_prime,1)
                    R1out = 1.5d0*ies(k,i,j,j_prime,2)
                    if (R0out.lt.0.0d0) stop "R0out should not be less than 0"
                    R0in = 0.5d0*ies(k,i,j_prime,j,1)
                    R1in = 1.5d0*ies(k,i,j_prime,j,2)

                    !shibata eq. 4.14, alpha=t , evaluate, time by 4*pi*dt*alp^2 and move to LHS for RF
                    ies_temp = -species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(1) - local_H(j,1))*R0in*local_J(j_prime) + &
                         local_H(j_prime,1)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_H(j_prime,:))*local_uup(1) + &
                         sum(local_Ltilde(j,1,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_J(j)*local_uup(1)+local_H(j,1))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime)
                    RF(j) = RF(j) + ies_temp

                    ies_sourceterms(j) = ies_sourceterms(j) - ies_temp

                    !dStdE
                    NL_jacobian(j,j) = NL_jacobian(j,j) - species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                         (-local_dJdE(j)*local_uup(1) - local_dHdE(j,1))*R0in*local_J(j_prime) + &
                         local_H(j_prime,1)*onethird*local_dJdE(j)*(R1out-R1in) - &
                         (sum(local_Hdown(j_prime,:)*local_dHdE(j,:))*local_uup(1) + &
                         sum(local_dLtildedE(j,1,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_dJdE(j)*local_uup(1)+local_dHdE(j,1))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime) 

                    !dStdEprime
                    NL_jacobian(j,j_prime) = NL_jacobian(j,j_prime) - &
                         species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(1) - local_H(j,1))*R0in*local_dJdE(j_prime) + &
                         local_dHdE(j_prime,1)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_dHdE(j_prime,:))*local_uup(1) + &
                         sum(local_Ltilde(j,1,:)*local_dHdE(j_prime,:)))*(R1in-R1out) + &
                         (local_J(j)*local_uup(1)+local_H(j,1))*local_dJdE(j_prime)*R0out)* &
                         nulibtable_inv_energies(j_prime)
                    
                    !dStdF
                    NL_jacobian(j,j+number_groups) = NL_jacobian(j,j+number_groups) - &
                         species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                         (-local_dJdF(j)*local_uup(1) - local_dHdF(j,1))*R0in*local_J(j_prime) + &
                         local_H(j_prime,1)*onethird*local_dJdF(j)*(R1out-R1in) - &
                         (sum(local_Hdown(j_prime,:)*local_dHdF(j,:))*local_uup(1) + &
                         sum(local_dLtildedF(j,1,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_dJdF(j)*local_uup(1)+local_dHdF(j,1))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime) 

                    !dStdFprime
                    NL_jacobian(j,j_prime+number_groups) = NL_jacobian(j,j_prime+number_groups) - &
                         species_factor*implicit_factor*dts*alp2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(1) - local_H(j,1))*R0in*local_dJdF(j_prime) + &
                         local_dHdF(j_prime,1)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_dHdF(j_prime,:))*local_uup(1) + &
                         sum(local_Ltilde(j,1,:)*local_dHdF(j_prime,:)))*(R1in-R1out) + &
                         (local_J(j)*local_uup(1)+local_H(j,1))*local_dJdF(j_prime)*R0out)* &
                         nulibtable_inv_energies(j_prime)

                    !shibata eq. 4.14, alpha=r , evaluate, time by 4*pi*dt*alp^2 and move to LHS for RF
                    ies_temp = -species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(2) - local_H(j,2))*R0in*local_J(j_prime) + &
                         local_H(j_prime,2)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_H(j_prime,:))*local_uup(2) + &
                         sum(local_Ltilde(j,2,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_J(j)*local_uup(2)+local_H(j,2))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime)

                    RF(j+number_groups) = RF(j+number_groups) + ies_temp
                    ies_sourceterms(j+number_groups) = ies_sourceterms(j+number_groups) - ies_temp

                    !dSrdE
                    NL_jacobian(j+number_groups,j) = NL_jacobian(j+number_groups,j) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                         (-local_dJdE(j)*local_uup(2) - local_dHdE(j,2))*R0in*local_J(j_prime) + &
                         local_H(j_prime,2)*onethird*local_dJdE(j)*(R1out-R1in) - &
                         (sum(local_Hdown(j_prime,:)*local_dHdE(j,:))*local_uup(2) + &
                         sum(local_dLtildedE(j,2,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_dJdE(j)*local_uup(2)+local_dHdE(j,2))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime) 
                    
                    !dSrdEprime
                    NL_jacobian(j+number_groups,j_prime) = NL_jacobian(j+number_groups,j_prime) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(2) - local_H(j,2))*R0in*local_dJdE(j_prime) + &
                         local_dHdE(j_prime,2)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_dHdE(j_prime,:))*local_uup(2) + &
                         sum(local_Ltilde(j,2,:)*local_dHdE(j_prime,:)))*(R1in-R1out) + &
                         (local_J(j)*local_uup(2)+local_H(j,2))*local_dJdE(j_prime)*R0out)* &
                         nulibtable_inv_energies(j_prime)                         

                    !dSrdF     
                    NL_jacobian(j+number_groups,j+number_groups) = &
                         NL_jacobian(j+number_groups,j+number_groups) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                         (-local_dJdF(j)*local_uup(2) - local_dHdF(j,2))*R0in*local_J(j_prime) + &
                         local_H(j_prime,2)*onethird*local_dJdF(j)*(R1out-R1in) - &
                         (sum(local_Hdown(j_prime,:)*local_dHdF(j,:))*local_uup(2) + &
                         sum(local_dLtildedF(j,2,:)*local_H(j_prime,:)))*(R1in-R1out) - &
                         (local_dJdF(j)*local_uup(2)+local_dHdF(j,2))*(nucubedprime-local_J(j_prime))*R0out)* &
                         nulibtable_inv_energies(j_prime) 

                    
                    !dSrdFprime
                    NL_jacobian(j+number_groups,j_prime+number_groups) = &
                         NL_jacobian(j+number_groups,j_prime+number_groups) - &
                         species_factor*implicit_factor*dts*onealp*X2*4.0d0*pi*( &
                         ((nucubed-local_J(j))*local_uup(2) - local_H(j,2))*R0in*local_dJdF(j_prime) + &
                         local_dHdF(j_prime,2)*onethird*((nucubed-local_J(j))*R1in + local_J(j)*R1out) - &
                         (sum(local_Hdown(j,:)*local_dHdF(j_prime,:))*local_uup(2) + &
                         sum(local_Ltilde(j,2,:)*local_dHdF(j_prime,:)))*(R1in-R1out) + &
                         (local_J(j)*local_uup(2)+local_H(j,2))*local_dJdF(j_prime)*R0out)* &
                         nulibtable_inv_energies(j_prime)                         
                    
                 enddo
                 ies_sourceterm(k,i,j,1) = ies_sourceterms(j)/(implicit_factor*dts*alp2)
                 ies_sourceterm(k,i,j,2) = ies_sourceterms(j+number_groups)/(implicit_factor*dts*X2)
              enddo
           endif

           NLenergyfluxterms = 0.0d0
           dNLenergyfluxtermsdx = 0.0d0           
           
           if (include_energycoupling_imp) then
              !Warning, not heavily tested, I work a lot with the
              !explicit version, this only gets used a bit near BH
              !formation
              
              !Mueller's method:
              if (GR) then
                 dmdr = 4.0d0*pi*x1(k)**2*(rho(k)*h*W(k)**2-press(k))
                 dmdt = -4.0d0*pi*x1(k)**2*alp(k)*rho(k)*h*W(k)**2*v(k)*invX
                 dXdr = X2*oneX*(dmdr*invr-mgrav(k)*invr2)
                 dXdt = X2*oneX*dmdt*invr
                 Kdownrr = -oneX*dXdt*invalp
                 dWdt = W2*oneW*onev*dvdt(k)
                 dWdr = W2*oneW*onev*div_v(k)
                 dWvuprdt = onev*invX*dWdt - oneW*onev*invX2*dXdt + oneW*invX*dvdt(k)
                 dWvuprdr = onev*invX*dWdr - oneW*onev*invX2*dXdr + oneW*invX*div_v(k)
              else
                 dmdr = 0.0d0
                 dmdt = 0.0d0
                 dXdr = 0.0d0
                 dXdt = 0.0d0
                 Kdownrr = 0.0d0
                 dWdt = W2*oneW*onev*dvdt(k)
                 dWdr = W2*oneW*onev*div_v(k)
                 dWvuprdt = onev*dWdt + oneW*dvdt(k)
                 dWvuprdr = onev*dWdr + oneW*div_v(k)
              endif

              Z = 0.0d0
              Yupr = 0.0d0
              Xuprr = 0.0d0
              Xupff = 0.0d0
              heatterm_NL = 0.0d0
              heattermff_NL = 0.0d0
              
              !g_rr v^r v^r = v^2 -> v = Xv^r -> v = v_r/X -> v_r = vx -> v^r = v/x
              if (.not.GR) stop "fix GR in here"

              Z(1) = 1.0d0/oneW !E
              Z(2) = onev/(oneX*oneW) !F
              Z(3) = v2/(X2*oneW) !prr
              Z(5) = v2*onev*invX2*invX !qrrr
              Yupr(2) = 1.0d0/(X2*oneW) !F
              Yupr(3) = onev/(X2*oneX*oneW)
              Yupr(5) = v2*invX2**2
              Xuprr(3) = 1.0d0/(X2*X2*oneW)
              Xuprr(5) = onev*invX2**2*invX
              Xupff(4) = 1.0d0/(x1(k)**2*oneW)
              Xupff(6) = onev*invX*invr2
              heatterm_NL(5) = 1.0d0 
              heattermff_NL(6) = 1.0d0

              !determine velocity coefficients
              !energy
              if (GR) then
                 velocity_coeffs(:,1) = onealp*(oneW*((Z(:)*onev*invX-Yupr(:))*dphidr(k) - &
                      Xuprr(:)*onev*dXdr - 2.0d0*Xupff(:)*onev*invX*x1(k) + &
                      Xuprr(:)*Kdownrr) + Z(:)*invalp*dWdt + Yupr(:)*dWdr - &
                      X2*Yupr(:)*invalp*dWvuprdt - X2*Xuprr(:)*dWvuprdr)
                 
                 !flux
                 velocity_coeffs(:,2) = onealp*(oneW*((Yupr(:)*onev*oneX-Xuprr(:)*X2)*dphidr(k) - &
                      heatterm_NL(:)*invX2**2*onev*dXdr - 2.0d0*heattermff_NL(:)/x1(k)**2*onev*invX*x1(k) + &
                      heatterm_NL(:)*invX2**2*Kdownrr) + Yupr(:)*X2*invalp*dWdt +  &
                      Xuprr(:)*X2*dWdr - Xuprr(:)*X2*X2*invalp*dWvuprdt - &
                      heatterm_NL(:)*invX2*dWvuprdr)
              else
                 stop "here too"
                 velocity_coeffs(:,1) = onealp*(oneW*(-2.0d0*Xupff(:)*onev*invX*x1(k)) + &
                      Z(:)*invalp*dWdt + Yupr(:)*dWdr - &
                      X2*Yupr(:)*invalp*dWvuprdt - X2*Xuprr(:)*dWvuprdr)
                 
                 !flux
                 velocity_coeffs(:,2) = onealp*(oneW*(-2.0d0*heattermff_NL(:)*invr2*onev*invX*x1(k)) + &
                      Yupr(:)*X2*invalp*dWdt +  &
                      Xuprr(:)*X2*dWdr - Xuprr(:)*X2*X2*invalp*dWvuprdt - &
                      heatterm_NL(:)*invX2*dWvuprdr)
              endif

              !interpolate interface value of F/E, p, heatterms/E
              do j=1,number_groups
                 !E/E, F/E,prr,pff,Qrrr/E,Qffr/E
                 littlefactors(j,1) = 1.0d0
                 littlefactors(j,2) = NLsolve_x(j+number_groups)/NLsolve_x(j)
                 littlefactors(j,3) = eddy(j)
                 littlefactors(j,4) = eddyff(j)
                 littlefactors(j,5) = heatterm(j)/NLsolve_x(j)
                 littlefactors(j,6) = heattermff(j)/NLsolve_x(j)
                 
                 velocity_center(j,1) = sum(littlefactors(j,:)*velocity_coeffs(:,1))
                 velocity_center(j,2) = sum(littlefactors(j,:)*velocity_coeffs(:,2))
              enddo

              !get geometric interpolated distributions functions, not
              !in the fluid frame.... wonder if this will be a
              !problem... its just for determiing xi, maybe not.
              logdistro = log(NLsolve_x(1:number_groups)*M1_moment_to_distro(:))
              do j=1,number_groups-1
                 loginterface_distroj(j) = logdistro(j) + &
                      (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
                      (logdistro(j+1) - logdistro(j))/ &
                      (nulibtable_logenergies(j+1)-nulibtable_logenergies(j))
                 interface_distroj(j) = exp(loginterface_distroj(j))
              enddo
              j=number_groups
              loginterface_distroj(j) = logdistro(j) + &
                   (nulibtable_logetop(j)-nulibtable_logenergies(j))* &
                   (logdistro(j) - logdistro(j-1))/ &
                   (nulibtable_logenergies(j)-nulibtable_logenergies(j-1))
              interface_distroj(j) = exp(loginterface_distroj(j))

              !find xi_i's
              xi(1) = 1.0d0
              do j=2,number_groups
                 xi(j) = interface_distroj(j)/(interface_distroj(j)+interface_distroj(j-1))
              enddo

              !unlike Mueller, there is no \Delta \epsilon because my
              !NLsolve_x are NOT per MeV, which is what this \Delta
              !\epsilon changed to term to
              do j=1,number_groups-1
                 FL(j) = velocity_center(j,1)*NLsolve_x(j)*xi(j)/ &
                      (1.0d0-nulibtable_energies(j)/nulibtable_energies(j+1)) 
                 FL(j+number_groups) = velocity_center(j,2)*NLsolve_x(j)*xi(j)/ &
                      (1.0d0-nulibtable_energies(j)/nulibtable_energies(j+1)) 
                 dFLdx(j) = FL(j)/NLsolve_x(j)
                 dFLdx(j+number_groups) = FL(j+number_groups)/NLsolve_x(j)
              enddo
              j=number_groups
              calculate_enext = nulibtable_etop(j) + &
                   0.5d0*(nulibtable_etop(j)-nulibtable_etop(j-1))**2/(nulibtable_etop(j-1)-nulibtable_etop(j-2))
              FL(j) = velocity_center(j,1)*NLsolve_x(j)*xi(j)/(1.0d0+nulibtable_energies(j)/calculate_enext)
              FL(j+number_groups) = velocity_center(j,2)*NLsolve_x(j)*xi(j)/ &
                   (1.0d0+nulibtable_energies(j)/calculate_enext)
              dFLdx(j) = FL(j)/NLsolve_x(j)
              dFLdx(j+number_groups) = FL(j+number_groups)/NLsolve_x(j)

              FR(1) = 0.0d0!never used
              dFRdx(1) = 0.0d0!never used 
              do j=2,number_groups
                 FR(j) = velocity_center(j,1)*NLsolve_x(j)*(1.0d0-xi(j))/ &
                      (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0) 
                 FR(j+number_groups) = velocity_center(j,2)*NLsolve_x(j)*(1.0d0-xi(j))/ &
                      (nulibtable_energies(j)/nulibtable_energies(j-1)-1.0d0) 
                 dFRdx(j) = FR(j)/NLsolve_x(j)
                 dFRdx(j+number_groups) = FR(j+number_groups)/NLsolve_x(j)

              enddo
              
              !no Delta epsilon on the bottom because we must convert
              !back to NOT per MeV by multiplying my \Delta \epsilon
              j=1
              NLenergyfluxterms(j) = FL(j) + FR(j+1)
              NLenergyfluxterms(j+number_groups) = FL(j+number_groups)+FR(j+1+number_groups)
              dNLenergyfluxtermsdx(j,1) = 0.0d0
              dNLenergyfluxtermsdx(j,2) = dFLdx(j)
              dNLenergyfluxtermsdx(j,3) = dFRdx(j+1)
              dNLenergyfluxtermsdx(j+number_groups,1) = 0.0d0
              dNLenergyfluxtermsdx(j+number_groups,2) = dFLdx(j+number_groups)
              dNLenergyfluxtermsdx(j+number_groups,3) = dFRdx(j+1+number_groups)

              if (trouble_brewing) then
                 dNLenergyfluxtermsdx(j,2) = (FL(j) + FR(j+1))/NLsolve_x(j)
                 dNLenergyfluxtermsdx(j,3) = 0.0d0
                 dNLenergyfluxtermsdx(j+number_groups,2) = (FL(j+number_groups) + &
                      FR(j+number_groups+1))/NLsolve_x(j)
                 dNLenergyfluxtermsdx(j+number_groups,3) = 0.0d0
              endif

              do j=2,number_groups-1
                 NLenergyfluxterms(j) = FL(j)+FR(j+1)-FL(j-1)-FR(j)
                 NLenergyfluxterms(j+number_groups) = FL(j+number_groups)+FR(j+1+number_groups)- &
                      FL(j-1+number_groups)-FR(j+number_groups)
                 dNLenergyfluxtermsdx(j,1) = -dFLdx(j-1)
                 dNLenergyfluxtermsdx(j,2) = dFLdx(j)-dFRdx(j)
                 dNLenergyfluxtermsdx(j,3) = dFRdx(j+1)
                 dNLenergyfluxtermsdx(j+number_groups,1) = -dFLdx(j-1+number_groups)
                 dNLenergyfluxtermsdx(j+number_groups,2) = dFLdx(j+number_groups)-dFRdx(j+number_groups)
                 dNLenergyfluxtermsdx(j+number_groups,3) = dFRdx(j+1+number_groups)
                 if (j.eq.2.and.trouble_brewing) then
                    dNLenergyfluxtermsdx(j,1) = 0.0d0
                    dNLenergyfluxtermsdx(j,2) = -FL(j-1)/NLsolve_x(j)+dFLdx(j)-dFRdx(j)
                    dNLenergyfluxtermsdx(j+number_groups,1) = 0.0d0
                    dNLenergyfluxtermsdx(j+number_groups,2) = -FL(j+number_groups-1)/NLsolve_x(j+number_groups)+ &
                         dFLdx(j+number_groups)-dFRdx(j+number_groups)
                 endif
              enddo

              j=number_groups
              NLenergyfluxterms(j) = FL(j)-FL(j-1)-FR(j) !assume FR(j+1)=0.0
              NLenergyfluxterms(j+number_groups) = FL(j+number_groups)-FL(j-1+number_groups)- &
                   FR(j+number_groups)
              dNLenergyfluxtermsdx(j,1) = -dFLdx(j-1)
              dNLenergyfluxtermsdx(j,2) = dFLdx(j)-dFRdx(j)
              dNLenergyfluxtermsdx(j,3) = 0.0d0
              dNLenergyfluxtermsdx(j+number_groups,1) = -dFLdx(j-1+number_groups)
              dNLenergyfluxtermsdx(j+number_groups,2) = dFLdx(j+number_groups)-dFRdx(j+number_groups)
              dNLenergyfluxtermsdx(j+number_groups,3) = 0.0d0
             
           endif

           !setup up redisual function and jacobian for source terms
           !and energy coupling.  RF is computed above for the any
           !implicit scattering/pair production innermost zones
           do j=1,number_groups

              RF(j) = RF(j) + (sourceS(j,1)+sourceG(j,1))*NLsolve_x(j) + &
                   (sourceS(j,2)+sourceG(j,2))*NLsolve_x(j+number_groups) - &
                   sourceS(j,3) - sourceG(j,3) + &
                   implicit_factor*dts*NLenergyfluxterms(j)

              RF(j+number_groups) = RF(j+number_groups) + &
                   (sourceS(j+number_groups,1) + sourceG(j+number_groups,1))*NLsolve_x(j) + &
                   (sourceS(j+number_groups,2) + sourceG(j+number_groups,2))* &
                   NLsolve_x(j+number_groups) - sourceS(j+number_groups,3) - &
                   sourceG(j+number_groups,3) + &
                   implicit_factor*dts*NLenergyfluxterms(j+number_groups)

              NL_jacobian(j,j) = NL_jacobian(j,j) + sourceS(j,1) + sourceG(j,1)
              NL_jacobian(j,j+number_groups) = NL_jacobian(j,j+number_groups) + &
                   sourceS(j,2) + sourceG(j,2)

              !add in energy couple implicit terms
              if (j.ne.1) then
                 NL_jacobian(j,j-1) = NL_jacobian(j,j-1) + &
                      implicit_factor*dts*dNLenergyfluxtermsdx(j,1)
                 NL_jacobian(j+number_groups,j-1) = &
                      NL_jacobian(j+number_groups,j-1) + &
                      implicit_factor*dts*dNLenergyfluxtermsdx(j+number_groups,1)
              endif

              NL_jacobian(j,j) = NL_jacobian(j,j) + &
                   implicit_factor*dts*dNLenergyfluxtermsdx(j,2)
              NL_jacobian(j+number_groups,j) = &
                   NL_jacobian(j+number_groups,j) + &
                   implicit_factor*dts*dNLenergyfluxtermsdx(j+number_groups,2)
              if (j.ne.number_groups) then
                 NL_jacobian(j,j+1) = NL_jacobian(j,j+1) + &
                      implicit_factor*dts*dNLenergyfluxtermsdx(j,3)
                 NL_jacobian(j+number_groups,j+1) = &
                      NL_jacobian(j+number_groups,j+1) + &
                      implicit_factor*dts*dNLenergyfluxtermsdx(j+number_groups,3)
              endif

              NL_jacobian(j+number_groups,j) = NL_jacobian(j+number_groups,j) + &
                   sourceS(j+number_groups,1) + sourceG(j+number_groups,1)
              NL_jacobian(j+number_groups,j+number_groups) = &
                   NL_jacobian(j+number_groups,j+number_groups) + &
                   sourceS(j+number_groups,2) + sourceG(j+number_groups,2)

           enddo

           !precondition assuming the dominant terms on the diagonals
           inverse = 0.0d0
           do j=1,number_groups
              det = NL_jacobian(j,j)*NL_jacobian(j+number_groups,j+number_groups) - &                        
                   NL_jacobian(j,j+number_groups)*NL_jacobian(j+number_groups,j)
              
              inverse(j,j) = NL_jacobian(j+number_groups,j+number_groups)/det
              inverse(j+number_groups,j+number_groups) = NL_jacobian(j,j)/det
              inverse(j,j+number_groups) = -NL_jacobian(j,j+number_groups)/det
              inverse(j+number_groups,j) = -NL_jacobian(j+number_groups,j)/det
           enddo

           !if things are implicit beyond E and F coupling in the same
           !group, then must use the full formalism, otherwise, do the
           !simple way, this can make a big difference if number of
           !groups is large
           new_NL_jacobian = 0.0d0
           
           if (include_Ielectron_imp.or.include_energycoupling_imp.or.include_epannihil_kernels) then
              do j=1,2*number_groups
                 do j_prime=1,2*number_groups
                    new_NL_jacobian(j,j_prime) = sum(inverse(j,:)*NL_jacobian(:,j_prime))
                 enddo
                 new_RF(j) = sum(inverse(j,:)*RF(:))
              enddo
           else
              do j=1,number_groups
                 new_NL_jacobian(j,j) = inverse(j,j)*NL_jacobian(j,j) + &
                      inverse(j,j+number_groups)*NL_jacobian(j+number_groups,j)
                 new_NL_jacobian(j,j+number_groups) = inverse(j,j)*NL_jacobian(j,j+number_groups) + &
                      inverse(j,j+number_groups)*NL_jacobian(j+number_groups,j+number_groups)

                 new_NL_jacobian(j+number_groups,j) = inverse(j+number_groups,j)*NL_jacobian(j,j) + &
                      inverse(j+number_groups,j+number_groups)*NL_jacobian(j+number_groups,j)
                 new_NL_jacobian(j+number_groups,j) = inverse(j+number_groups,j)*NL_jacobian(j,j) + &
                      inverse(j+number_groups,j+number_groups)*NL_jacobian(j+number_groups,j)
                 
                 new_RF(j) = inverse(j,j)*RF(j)+inverse(j,j+number_groups)*RF(j+number_groups)
                 new_RF(j+number_groups) = inverse(j+number_groups,j)*RF(j) + &
                      inverse(j+number_groups,j+number_groups)*RF(j+number_groups)
              enddo
           endif

           old_jacobian = NL_jacobian
           old_RF = RF
           NL_jacobian = new_NL_jacobian
           RF = new_RF

           oldx = NLsolve_x
           !RF and NL-jacobian setup, now invert to get dx
           
#if HAVE_LAPACK
           RF = -RF
           call dgesv(2*number_groups,1,NL_jacobian,2*number_groups,pivot,RF,2*number_groups,info)
#else
           stop "You need to have matrix inversion software" 
#endif

           if (sum(RF).ne.sum(RF)) then
              write(*,*) k,i,nt
              write(*,*) NL_jacobian
              write(*,*) 
              write(*,*) new_NL_jacobian
              write(*,*)
              write(*,*) old_jacobian
              write(*,*)
              write(*,*) oldx
              write(*,*) 
              write(*,*) count, "RF:", RF
              write(*,*) 
              write(*,*) count, "old RF:",old_RF
              write(*,*) 
              write(*,*) count, "old RF, post pre:",new_RF
              write(*,*) 
              write(*,*) count, "NL:", NLsolve_x!+RF                 
              stop
           endif
        
           NLsolve_x = NLsolve_x + RF

           if (minval(NLsolve_x(1:number_groups)).lt.0.0d0) then
              myloc = minloc(NLsolve_x(1:number_groups))
              if (problem_zone.ne.0) then
                 !set we have a new problem.
                 problem_fixing = .false.
              endif
              if (problem_fixing) then
                 write(*,*) "problem zone ", problem_zone,k,i,nt
                 write(*,*) "a = ",new_NL_jacobian(problem_zone,:)
                 write(*,*) "c = ",new_NL_jacobian(problem_zone+number_groups,:)
                 write(*,*) "original RF(PZ) = ",old_RF(problem_zone)
                 write(*,*) "original RF(PZ+ng) = ",old_RF(problem_zone+number_groups)
                 write(*,*) "explicit en flux, B,C,D", B_M1(k,i,problem_zone,1), &
                      C_M1(k,i,problem_zone,1),D_M1(k,i,problem_zone,1), &
                      B_M1(k,i,problem_zone,1)+C_M1(k,i,problem_zone,1)+ &
                      D_M1(k,i,problem_zone,1)
                 write(*,*) "explicit flux flux, B,C,D", B_M1(k,i,problem_zone,2), &
                      C_M1(k,i,problem_zone,2),D_M1(k,i,problem_zone,2), &
                      B_M1(k,i,problem_zone,2)+C_M1(k,i,problem_zone,2)+ &
                      D_M1(k,i,problem_zone,2)
                 write(*,*) "Se = ",new_RF(problem_zone)
                 write(*,*) "Sf = ",new_RF(problem_zone+number_groups)
                 
                 write(*,*) oldx(problem_zone),oldx(problem_zone+number_groups)
                 write(*,*) NLsolve_x(problem_zone), &
                      NLsolve_x(problem_zone)-RF(problem_zone),RF(problem_zone)
                 write(*,*) NLsolve_x(problem_zone+number_groups), &
                      NLsolve_x(problem_zone+number_groups)-RF(problem_zone+number_groups), &
                      RF(problem_zone+number_groups)
                 
                 stop "can't fix it"

              endif

              problem_fixing = .true.
              
              do j=1,number_groups
                 if (NLsolve_x(j).lt.0.0d0) then
                    problem_zone = j
                    !here's a rough sence of how the terms empirically
                    !go at the origin.  The energy coupling terms, for
                    !the energy equation, are equally dominated by the
                    !pressure tensor terms
                    !(prr=pthetatheta=pphiphi). Therefore, if problems
                    !arise with a given energy bin, assume the energy
                    !flux is equal to the (F/E^(n)) * E^(N+1),
                    !i.e. make it an implicit term, add it to
                    !sourceG(:,1) and substract from sourceG(:,3) (but
                    !note these are on different sides of the equation
                    !we are solving).  The spatial flux for the energy
                    !will be dominated by F
                    if (nothappenyet2) then
                       write(*,*) "transferring terms to implicit",i,k,j,eddy(j),eddytt(j)
                       nothappenyet2 = .false.
                    endif

                    !cases, low energy, energy coupling is large
                    !because difference between neighbouring bins is
                    !large, you can end up sending too much flux up in
                    !energy then is in the bin (hence, why we are
                    !here) Solution, make energy coupling term implicit.
                    if (j.le.3) then
                       sourceG(problem_zone,3) = q_M1_old(k,i,problem_zone,1) + &
                            B_M1(k,i,j,1) + D_M1(k,i,j,1)
                       sourceG(problem_zone,1) = sourceG(problem_zone,1) - &
                         C_M1(k,i,problem_zone,1)/q_M1_old(k,i,problem_zone,1)
                       trouble_brewing = .true.
                    else if (j.ge.number_groups-3) then
                       !sometimes the highest energy bin want to send
                       !out more energy flux than it has, fix that!
                       if (M1en_Exp_term(j).lt.0.0d0) then 
                          q_M1_old(k,i,j,1) = 2.0d0*abs(M1en_Exp_term(j))
                       endif
                       sourceG(j,3) = M1en_Exp_term(j) + q_M1_old(k,i,j,1)                       
                       !be wary of this, ensure energy conservation is not too violated

                    endif
                 endif
              enddo
              !reset all variables to original
              NLsolve_x(1:number_groups) = q_M1(k,i,:,1)
              NLsolve_x(number_groups+1:2*number_groups) = q_M1(k,i,:,2)

           endif

           RF = RF/NLsolve_x

           if (count.gt.90) then
              write(*,*) k,count, maxval(abs(RF)),maxloc(abs(RF))
              write(*,*) RF(14)*NLsolve_x(14),RF(32)*NLsolve_x(32)
              write(*,*) RF
              write(*,*) NLsolve_x
           endif

           if (maxval(abs(RF)).lt.1.0d-7) then
              !we have our solution
              stillneedconvergence = .false.
           endif
           
           if (count.gt.20) then
              !check if potential non convergence is because flux is
              !much smaller than energy density, round off error
              !prevents convergence
              maxRF = maxval(abs(RF))
              myloc = maxloc(abs(RF))
              if (myloc(1).gt.number_groups) then
                 !issue in flux
                 if (maxRF*abs(NLsolve_x(myloc(1)))/NLsolve_x(myloc(1)-number_groups).lt.1.0d-14) then
                    !trying to get tolerance on flux to less than machine precision
                    stillneedconvergence = .false.
                 endif
              endif

              !also check if its because lowest energy group is giving trouble
              if (myloc(1)==1.and.include_energycoupling_imp) then              
                 trouble_brewing = .true.
              endif
           endif
              

           if (count.gt.100.and.maxval(abs(RF)).lt.1.0d-5) then
              write(*,*) "warning, low tolerance after 100 iterations", k,nt,i
              stillneedconvergence = .false.
           else if (count.gt.100.and.maxval(abs(RF)).gt.1.0d-5) then
              write(*,*) "warning, no tolerance after 100 iterations", k
              stop

              !it could be that the explicit flux calculation is
              !strongly balenced by a geometric source term (at or
              !near the center), test for this, if true then reset
              !explicit flux source term to give convergence
           endif
           
           count = count + 1

        enddo
        problem_fixing = .false.
        trouble_brewing = .false.
        problem_zone = 0

        !set and check new neutrino variables.
        do j=1,number_groups
           q_M1(k,i,j,1) = NLsolve_x(j)
           q_M1(k,i,j,2) = NLsolve_x(j+number_groups)

           if (q_M1(k,i,j,1).le.0.0d0) then
              write(*,*) k,i,j,nt
              stop "do_implicit_step: new en RHS is < 0.0"
           endif

           if (abs(q_M1(k,i,j,2)/oneX).gt.q_M1(k,i,j,1)) then
              if (nothappenyet1.and.k.lt.M1_imaxradii) then
                 nothappenyet1 = .false.
                 !this will happen a lot.  output is supressed to at most once per 10 time steps.
                 if (mod(nt,10).eq.0) write(*,*) "warning: do_implicit_step: flux>en",i,j,k,nt,q_M1(k,i,j,2)/oneX,q_M1(k,i,j,1)
              endif
              !fix it
              q_M1(k,i,j,2) = oneX*q_M1(k,i,j,2) / abs((1.0d0+1.0d-8)*q_M1(k,i,j,2)/q_M1(k,i,j,1))
              
           endif

           if (q_M1(k,i,j,1).ne.q_M1(k,i,j,1)) then
              write(*,*) k,i,j,nt
              stop "do_implicit_step: new en RHS is NaNing"
           endif

           if (q_M1(k,i,j,2).ne.q_M1(k,i,j,2)) then
              write(*,*) k,i,j
              stop "do_implicit_step: new flux RHS is NaNing"
           endif

           q_M1_fluid(k,i,j,1) = q_M1(k,i,j,1)*W2 - &
                2.0d0*q_M1(k,i,j,2)*W2*onev*invX + &
                q_M1(k,i,j,3)*q_M1(k,i,j,1)*W2*v2*invX2
           
           q_M1_fluid(k,i,j,2) = -(q_M1(k,i,j,1)*oneW - &
                q_M1(k,i,j,2)*oneW*onev/oneX)*W2*onev*invX + &
                W2*oneW*q_M1(k,i,j,2)*invX2 - &
                q_M1(k,i,j,3)*q_M1(k,i,j,1)*invX2**2*W2*oneW*onev*oneX

        enddo

        !determine matter source terms
        Sr = 0.0d0
        Stnalpha = 0.0d0
        Stnum = 0.0d0
        do j=1,number_groups

           oneM1en =  q_M1(k,i,j,1)
           oneM1flux =  q_M1(k,i,j,2)
           oneeddy = eddy(j)

           !note, be careful with signs and alphas here
           Stzone = (-sourceS(j,1)*oneM1en - sourceS(j,2)*oneM1flux + &
                sourceS(j,3))/(implicit_factor*dts*alp2)
           
           if (include_Ielectron_imp.or.include_Ielectron_exp) then
              Stnalpha = Stnalpha - onealp*ies_sourceterm(k,i,j,1)
           endif
           if (include_epannihil_kernels) then
              Stnalpha = Stnalpha - onealp*epannihil_sourceterm(k,i,j,1)
           endif
           Stnalpha = Stnalpha - Stzone*onealp
           
           Srzone = (-sourceS(j+number_groups,1)*oneM1en - &
                sourceS(j+number_groups,2)*oneM1flux + &
                sourceS(j+number_groups,3))/(implicit_factor*dts*onealp*X2)
              
           if (include_Ielectron_imp.or.include_Ielectron_exp) then
              Sr = Sr + ies_sourceterm(k,i,j,2)
           endif
           if (include_epannihil_kernels) then
              Sr = Sr + epannihil_sourceterm(k,i,j,2)
           endif
           Sr = Sr + Srzone
           
           if (i.eq.1) sign_one = 1.0d0
           if (i.eq.2) sign_one = -1.0d0
           if (i.gt.2) sign_one = 0.0d0
              
           !find source term in fluid frame / energy
           Stnum = Stnum - sign_one*(eas(k,i,j,1) - eas(k,i,j,2)*q_M1_fluid(k,i,j,1))*nulibtable_inv_energies(j)

           !$OMP CRITICAL
           ynu(k) = ynu(k) + sign_one*q_M1_fluid(k,i,j,1)*4.0d0*pi/rho(k)* &
                nulibtable_inv_energies(j)*(amu_cgs*mass_gf)

           press_nu(k) = press_nu(k) + oneeddy*oneM1en*4.0d0*pi*invX2**2
           energy_nu(k) = energy_nu(k) + oneM1en*4.0d0*pi
           !$OMP END CRITICAL        

        enddo
           
        !$OMP CRITICAL
        !momentum term, units of reduced(energy/cm^2/s^2/srad)?
        !Shibata 7.7, kinda, there variable is slight different
        M1_matter_source(k,2) = M1_matter_source(k,2) - onealp*oneX*Sr
        !energy term, units of reduced(energy/cm^3/s/srad), Shibata 7.8
        M1_matter_source(k,3) = M1_matter_source(k,3) + onealp*Stnalpha
        !ye term, units of reduced(#/cm^3/s/srad), Shibata 7.9
        M1_matter_source(k,4) =  M1_matter_source(k,4) + onealp*Stnum
        !$OMP END CRITICAL

     enddo
  enddo
  !$OMP END PARALLEL DO! end do

end subroutine M1_implicitstep
