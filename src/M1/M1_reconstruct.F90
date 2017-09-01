!-*-f90-*- 
!This routine takes care of the spatial reconstruction of
!the first two neutrino moments.

subroutine M1_reconstruct

  use GR1D_module
  use ppm
  implicit none

  integer :: i,j,k
  real*8 :: M1en_space(n1),M1flux_space(n1)
  real*8 :: M1en_space_plus(n1),M1flux_space_plus(n1)
  real*8 :: M1en_space_minus(n1),M1flux_space_minus(n1)
  real*8 :: M1_testcase_travelling_pulse,M1_testcase_diffusion_wave !function calls for test case boundary conditions

  !$OMP PARALLEL DO PRIVATE(i,k,M1en_space,M1en_space_plus,M1en_space_minus, &
  !$OMP M1flux_space,M1flux_space_plus,M1flux_space_minus) 
  do j=1,number_groups
     do i=1,number_species_to_evolve

        M1en_space = q_M1(:,i,j,1)
        M1flux_space = q_M1(:,i,j,2)

        !boundaries
        do k=1,ghosts1
           M1en_space(ghosts1+1-k) = M1en_space(ghosts1+k)
           M1flux_space(ghosts1+1-k) = -M1flux_space(ghosts1+k)
           M1en_space(M1_imaxradii+k) = M1en_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1(M1_imaxradii+k))**2
           M1flux_space(M1_imaxradii+k) = M1flux_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1(M1_imaxradii+k))**2
        enddo

        !some testcases have their own boundary conditions.
        if (M1_testcase_number.eq.3) then
           do k=1,ghosts1
              M1en_space(ghosts1+1-k) = M1_testcase_travelling_pulse(x1(ghosts1+1-k),time*time_gf,0)
              q_M1(ghosts1+1-k,i,j,1) = M1en_space(ghosts1+1-k)
              M1flux_space(ghosts1+1-k) = M1_testcase_travelling_pulse(x1(ghosts1+1-k),time*time_gf,1)
              q_M1(ghosts1+1-k,i,j,2) = M1flux_space(ghosts1+1-k)
              M1en_space(M1_imaxradii+k) = M1_testcase_travelling_pulse(x1(M1_imaxradii+k),time*time_gf,0)
              q_M1(M1_imaxradii+k,i,j,1) = M1en_space(M1_imaxradii+k)
              M1flux_space(M1_imaxradii+k) = M1_testcase_travelling_pulse(x1(M1_imaxradii+k),time*time_gf,1)
              q_M1(M1_imaxradii+k,i,j,2) = M1flux_space(M1_imaxradii+k)
           enddo
        else if (M1_testcase_number.eq.8.or.M1_testcase_number.eq.9) then
           do k=1,ghosts1
              M1en_space(ghosts1+1-k) = M1_testcase_diffusion_wave(x1(ghosts1+1-k),time*time_gf, &
                   eas(ghosts1+1-k,i,j,2),0,M1_testcase_number)
              q_M1(ghosts1+1-k,i,j,1) = M1en_space(ghosts1+1-k)
              M1flux_space(ghosts1+1-k) = M1_testcase_diffusion_wave(x1(ghosts1+1-k),time*time_gf, &
                   eas(ghosts1+1-k,i,j,2),1,M1_testcase_number)
              q_M1(ghosts1+1-k,i,j,2) = M1flux_space(ghosts1+1-k)
              M1en_space(M1_imaxradii+k) = M1_testcase_diffusion_wave(x1(M1_imaxradii+k),time*time_gf, &
                   eas(M1_imaxradii+k,i,j,2),0,M1_testcase_number)
              q_M1(M1_imaxradii+k,i,j,1) = M1en_space(M1_imaxradii+k)
              M1flux_space(M1_imaxradii+k) = M1_testcase_diffusion_wave(x1(M1_imaxradii+k),time*time_gf, &
                   eas(M1_imaxradii+k,i,j,2),1,M1_testcase_number)
              q_M1(M1_imaxradii+k,i,j,1) = M1flux_space(M1_imaxradii+k)
           enddo
        endif

        !set F to F/E
        M1flux_space = M1flux_space/M1en_space

        !reconstruct
        if (M1_reconstruction_method.eq.'tvd') then
           call tvd_reconstruction(n1,ghosts1,M1en_space,M1en_space_plus,M1en_space_minus,'minmod')
           call tvd_reconstruction(n1,ghosts1,M1flux_space,M1flux_space_plus,M1flux_space_minus,'minmod') 
        else if (M1_reconstruction_method.eq.'ppm') then
           call ppm_interpolate(M1en_space,M1en_space_plus,M1en_space_minus)
           call ppm_interpolate(M1flux_space,M1flux_space_plus,M1flux_space_minus)
           call ppm_monotonize(M1en_space,M1en_space_plus,M1en_space_minus)
           call ppm_monotonize(M1flux_space,M1flux_space_plus,M1flux_space_minus)       
        else
           stop "Implement this reconstruction method in M1"
        endif

        !reset F/E to F
        M1flux_space = M1flux_space*M1en_space
        M1flux_space_plus = M1flux_space_plus*M1en_space_plus
        M1flux_space_minus = M1flux_space_minus*M1en_space_minus

        !boundaries for reconstructed variables
        do k=1,ghosts1
           M1en_space_plus(ghosts1+1-k) = M1en_space_minus(ghosts1+k)
           M1flux_space_plus(ghosts1+1-k) = -M1flux_space_minus(ghosts1+k)
           M1en_space_minus(ghosts1+1-k) = M1en_space_plus(ghosts1+k)
           M1flux_space_minus(ghosts1+1-k) = -M1flux_space_plus(ghosts1+k)

           M1en_space_plus(M1_imaxradii+k-1) = M1en_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1i(M1_imaxradii+k))**2
           M1flux_space_plus(M1_imaxradii+k-1) = M1flux_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1i(M1_imaxradii+k))**2
           M1en_space_minus(M1_imaxradii+k) = M1en_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1i(M1_imaxradii+k))**2
           M1flux_space_minus(M1_imaxradii+k) = M1flux_space(M1_imaxradii)* &
                (x1(M1_imaxradii)/x1i(M1_imaxradii+k))**2

        enddo

        !some testcases have their own boundary conditions.
        if (M1_testcase_number.eq.3) then
           do k=1,ghosts1
              M1en_space_plus(ghosts1+1-k) = M1_testcase_travelling_pulse(x1i(ghosts1+1-k+1),time*time_gf,0)
              M1flux_space_plus(ghosts1+1-k) = M1_testcase_travelling_pulse(x1i(ghosts1+1-k+1),time*time_gf,1)
              M1en_space_minus(ghosts1+1-k) = M1_testcase_travelling_pulse(x1i(ghosts1+1-k),time*time_gf,0)
              M1flux_space_minus(ghosts1+1-k) = M1_testcase_travelling_pulse(x1i(ghosts1+1-k),time*time_gf,1)
              M1en_space_plus(M1_imaxradii+k-1) = M1_testcase_travelling_pulse(x1i(M1_imaxradii+k),time*time_gf,0)
              M1flux_space_plus(M1_imaxradii+k-1) = M1_testcase_travelling_pulse(x1i(M1_imaxradii+k),time*time_gf,1)
              M1en_space_minus(M1_imaxradii+k) = M1_testcase_travelling_pulse(x1i(M1_imaxradii+k),time*time_gf,0)
              M1flux_space_minus(M1_imaxradii+k) = M1_testcase_travelling_pulse(x1i(M1_imaxradii+k),time*time_gf,1)
           enddo
        else if (M1_testcase_number.eq.8.or.M1_testcase_number.eq.9) then
           do k=1,ghosts1
              M1en_space_plus(ghosts1+1-k) = M1_testcase_diffusion_wave(x1i(ghosts1+1-k+1), &
                   time*time_gf,eas(ghosts1+1-k+1,i,j,2),0,M1_testcase_number)
              M1flux_space_plus(ghosts1+1-k) = M1_testcase_diffusion_wave(x1i(ghosts1+1-k+1), &
                   time*time_gf,eas(ghosts1+1-k+1,i,j,2),1,M1_testcase_number)
              M1en_space_minus(ghosts1+1-k) = M1_testcase_diffusion_wave(x1i(ghosts1+1-k), &
                   time*time_gf,eas(ghosts1+1-k,i,j,2),0,M1_testcase_number)
              M1flux_space_minus(ghosts1+1-k) = M1_testcase_diffusion_wave(x1i(ghosts1+1-k), &
                   time*time_gf,eas(ghosts1+1-k,i,j,2),1,M1_testcase_number)
              M1en_space_plus(M1_imaxradii+k-1) = M1_testcase_diffusion_wave(x1i(M1_imaxradii+k), &
                   time*time_gf,eas(M1_imaxradii+k,i,j,2),0,M1_testcase_number)
              M1flux_space_plus(M1_imaxradii+k-1) = M1_testcase_diffusion_wave(x1i(M1_imaxradii+k), &
                   time*time_gf,eas(M1_imaxradii+k,i,j,2),1,M1_testcase_number)
              M1en_space_minus(M1_imaxradii+k) = M1_testcase_diffusion_wave(x1i(M1_imaxradii+k), &
                   time*time_gf,eas(M1_imaxradii+k,i,j,2),0,M1_testcase_number)
              M1flux_space_minus(M1_imaxradii+k) = M1_testcase_diffusion_wave(x1i(M1_imaxradii+k), &
                   time*time_gf,eas(M1_imaxradii+k,i,j,2),1,M1_testcase_number)
           enddo
        endif

        !set reconstructed values for later
        q_M1p(:,i,j,1,1) = M1en_space_plus
        q_M1m(:,i,j,1,1) = M1en_space_minus

        q_M1p(:,i,j,2,1) = M1flux_space_plus
        q_M1m(:,i,j,2,1) = M1flux_space_minus

     enddo
  enddo
  !$OMP END PARALLEL DO! end do

end subroutine M1_reconstruct
