!======================================================================
!
!======================================================================
subroutine solmultiphasethermofluid_stag(config, istep)

  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  use configuration
  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: istep

  integer :: inewt, i, NL_max
  ! real(8) :: momres0, conres0, convres0, meshres0, lsres0
  real(8) :: residual0(4), residual(4)
  integer :: converged(4)

  real(8) :: sol(NNODE, 6)
  real(8) :: t1, t2
  real(8) :: utol(4)
  
  ! utol = (/NS_NL_UTOL, NS_NL_PTOL, LSC_NL_TOL, TEM_NL_TOL/)
  utol(:) = config%newton_raphson%rtol(:)
  inewt = 0
  ! momres0 = -1.0d0
  ! conres0 = -1.0d0
  ! convres0 = -1.0d0
  ! lsres0 = -1.0d0

  residual(:) = 0d0
  residual0(:) = 0d0
  converged(:) = 0

  IBC(:, :) = 0
  call setBCs_NSVOF()
  call setBCs_Tem()
  !call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
  !                          ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
  call assembleQuenching(config, ASSEMBLE_TENSOR_VEC, &
                         ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

  ! write(*,*) "myid=", myid, "GET RHS0"
  call calculate_residual(residual0, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
  residual(:) = residual0(:)
  call check_convergence(converged, residual, residual0, utol, &
                         ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
  call print_residual(residual, residual0, utol, &
                      ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                      inewt)

  if(istep <= -1) then 
    NL_max = 1
  else 
    NL_max = NS_NL_itermax
  endif
  do inewt = 1, NL_max

    !---------------------------
    ! Solve NavStoVOF
    !---------------------------
    ! IBC(:, :) = 0
    ! call setBCs_NSVOF()
    ! IBC(:, 5) = 1
    ! IBC(:, 6:8) = 1
    ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
    !                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    call assembleQuenching(config, ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    !call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
    !call check_convergence(converged, residual, residual0, utol, &
    !                      ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    !call print_residual(residual, residual0, utol, &
    !                    ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF, &
    !                    inewt)

    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    sol(:, 1:5) = 0d0
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, config%ksp%max_iter(1), &
                        config%ksp%min_iter(1), &
                        NNODE, maxNSHL, NSD, &
                        sol(:, 5), lhsLS, &
                        lhsLSu, lhsUls, lhsPls, rhsGls)

    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve NSVOF GMRES:", t2 - t1, "seconds"
    end if
    acg = acg + sol(:, 1:3)
    ug = ug + gami*Delt*sol(:, 1:3)
    pg = pg + alfi*gami*Delt*sol(:, 4)

    ! IBC(:, :) = 0
    ! call setBCs_NSVOF()
    ! IBC(:, 1:4) = 1
    ! IBC(:, 6:8) = 1

    ! call assembleQuenching(config, ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
    !                        ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    ! sol(:, 1:5) = 0d0
    ! call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
    !                     rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
    !                     lhsK11, lhsG, lhsD1, lhsM, icnt, &
    !                     NS_GMRES_tol, config%ksp%max_iter(1), &
    !                     config%ksp%min_iter(1), &
    !                     NNODE, maxNSHL, NSD, &
    !                     sol(:, 5), lhsLS, &
    !                     lhsLSu, lhsUls, lhsPls, rhsGls)

    ! if(sum(sol(:, 1:4)**2) > 0) then
    !   write(*,*) "BC NS not satisfied", myid, "WWWWWW"
    ! endif
    rphig = rphig + sol(:, 5)
    phig = phig + gami*Delt*sol(:, 5)

    !-----------------------------
    ! Solve Temperature
    !-----------------------------
    if (istep > 0) then
      ! IBC(:, :) = 0
      ! call setBCs_Tem()
      ! IBC(:, 1:5) = 1
      ! IBC(:, 7:8) = 1
      ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
      !                           ASSEMBLE_FIELD_TEM)
      call assembleQuenching(config, ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                             ASSEMBLE_FIELD_TEM)
      call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
      call check_convergence(converged, residual, residual0, utol, &
                             ASSEMBLE_FIELD_TEM)
      call print_residual(residual, residual0, utol, &
                          ASSEMBLE_FIELD_TEM, &
                          inewt)

      !write(*,*) "Residual T = ", residual0(4), residual(4)
      if (myid .eq. 0) then
        call CPU_TIME(t1)
      end if

      sol(:, 6) = 0d0
      call SparseGMRES_tem(LHStem, config%ksp%rtol(4), col, row, &
                           rhsgtem, sol(:, 6), config%ksp%max_iter(4), config%ksp%min_iter(4), &
                           NNODE, maxNSHL, icnt, NSD)
      if (myid .eq. 0) then
        call CPU_TIME(t2)
        write (*, *) "Total time solve TEM GMRES:", t2 - t1, "seconds"
      end if
      !write(*,*) "solT =", sum(sol(:, 6)**2) 
      rTg = rTg + sol(:, 6)
      Tg = Tg + gami*Delt*sol(:, 6)

    end if
    ! IBC(:, :) = 0
    ! call setBCs_NSVOF()
    ! call setBCs_Tem()
    ! call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
    !                        ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
    call assembleQuenching(config, ASSEMBLE_TENSOR_VEC, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

    call calculate_residual(residual, RHSGu, RHSGp, RHSGls, RHSGtem, NNODE, NSD)
    call check_convergence(converged, residual, residual0, utol, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)
    call print_residual(residual, residual0, utol, &
                        ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                        inewt)

    if(ismaster ) write (*, *) "Convergence:", converged
    if (size(converged) == sum(converged)) exit

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  end do

end subroutine solmultiphasethermofluid_stag

