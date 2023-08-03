!======================================================================
!
!======================================================================
subroutine solmultiphasethermofluid_stag(istep)

  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

  integer, intent(in) :: istep

  integer :: inewt, i
  ! real(8) :: momres0, conres0, convres0, meshres0, lsres0
  real(8) :: residual0(4), residual(4)
  integer :: converged(4)
  logical :: NSConverged, LSConverged

  integer :: assemble_field_flag
  real(8) :: sol(NNODE, 6)
  real(8) :: t1, t2
  ! momres0 = -1.0d0
  ! conres0 = -1.0d0
  ! convres0 = -1.0d0
  ! lsres0 = -1.0d0

  residual(:) = 0d0
  residual0(:) = 0d0
  converged(:) = 0
  call setBCs_NSVOF()
  call setBCs_Tem()
  call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
                            ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

  call calculate_residual(residual0)
  residual(:) = residual0(:)
  call check_convergence(residual0, residual, 1, 0, &
                         ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                         converged)

  do inewt = 1, NS_NL_itermax

    !---------------------------
    ! Solve NavStoVOF
    !---------------------------
    call setBCs_NSVOF()
    IBC(:, 6:8) = 1
    call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                              ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF)
    call calculate_residual(residual)
    call check_convergence(residual0, residual, 1, inewt, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF, &
                           converged)
    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    sol(:, 1:5) = 0d0
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, sol(:, 1:3), sol(:, 4), &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
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

    rphig = rphig + sol(:, 5)
    phig = phig + gami*Delt*sol(:, 5)

    !-----------------------------
    ! Solve Temperature
    !-----------------------------
    if (istep > 10) then
      call setBCs_Tem()
      IBC(:, 1:5) = 1
      IBC(:, 7:8) = 1
      call assembleNavStoVOFTem(ASSEMBLE_TENSOR_MAT + ASSEMBLE_TENSOR_VEC, &
                                ASSEMBLE_FIELD_TEM)
      call calculate_residual(residual)
      call check_convergence(residual0, residual, 1, inewt, &
                            ASSEMBLE_FIELD_TEM, &
                            converged)
      !write(*,*) "Residual T = ", residual0(4), residual(4)
      if (myid .eq. 0) then
        call CPU_TIME(t1)
      end if

      sol(:, 6) = 0d0
      call SparseGMRES_tem(LHStem, NS_GMRES_tol, col, row, &
                           rhsgtem, sol(:, 6), NS_GMRES_itermax, NS_GMRES_itermin, &
                           NNODE, maxNSHL, icnt, NSD)
      if (myid .eq. 0) then
        call CPU_TIME(t2)
        write (*, *) "Total time solve TEM GMRES:", t2 - t1, "seconds"
      end if
      !write(*,*) "solT =", sum(sol(:, 6)**2) 
      rTg = rTg + sol(:, 6)
      Tg = Tg + gami*Delt*sol(:, 6)

    end if
    call assembleNavStoVOFTem(ASSEMBLE_TENSOR_VEC, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM)

    call calculate_residual(residual)
    call check_convergence(residual0, residual, 1, inewt, &
                           ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF + ASSEMBLE_FIELD_TEM, &
                           converged)
    if(ismaster ) write (*, *) "Convergence:", converged
    if (size(converged) == sum(converged)) exit

    if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  end do

end subroutine solmultiphasethermofluid_stag

!======================================================================
!
!======================================================================
subroutine solflow_stag(istep, SH, NM)

  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell

  implicit none

  type(shell_bld), intent(inout) :: SH
  type(shell_nmb), intent(inout) :: NM
  integer :: istep

  integer :: inewt
  real(8) :: momres0, conres0, convres0, meshres0, lsres0
  logical :: NSConverged, LSConverged

  momres0 = -1.0d0
  conres0 = -1.0d0
  convres0 = -1.0d0
  lsres0 = -1.0d0

  !------------------------------------------------
  ! Solve Navier-Stokes
  !------------------------------------------------
  do inewt = 1, NS_NL_itermax

    if (nonmatch) then
      call solveNavStoDG(inewt, momres0, conres0, NSConverged, lsres0, LSConverged, NM)
    else
      call solveNavSto(inewt, momres0, conres0, NSConverged, lsres0, LSConverged)
    end if

    if (NSConverged) exit

    if (shel) then
      call solveKLShell(SH, inewt, istep)
      call SolveMesh(inewt, meshres0)
    end if

  end do

end subroutine solflow_stag

!======================================================================
!
!======================================================================
subroutine solflow_mono(istep)

  use aAdjKeep
  use mpi
  use commonvars

  implicit none

  integer :: istep
  integer :: inewt
  real(8) :: momres0, conres0, convres0, meshres0
  real(8) :: normFb0, normMb0
  logical :: NSConverged, LSConverged, RBMConverged

  write (*, *) "Are you sure you want to call solflow_mono???"
  stop

  momres0 = -1.0d0
  conres0 = -1.0d0
  convres0 = -1.0d0

  !------------------------------------------------
  ! Get prediction move
  !------------------------------------------------
  if (move) then
    call moveRBMesh(0, normFb0, normMb0, meshres0, RBMConverged)
  else
    RBMConverged = .true.
  end if

  !------------------------------------------------

  !------------------------------------------------
  ! Newton loop
  !------------------------------------------------
  do inewt = 1, max(NS_NL_itermax, LSC_NL_itermax)
    if (mono_iter == 0) then
      call solveNavSto(inewt, momres0, conres0, NSConverged)
    else if (mono_iter == 1) then
      call solveNavStoLS(inewt, momres0, conres0, NSConverged, &
                         convres0, LSConverged)
    else
      write (*, *) "Undefined iteration procedure"
      write (*, *) "   mono_iter =", mono_iter
    end if

    if (move) call moveRBMesh(inewt, normFb0, normMb0, &
                              meshres0, RBMConverged)

    if (NSConverged .and. LSConverged .and. RBMConverged) exit

  end do

  !------------------------------------------------
  ! Move rigid body and mesh - if not converged
  !------------------------------------------------
  if (move .and. (.not. (NSConverged .and. RBMConverged))) then
    call moveRBMesh(inewt, normFb0, normMb0, meshres0, RBMConverged)
  end if

end subroutine solflow_mono
