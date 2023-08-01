
subroutine calculate_residual(residual)
  use mpi
  use aAdjKeep
  implicit none
  integer, parameter :: NRES = 4

  real(8), intent(out) :: residual(4)

  real(8) :: rl(NRES)
  integer :: i

  rl(1) = sum(RHSGu(:, :)*RHSGu(:, :))
  rl(2) = sum(RHSGp(:)*RHSGp(:))
  rl(3) = sum(RHSGls(:)*RHSGls(:))
  rl(4) = sum(RHSGtem(:)*RHSGtem(:))
  if (numnodes > 1) then
    call MPI_ALLREDUCE(rl, residual, NRES, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  else
    residual(:) = rl(:)
  end if

  do i = 1, NRES
    residual(i) = sqrt(residual(i))
  end do

end subroutine calculate_residual

subroutine check_convergence(r0, r, verbose, inewt, assemble_field_flag, converged)
  use mpi
  use aAdjKeep
  use commonvars
  use commonpars
  implicit none

  integer, parameter :: NRES = 4

  real(8), intent(in) :: r0(NRES), r(NRES)
  integer, intent(in) :: inewt
  integer, intent(in) :: verbose
  integer, intent(in) :: assemble_field_flag

  integer, intent(out) :: converged(NRES)
  character(len=80) :: fomt

  converged(:) = 0

  if (r0(1)*NS_NL_UTOL > r(1)) converged(1) = 1
  if (r0(2)*NS_NL_PTOL > r(2)) converged(2) = 1
  if (r0(3)*LSC_NL_TOL > r(3)) converged(3) = 1
  if (r0(4)*LSC_NL_TOL > r(4)) converged(4) = 1
  if (ismaster .and. verbose > 0) then
    fomt = "(I3,a,x,ES13.6,x,F12.6)"
    if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
      write (*, fomt) inewt, ") Total Mom. Res. Norm = ", &
        r(1), 1.0d2*r(1)/r0(1)
      write (*, fomt) inewt, ") Continuity Res. Norm = ", &
        r(2), 1.0d2*r(2)/r0(2)
    end if

    if (iand(assemble_field_flag, ASSEMBLE_FIELD_LS) > 0) then
      write (*, fomt) inewt, ") LS Res. Norm = ", &
        r(3), 1.0d2*r(3)/r0(3)
    end if
    write (*, *)
  end if

end subroutine check_convergence

!======================================================================
!
!======================================================================
subroutine assembleNavStoVOFTem(assemble_tensor_flag, assemble_field_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: assemble_field_flag ! assemble NS + LS/VOF + Tem

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE), &
             rTgAlpha(NNODE), TgAlpha(NNODE)

  real(8) :: t1, t2
  integer :: FFlag

  !---------------------------------
  ! set Dirichlet BCs
  !---------------------------------
  call setBCs_NSVOF()

  !---------------------------------
  ! prediction stage
  !---------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  TgAlpha = Tgold + alfi*(Tg - Tgold)
  rTgAlpha = rTgold + almi*(rTg - rTgold)

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    RHSGu = 0.0d0
    RHSGp = 0.0d0
    RHSGls = 0.0d0
    RHSGtem = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    LHSK11 = 0.0d0
    LHSG = 0.0d0
    LHSD1 = 0.0d0
    LHSM = 0.0d0

    LHSLS = 0.0d0
    LHSULS = 0.0d0
    LHSLSU = 0.0d0
    LHSPLS = 0.0d0

    LHSTem = 0.0d0
  endif

  if (myid .eq. 0) then
    call CPU_TIME(t1)
  endif

  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF) > 0) then
    call IntElmAss_NSVOF(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                            TgAlpha, rTgAlpha, &
                            assemble_tensor_flag)

    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) FFLAG = 0
    call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                              acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                              FFLAG)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    call IntElmAss_Tem(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                       acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                       TgAlpha, rTgAlpha, &
                       assemble_tensor_flag)

  end if

  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    call commu(RHSGp, 1, 'in ')
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGls, 1, 'in ')
    call commu(RHSGTem, 1, 'in ')
    ! call commu(rhsgq, 1, 'in ')
    ! call commu(rhsgq, 1, 'out ')
    ! call commu(lhsgq, 1, 'in ')
    ! call commu(lhsgq, 1, 'out ')
  end if

end subroutine assembleNavStoVOFTem
!======================================================================
!
!======================================================================
! subroutine solveNavStoVOF(inewt, MomRes0, ConRes0, NSConverged, lsres0, LSConverged)
!
!   use aAdjKeep
!   use mpi
!   use commonvars
!
!   implicit none
!
!   integer, intent(in)    :: inewt
!   real(8), intent(inout) :: MomRes0, ConRes0, LSRes0
!   logical, intent(out)   :: NSConverged, LSConverged
!
!   integer :: i, j, b, p
!
!   real(8) :: dgAlpha(NNODE, NSD), &
!              ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
!              acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
!              pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE)
!
!   real(8) :: MomRes, ConRes, LSRes, MomResl, ConResl, LSResl, gfor(NSD), gmom(NSD)
!   real(8) :: tmp1(NNODE, NSD), vtmp(NSD), TRACT(NNODE, NSD)
!   real(8) :: t1, t2
!   character(len=80) :: fomt
!
!   call assembleNavStoVOF()
!   MomResl = sum(RHSGu(:, 1)**2) + sum(RHSGu(:, 2)**2) + sum(RHSGu(:, 3)**2)
!   ConResl = sum(RHSGp**2)
!   LSResl = sum(RHSGls**2)
!
!   MomRes = MomResl
!   ConRes = ConResl
!   LSRes = LSResl
!
!   if (numnodes > 1) then
!     call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
!                        MPI_SUM, MPI_COMM_WORLD, mpi_err)
!     call MPI_ALLREDUCE(ConResl, ConRes, 1, MPI_DOUBLE_PRECISION, &
!                        MPI_SUM, MPI_COMM_WORLD, mpi_err)
!     call MPI_ALLREDUCE(LSResl, LSRes, 1, MPI_DOUBLE_PRECISION, &
!                        MPI_SUM, MPI_COMM_WORLD, mpi_err)
!   end if
!
!   MomRes = sqrt(MomRes)
!   ConRes = sqrt(ConRes)
!   LSRes = sqrt(LSRes)
!
!   ! Print residual
!   if (inewt == 1) then
!     MomRes0 = MomRes
!     ConRes0 = ConRes
!     LSRes0 = LSRes
!   end if
!
!   if (ismaster) then
!     fomt = "(I3,a,x,ES13.6,x,F12.6)"
!     write (*, fomt) inewt, ") Total Mom. Res. Norm = ", &
!       MomRes, 1.0d2*MomRes/MomRes0
!
!     write (*, fomt) inewt, ") Continuity Res. Norm = ", &
!       ConRes, 1.0d2*ConRes/ConRes0
!
!     write (*, fomt) inewt, ") LS Res. Norm = ", &
!       LSRes, 1.0d2*LSRes/LSRes0
!     write (*, *)
!   end if
!
!   if (isnan(MomRes + ConRes + LSRes)) then
!     write (*, *) myid, inewt, '!=== Norm is NaN ===!'
!     stop
!   end if
!
!   ! Check convergence and solve linear system
!   if (((MomRes/MomRes0) < NS_NL_Utol) .and. &
!       ((ConRes/ConRes0) < NS_NL_Ptol) .and. &
!       ((LSRes/LSRes0) < LSC_NL_tol)) then
!     NSConverged = .true.
!     if (ismaster) then
!       write (*, *)
!       write (*, *) "    Converged: skip linear solve"
!       write (*, *)
!     end if
!
!   else
!     NSConverged = .false.
!     acgAlpha = 0.0d0; pgAlpha = 0.0d0; rphigAlpha = 0.0d0
! !    lhsLSu = 0.0d0!; lhsUls = 0.0d0; lhsPls = 0.0d0
!     if (myid .eq. 0) then
!       call CPU_TIME(t1)
!     end if
!     call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
!                         rhsGu, rhsGp, acgAlpha, pgAlpha, &
!                         lhsK11, lhsG, lhsD1, lhsM, icnt, &
!                         NS_GMRES_tol, NS_GMRES_itermax, &
!                         NS_GMRES_itermin, &
!                         NNODE, maxNSHL, NSD, &
!                         rphigAlpha, lhsLS, &
!                         lhsLSu, lhsUls, lhsPls, rhsGls)
!
!     if (myid .eq. 0) then
!       call CPU_TIME(t2)
!       write (*, *) "Total time solve:", t2 - t1, "seconds"
!     end if
!
!     acg = acg + acgAlpha
!     ug = ug + gami*Delt*acgAlpha
!     pg = pg + alfi*gami*Delt*pgAlpha
!
!     rphig = rphig + rphigAlpha
!     phig = phig + gami*Delt*rphigAlpha
!   end if
!
!   !--------------------------------------------------------------------
!   ! compute torque (weak BC/Surface Int.)
!   !--------------------------------------------------------------------
!   acgAlpha = acgold + almi*(acg - acgold)
!   acgmAlpha = acgmold + almi*(acgm - acgmold)
!   ugAlpha = ugold + alfi*(ug - ugold)
!   ugmAlpha = ugmold + alfi*(ugm - ugmold)
!   dgAlpha = dgold + alfi*(dg - dgold)
!   pgAlpha = pg
!   phigAlpha = phigold + alfi*(phig - phigold)
!   rphigAlpha = rphigold + alfi*(rphig - rphigold)
!
!   if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
!
! end subroutine solveNavStoVOF

