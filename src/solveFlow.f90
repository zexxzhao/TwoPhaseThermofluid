!======================================================================
!
!======================================================================
subroutine solmultiphasethermofluid_stag(istep)

  use aAdjKeep
  use mpi
  use commonvars

  implicit none

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

    call solveNavSto(inewt, momres0, conres0, NSConverged, lsres0, LSConverged)

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
