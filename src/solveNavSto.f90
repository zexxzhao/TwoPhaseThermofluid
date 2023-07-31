!======================================================================
!
!======================================================================
subroutine solveNavSto(inewt, MomRes0, ConRes0, NSConverged, lsres0, LSConverged)

  use aAdjKeep
  use mpi
  use commonvars

  implicit none

  integer, intent(in)    :: inewt
  real(8), intent(inout) :: MomRes0, ConRes0, LSRes0
  logical, intent(out)   :: NSConverged, LSConverged

  integer :: i, j, b, p

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE)

  real(8) :: MomRes, ConRes, LSRes, MomResl, ConResl, LSResl, gfor(NSD), gmom(NSD)
  real(8) :: tmp1(NNODE, NSD), vtmp(NSD), TRACT(NNODE, NSD)
  real(8) :: t1, t2
  character(len=80) :: fomt

  call setBCs_CFD()
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  ! Assemble Navier-Stokes LHS/RHS
  RHSGu = 0.0d0
  RHSGp = 0.0d0
  RHSGls = 0.0d0

  LHSK11 = 0.0d0
  LHSG = 0.0d0
  LHSD1 = 0.0d0
  LHSM = 0.0d0

  LHSLS = 0.0d0
  LHSULS = 0.0d0
  LHSLSU = 0.0d0
  LHSPLS = 0.0d0
  lhsgq = 0.0d0
  rhsgq = 0.0d0
!!!  call FaceAssembly_NS_conv(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!!                            acgmAlpha, pgAlpha, phigAlpha)

!!!  call FaceAssembly_NS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!!                       acgmAlpha, pgAlpha, phigAlpha)
  if (myid .eq. 0) then
    call CPU_TIME(t1)
  end if
  call IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                    acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, 0)

  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, 0)
  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if
  ! Compute residuals
  if (numnodes > 1) then
    call commu(RHSGp, 1, 'in ')
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGls, 1, 'in ')
    call commu(rhsgq, 1, 'in ')
    call commu(rhsgq, 1, 'out ')
    call commu(lhsgq, 1, 'in ')
    call commu(lhsgq, 1, 'out ')
  end if

  MomResl = sum(RHSGu(:, 1)**2) + sum(RHSGu(:, 2)**2) + sum(RHSGu(:, 3)**2)
  ConResl = sum(RHSGp**2)
  LSResl = sum(RHSGls**2)

  MomRes = MomResl
  ConRes = ConResl
  LSRes = LSResl

  if (numnodes > 1) then
    call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    call MPI_ALLREDUCE(ConResl, ConRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
    call MPI_ALLREDUCE(LSResl, LSRes, 1, MPI_DOUBLE_PRECISION, &
                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
  end if

  MomRes = sqrt(MomRes)
  ConRes = sqrt(ConRes)
  LSRes = sqrt(LSRes)

  ! Print residual
  if (inewt == 1) then
    MomRes0 = MomRes
    ConRes0 = ConRes
    LSRes0 = LSRes
  end if

  if (ismaster) then
    fomt = "(I3,a,x,ES13.6,x,F12.6)"
    write (*, fomt) inewt, ") Total Mom. Res. Norm = ", &
      MomRes, 1.0d2*MomRes/MomRes0

    write (*, fomt) inewt, ") Continuity Res. Norm = ", &
      ConRes, 1.0d2*ConRes/ConRes0

    write (*, fomt) inewt, ") LS Res. Norm = ", &
      LSRes, 1.0d2*LSRes/LSRes0
    write (*, *)
  end if

  if (isnan(MomRes + ConRes + LSRes)) then
    write (*, *) myid, inewt, '!=== Norm is NaN ===!'
    stop
  end if

  ! Check convergence and solve linear system
  if (((MomRes/MomRes0) < NS_NL_Utol) .and. &
      ((ConRes/ConRes0) < NS_NL_Ptol) .and. &
      ((LSRes/LSRes0) < LSC_NL_tol)) then
    NSConverged = .true.
    if (ismaster) then
      write (*, *)
      write (*, *) "    Converged: skip linear solve"
      write (*, *)
    end if

  else
    NSConverged = .false.
    acgAlpha = 0.0d0; pgAlpha = 0.0d0; rphigAlpha = 0.0d0
!    lhsLSu = 0.0d0!; lhsUls = 0.0d0; lhsPls = 0.0d0
    if (myid .eq. 0) then
      call CPU_TIME(t1)
    end if
    call SparseGMRES_up(col, row, IBC, IPER, D_FLAG, P_FLAG, &
                        rhsGu, rhsGp, acgAlpha, pgAlpha, &
                        lhsK11, lhsG, lhsD1, lhsM, icnt, &
                        NS_GMRES_tol, NS_GMRES_itermax, &
                        NS_GMRES_itermin, &
                        NNODE, maxNSHL, NSD, &
                        rphigAlpha, lhsLS, &
                        lhsLSu, lhsUls, lhsPls, rhsGls)

    if (myid .eq. 0) then
      call CPU_TIME(t2)
      write (*, *) "Total time solve:", t2 - t1, "seconds"
    end if

    acg = acg + acgAlpha
    ug = ug + gami*Delt*acgAlpha
    pg = pg + alfi*gami*Delt*pgAlpha

    rphig = rphig + rphigAlpha
    phig = phig + gami*Delt*rphigAlpha
  end if

  !--------------------------------------------------------------------
  ! compute torque (weak BC/Surface Int.)
  !--------------------------------------------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg
  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + alfi*(rphig - rphigold)

!!$  ! compute RHS and the torque for individual blade
!!$  call compute_RHS_PhVI(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!$                        acgmAlpha, pgAlpha, phigAlpha)

  ! compute RHS and the torque (RHSGu is used afterwards for FSI)
!  call compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!                   acgmAlpha, pgAlpha, phigAlpha, 4,     &  ! 1-nrel; 4-5mw
!                   torque3, torque4, force_trac, force_cons)

!!$!!$  ! output the force to compare with shell
!!$!!$  RHSGm = 0.0d0
!!$!!$  do b = 1, NBOUND
!!$!!$    if (bound(b)%Face_ID == 4) then
!!$!!$      do i = 1, bound(b)%NNODE
!!$!!$        j = bound(b)%BNODES(i)
!!$!!$        RHSGm(j,:) = RHSGu(j,:)
!!$!!$      end do
!!$!!$    end if
!!$!!$  end do
!!$!!$
!!$!!$  MomResl = sum(RHSGm(:,1)**2) + sum(RHSGm(:,2)**2) + sum(RHSGm(:,3)**2)
!!$!!$  MomRes  = MomResl
!!$!!$
!!$!!$  if (numnodes > 1) then
!!$!!$    call MPI_ALLREDUCE(MomResl, MomRes, 1, MPI_DOUBLE_PRECISION, &
!!$!!$                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
!!$!!$  end if
!!$!!$
!!$!!$  MomRes = sqrt(MomRes)
!!$!!$
!!$!!$  if (ismaster) then
!!$!!$    write(*,*) "*) Total Force Norm = ", MomRes
!!$!!$  end if
!!$
!!$  !-- compute conservative force and torque (Strong BC/Volume Int.) ---
!!$  if (BCugType(4,1) == 3 .or. BCugType(5,1) == 3 .or. &
!!$      BCugType(6,1) == 3 .or. BCugType(7,1) == 3) then
!!$
!!$    RHSGu = 0.0d0
!!$    call IntElmAss_3D(inewt, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
!!$                      acgmAlpha, pgAlpha, phigAlpha, 1)
!!$
!!$    if (numnodes > 1) call commu(RHSGu, NSD, 'in ')
!!$
!!$    ! only leave vaules on the boundary
!!$    RHSGm = 0.0d0
!!$    do b = 1, NBOUND
!!$      if (bound(b)%Face_ID == 21 .or. bound(b)%Face_ID == 22 .or. &
!!$          bound(b)%Face_ID == 23 .or. bound(b)%Face_ID == 24) then
!!$        do i = 1, bound(b)%NNODE
!!$          j = bound(b)%BNODES(i)
!!$          RHSGm(j,:) = RHSGu(j,:)
!!$        end do
!!$      end if
!!$    end do
!!$
!!$    call compute_moment3(NNODE, NSD, RHSGm, xg, dgAlpha, torque2)
!!$  end if
!!$  !-- end -----------------------------------------------------------

  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

end subroutine solveNavSto

!======================================================================
! Note: RHSGu is used in solving the Shell problem. Keep it!!
!======================================================================
subroutine compute_RHS(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                       acgmAlpha, pgAlpha, phigAlpha, flag, &
                       tq_trac, tq_cons, f_trac, f_cons)
  use aAdjKeep
  use mpi
  use commonvars
  implicit none

  integer, intent(in) :: flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), &
                         ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
                         pgAlpha(NNODE), phigAlpha(NNODE)
  real(8), intent(out) :: tq_trac, tq_cons, f_trac(NSD), f_cons(NSD)
  integer :: i

  ! usd RHSGu and RHSGm to store integrated traction vector
  RHSGu = 0.0d0; RHSGm = 0.0d0

  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, flag)

  if (numnodes > 1) then
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGm, NSD, 'in ')
  end if

!  do i = 1, NNODE
!  if(abs((RHSGu(i,1))-(RHSGm(i,1))) > 1.0d-6) then
!    write(*,*) 'X:', myid, i, (RHSGu(i,1)), (RHSGm(i,1))
!  end if
!  if(abs((RHSGu(i,2))-(RHSGm(i,2))) > 1.0d-6) then
!    write(*,*) 'Y:', myid, i, (RHSGu(i,2)), (RHSGm(i,2))
!  end if
!  if(abs((RHSGu(i,3))-(RHSGm(i,3))) > 1.0d-6) then
!    write(*,*) 'Z:', myid, i, (RHSGu(i,3)), (RHSGm(i,3))
!  end if
!  end do
  ! pressure and viscous force
!  call compute_moment3(NNODE, NSD, RHSGm, xg, dgAlpha, tq_trac)

  ! conservative torque for weak BC
!  call compute_moment3(NNODE, NSD, RHSGu, xg, dgAlpha, tq_cons)

  !===================================================================
  !Compute x,y,z components of the forces on the sphere
  !X:
  ! pressure and viscous force
  call compute_force1(NNODE, NSD, RHSGm, xg, dgAlpha, f_trac(1))

  ! conservative torque for weak BC
  call compute_force1(NNODE, NSD, RHSGu, xg, dgAlpha, f_cons(1))
  !Y:
  ! pressure and viscous force
  call compute_force2(NNODE, NSD, RHSGm, xg, dgAlpha, f_trac(2))

  ! conservative torque for weak BC
  call compute_force2(NNODE, NSD, RHSGu, xg, dgAlpha, f_cons(2))
  !Z:
  ! pressure and viscous force
  call compute_force3(NNODE, NSD, RHSGm, xg, dgAlpha, f_trac(3))

  ! conservative torque for weak BC
  call compute_force3(NNODE, NSD, RHSGu, xg, dgAlpha, f_cons(3))

end subroutine compute_RHS

!======================================================================
!
!======================================================================
subroutine compute_RHS_PhVI(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha)
  use aAdjKeep
  use mpi
  use commonvars
  implicit none

  real(8), intent(in) :: dgAlpha(NNODE, NSD), &
                         ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
                         pgAlpha(NNODE), phigAlpha(NNODE)
  integer :: i
  real(8) :: vtmp(NSD), TRACT(NNODE, NSD), zero(NNODE, NSD), &
             root(NNODE, NSD)

  !------------------------------------------------------------------
  ! Compute the other two cross-product components for a single blade
  RHSGu = 0.0d0
  RHSGm = 0.0d0
  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, 2)
  if (numnodes > 1) then
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGm, NSD, 'in ')
  end if

  call compute_moment3(NNODE, NSD, RHSGu, xg, dgAlpha, torque1)

  ! rotate the force back to reference config.
  TRACT = 0.0d0
  do i = 1, NNODE
    vtmp = 0.0d0
    call rot_vec_z(NSD, theta, RHSGm(i, :), vtmp)
    TRACT(i, :) = vtmp(:)
  end do

  zero = 0.0d0
  root = 0.0d0; root(:, 1) = 0.508d0  ! centre of the root

  call compute_moment3(NNODE, NSD, TRACT, xg, zero, torque3)
!!$  call compute_moment1(NNODE, NSD, TRACT, xg , acgAlpha, moment1) ! reference
  call compute_moment1(NNODE, NSD, TRACT, xg, root, moment1) ! root
!!$  if (ismaster) write(*,*)
!!$  call compute_moment2(NNODE, NSD, TRACT, xg, acgmAlpha, moment2) ! root

  !------------------------------------------------------------------
  ! Compute the other two cross-product components for a single blade
  RHSGu = 0.0d0
  RHSGm = 0.0d0
  call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                            acgmAlpha, pgAlpha, phigAlpha, 3)
  if (numnodes > 1) then
    call commu(RHSGu, NSD, 'in ')
    call commu(RHSGm, NSD, 'in ')
  end if

  call compute_moment3(NNODE, NSD, RHSGu, xg, dgAlpha, torque2)

  ! rotate the force back to reference config.
  TRACT = 0.0d0
  do i = 1, NNODE
    vtmp = 0.0d0
    call rot_vec_z(NSD, theta, RHSGm(i, :), vtmp)
    TRACT(i, :) = vtmp(:)
  end do

  zero = 0.0d0
  root = 0.0d0; root(:, 1) = -0.508d0   ! centre of the root

  call compute_moment3(NNODE, NSD, TRACT, xg, zero, torque4)
  call compute_moment1(NNODE, NSD, TRACT, xg, root, moment2) ! root

  !-- end --------------------

  if (ismaster) then
    write (*, *)
    write (*, *) "    torque  =", torque1 + torque2, torque3 + torque4
  end if
end subroutine compute_RHS_PhVI

!======================================================================
!
!======================================================================
subroutine compute_force1(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = force(i, 1)

  end do
  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    force-x =", moment
  end if
end subroutine compute_force1

!======================================================================
!
!======================================================================
subroutine compute_force2(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = force(i, 2)

  end do
  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    force-y =", moment
  end if
end subroutine compute_force2

!======================================================================
!
!======================================================================
subroutine compute_force3(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = force(i, 3)

  end do
  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    force-z =", moment
  end if
end subroutine compute_force3

!======================================================================
!
!======================================================================
subroutine compute_moment3(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl, Center(3)

  Center = 0.0d0
  ! Center(3) = 90.0d0

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = (xg(i, 2) - Center(2) + dg(i, 2))*force(i, 3) - &
               (xg(i, 3) - Center(3) + dg(i, 3))*force(i, 2)

  end do

  momentl = sum(force(:, 1))
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

!  if (ismaster) then
!    write(*,*) "    force X =", moment
!  end if

  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    moment3 =", moment
  end if
end subroutine compute_moment3

!======================================================================
!
!======================================================================
subroutine compute_moment2(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = (xg(i, 3) + dg(i, 3))*force(i, 1) - &
               (xg(i, 1) + dg(i, 1))*force(i, 3)

  end do
  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    moment2 =", moment
  end if
end subroutine compute_moment2

!======================================================================
!
!======================================================================
subroutine compute_moment1(NNODE, NSD, force, xg, dg, moment)
  use mpi

  implicit none

  integer, intent(in)  :: NNODE, NSD
  real(8), intent(in)  :: xg(NNODE, NSD), dg(NNODE, NSD), force(NNODE, NSD)
  real(8), intent(out) :: moment

  integer :: i
  real(8) :: cprod(NNODE), momentl

  cprod = 0.0d0
  do i = 1, NNODE
    cprod(i) = (xg(i, 2) + dg(i, 2))*force(i, 3) - &
               (xg(i, 3) + dg(i, 3))*force(i, 2)

  end do
  momentl = sum(cprod)
  moment = 0.0d0
  call MPI_ALLREDUCE(momentl, moment, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *) "    moment1 =", moment
  end if
end subroutine compute_moment1
