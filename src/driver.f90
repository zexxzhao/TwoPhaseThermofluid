!======================================================================
! Main routine to call all the subroutines
!======================================================================
program NURBScode

  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  implicit none

  type(shell_bld) :: SH       ! shell for the blade
  type(shell_nmb) :: NM       ! shell for non-matching (zones)

  integer :: i, j, k, ii, istep, Rstep, nn, dd, ibld, avgstepold, avgstep
  real(8) :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3), &
             RmatOld(3, 3), RdotOld(3, 3), RddtOld(3, 3), &
             utmp(3), umtmp(3), Rtmp(3, 3), Tacc, ForceTemp(3)
  real(8), allocatable :: dshalpha(:, :)
  real(8), allocatable :: NRmat(:, :, :), NRdot(:, :, :), NRddt(:, :, :), &
                          NRmatOld(:, :, :), NRdotOld(:, :, :), NRddtOld(:, :, :)

  character(len=30) :: fname, iname, cname

  ! Initialize MPI
  call MPI_INIT(mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_err)
  ismaster = myid .eq. mpi_master

!A  Tacc = 2.0d0

  Center_Rot(1) = 7.5d0
  Center_Rot(2) = 0.0d0
  Center_Rot(3) = 90.0d0

!  solshell = (myid.eq.0)!.or.(myid.eq.1).or.(myid.eq.2)
  solshell = .false.
  NBlade = 0
!!!  solshell = myid.eq.0

  ! flag for non-matching computation
  nonmatch = .false.
  if (ismaster) write (*, *) "Get run parameters"
  call getparam()
  ! Read mesh and MPI-communication Data
  if (ismaster) write (*, *) "Read mesh and communication data"
  call input(myid + 1)
  if (numnodes > 1) call ctypes()

  !------------------------------------------------------------
  ! Only the first three processors are used to solve the blade
  ! problem (there are three blades...)
  if (solshell) then
    ! Read shell mesh
    if (ismaster) write (*, *) "Read shell mesh"
    call input_shell_blade(NSD, SH)
  end if
  ! For the first three processors, get the number of
  ! fem shell nodes of the blades
  blade(:)%NNODE = 0
  if (solshell) then
    blade(myid + 1)%NNODE = SH%FEM%NNODE
  end if

  ! Broadcast them to all processors for the purpose of
  ! allocating arrays
  do ibld = 1, NBlade
    call MPI_BCAST(blade(ibld)%NNODE, 1, MPI_INTEGER, &
                   ibld - 1, MPI_COMM_WORLD, mpi_err)
  end do

  !----------------------------------------------------
  ! find how many surfaces numbered before the object
  !----------------------------------------------------
  SH%bmap = 0
  do i = 1, NBOUND
    if (bound(i)%Face_ID <= 20) SH%bmap = SH%bmap + 1
  end do
  !-------------------------------------------------------------
  ! non-matching
  !-------------------------------------------------------------
  if (nonmatch) then
    call input_shell_nmb(NSD, NM)
  end if
  ! Get run parameters


  ! Get element gauss points (this needs to be after "getparam")
!  allocate(ELMNGAUSS(NELEM))
!  do i = 1, NELEM
!    if (ELMNSHL(i) == 4) then
!      ELMNGAUSS(i) = NGAUSSTET
!    else if (ELMNSHL(i) == 6) then
!      ELMNGAUSS(i) = NGAUSSPRI
!    else
!      write(*,*) "ERROR: Undefined ELMNSHL for ELMNGAUSS", ELMNSHL
!    end if
!  end do

  ! parameter for shell
  allocate (SH%Nnewt(NS_NL_itermax))
  SH%Nnewt = 1

  ! Generate Sparse Structures
  if (ismaster) write (*, *) "Generating sparse structure"
  call genSparStruc()

  ! Allocate Matrices and Vectors
  if (ismaster) write (*, *) "Allocating matrices and vectors"
  call allocMatVec(SH, NM)

!!!!!  call MPI_BARRIER (MPI_COMM_WORLD, mpi_err)
!!!!!  write(*,*) 'here', myid

!!$  ! Compute mass matrix and mesh data
!!$  if (ismaster) write(*,*) "Compute mass matrix and obtain meshsize"
!!$  call IntElmAss_mass
!!$  call FaceAssembly_area

  ! Read in restart files
  call readStep(Rstep)
  ! Get initial condition
  if (Rstep == 0) then
    ! x_inflow = 0.0d0
!    call readinlet()
    call generateIC(SH)
    call writeSol(Rstep)
    ! x_inflow = 66d0
!!!    call writeRB(Rstep)

  else
    call readSol(Rstep)
!    call readAvgSol(avgstepold)

    if (nonmatch) call readShellSol_NM(Rstep, NM)

    if (solshell) then
!      SH%TSP%dshOld = 0.0d0; SH%FEM%dshOld = 0.0d0
!      SH%TSP%ushOld = 0.0d0; SH%FEM%ushOld = 0.0d0
!      SH%TSP%ashOld = 0.0d0; SH%FEM%ashOld = 0.0d0
!      SH%TSP%dsh = SH%TSP%dshOld; SH%FEM%dsh = SH%FEM%dshOld
!      SH%TSP%ush = SH%TSP%ushOld; SH%FEM%ush = SH%FEM%ushOld
!      SH%TSP%ash = SH%TSP%ashOld; SH%FEM%ash = SH%FEM%ashOld

      call readShellSol(Rstep, SH)
    end if
  end if

  !------------------------------------------
  ! Loop over time steps
  !------------------------------------------
  avgstep = 0
  do istep = Rstep + 1, Nstep

    avgstep = avgstep + 1

    time = time + Delt

    if (ismaster) then
      write (*, '(60("="))')
      write (*, "(a,x,I8,x,ES14.6)") "Time Step Number:", istep, time
      write (*, '(60("="))')
    end if

    ! Predictor: Same Velocity
    ug = ugold
    acg = (gami - 1.0d0)/gami*acgold
    pg = pgold

    phig = phigold
    rphig = (gami - 1.0d0)/gami*rphigold

    vbn1 = vbn0
    dbn1 = dbn0 + Delt*vbn0
    wbn1 = wbn0

    ! solve flow
    call solmultiphasethermofluid_stag(istep)
    ! Flags
    ! move = .false.!(time >= move_time)
    ! mono = .false.!(time >= mono_time).or.move
    ! conv = .false.!(time >= conv_time).or.mono
    ! nonmatch = .false.
    ! shel = .false.!(time >= shel_time)

    ! if (mono) then
    !   call solflow_mono(istep)
    ! else
    !   call solflow_stag(istep, SH, NM)
    ! end if

    !--------------------------------------------
    ! Output traction: rotate back to reference
    !--------------------------------------------

    if (mod(istep, ifq) == 0) then
      call writeSol(istep)
!!!      call writeRB (istep)
    end if
    !call writeVelocity(istep)
    !--------------------------------------------
    ! Update Old Quantities
    !--------------------------------------------
    ! call update_sol(SH, NM)

!    !--------------------------------------------
!    ! Get the averaged (in time) solutions
!    !--------------------------------------------
!    ! averaged relative velocity and pressure
    do i = 1, NNODE
!      utmp = 0.0d0; umtmp = 0.0d0
!      call rot_vec_z(NSD, theta,  ug(i,:),  utmp)
!      call rot_vec_z(NSD, theta, ugm(i,:), umtmp)
      uavg(i, :) = uavg(i, :) + ug(i, :)
    end do
    pavg = pavg + phig
!    ! use the same ifq as solution outputs so that when restarted
!    ! you won't average the same solution again
    if (mod(istep, ifq) == 0) then
!    call writeAvgSol(avgstepold+avgstep)
    end if

  end do

  ! Deallocate Matrices and Vectors
  if (ismaster) write (*, *) "Deallocating matrices and vectors"
  call deallocMatVec

  ! Finalize MPI
  call MPI_FINALIZE(mpi_err)

end program NURBScode

!======================================================================
! subroutine to rotate the vector back to reference
!======================================================================
subroutine rot_vec_z(nsd, theta, sol, rot)
  implicit none
  integer, intent(in)  :: nsd
  real(8), intent(in)  :: theta, sol(nsd)
  real(8), intent(out) :: rot(nsd)

  real(8) :: tmpx, tmpy

  tmpx = sol(1)
  tmpy = sol(2)
  rot(1) = cos(-theta)*tmpx - sin(-theta)*tmpy
  rot(2) = sin(-theta)*tmpx + cos(-theta)*tmpy
  rot(3) = sol(3)
end subroutine rot_vec_z

!======================================================================
! subroutine to update the solution
!======================================================================
subroutine update_sol(SH, NM)
  use aAdjKeep
  use commonvars
  use mpi
  use defs_shell
  implicit none

  type(shell_bld), intent(inout) :: SH
  type(shell_nmb), intent(inout) :: NM
  integer :: i

  ! Update Old Quantities
  dgold = dg
  ugold = ug
  acgold = acg
  ugmold = ugm
  acgmold = acgm
  pgold = pg
  rphigold = rphig
  phigold = phig

  vbn0 = vbn1
  dbn0 = dbn1
  wbn0 = wbn1
  Rn0 = Rn1

  thetaOld = theta
  thetdOld = thetd
  theddOld = thedd

  if (solshell) then
!    SH%TSP%dshOld = SH%TSP%dsh
!    SH%TSP%ushOld = SH%TSP%ush
!    SH%TSP%ashOld = SH%TSP%ash

    SH%NRB%dshOld = SH%NRB%dsh
    SH%NRB%ushOld = SH%NRB%ush
    SH%NRB%ashOld = SH%NRB%ash

    SH%FEM%dshOld = SH%FEM%dsh
    SH%FEM%ushOld = SH%FEM%ush
    SH%FEM%ashOld = SH%FEM%ash
  end if

  if (nonmatch) then
    do i = 1, 2
      NM%FEM(i)%dshOld = NM%FEM(i)%dsh
      NM%FEM(i)%ushOld = NM%FEM(i)%ush
      NM%FEM(i)%ashOld = NM%FEM(i)%ash
    end do
  end if
end subroutine update_sol

!======================================================================
! apply the exact rotation to the predictor
!======================================================================
subroutine predictor_rot_shell(SHL, Rmat, Rdot, Rddt, RmatOld, RdotOld, &
                               RddtOld, Delt, beti, gami)
  use defs_shell
  implicit none

  type(mesh), intent(inout) :: SHL
  real(8), intent(in)    :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3), &
                            RmatOld(3, 3), RdotOld(3, 3), RddtOld(3, 3), &
                            Delt, beti, gami
  integer :: nn, dd

  forall (nn=1:SHL%NNODE, dd=1:3)
    SHL%ush(nn, dd) = SHL%ushold(nn, dd) + &
                      sum((Rdot(dd, :) - RdotOld(dd, :))*SHL%B_NET(nn, 1:3))

    SHL%ash(nn, dd) = sum(Rddt(dd, :)*SHL%B_NET(nn, 1:3)) + &
                      (gami - 1.0d0)/gami*(SHL%ashold(nn, dd) - &
                                           sum(RddtOld(dd, :)*SHL%B_NET(nn, 1:3)))

    SHL%dsh(nn, dd) = SHL%dshold(nn, dd) + &
                      sum((Rmat(dd, :) - RmatOld(dd, :))*SHL%B_NET(nn, 1:3)) + &
                      Delt*(SHL%ushold(nn, dd) - sum(RdotOld(dd, :)*SHL%B_NET(nn, 1:3))) + &
                      Delt*Delt/2.0d0*((1.0d0 - 2.0d0*beti)* &
                                       (SHL%ashold(nn, dd) - sum(RddtOld(dd, :)*SHL%B_NET(nn, 1:3))) + &
                                       2.0d0*beti*(SHL%ash(nn, dd) - sum(Rddt(dd, :)*SHL%B_NET(nn, 1:3))))
  end forall
end subroutine predictor_rot_shell

