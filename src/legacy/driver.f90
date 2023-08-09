!======================================================================
! Main routine to call all the subroutines
!======================================================================
program NURBScode

  use aAdjKeep
  use mpi
  use commonvars
  use defs_shell
  use configuration

  implicit none

  integer :: i, j, k, ii, istep, Rstep, nn, dd, ibld, avgstepold, avgstep
  real(8) :: Rmat(3, 3), Rdot(3, 3), Rddt(3, 3), &
             RmatOld(3, 3), RdotOld(3, 3), RddtOld(3, 3), &
             utmp(3), umtmp(3), Rtmp(3, 3), Tacc, ForceTemp(3)
  real(8), allocatable :: dshalpha(:, :)
  real(8), allocatable :: NRmat(:, :, :), NRdot(:, :, :), NRddt(:, :, :), &
                          NRmatOld(:, :, :), NRdotOld(:, :, :), NRddtOld(:, :, :)

  character(len=30) :: fname, iname, cname
  type(ConfigType) :: config

  ! Initialize MPI
  call MPI_INIT(mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numnodes, mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, mpi_err)
  ismaster = myid .eq. mpi_master

!A  Tacc = 2.0d0


!!!  solshell = myid.eq.0

  ! flag for non-matching computation
  nonmatch = .false.
  if (ismaster) write (*, *) "Get run parameters"
  call getparam()
  call init_config(config)
  ! Read mesh and MPI-communication Data
  if (ismaster) write (*, *) "Read mesh and communication data"
  call input(myid + 1)
  if (numnodes > 1) call ctypes()

  ! Generate Sparse Structures
  if (ismaster) write (*, *) "Generating sparse structure"
  call genSparStruc()

  ! Allocate Matrices and Vectors
  if (ismaster) write (*, *) "Allocating matrices and vectors"
  call allocMatVec()

  ! Read in restart files
  call readStep(Rstep)
  ! Get initial condition
  if (Rstep == 0) then
    call generateIC()
    call writeSol(Rstep)
  else
    call readSol(Rstep)
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

    !--------------------------------------------
    ! Prediction
    !--------------------------------------------

    ug = ugold
    acg = (gami - 1.0d0)/gami*acgold
    pg = pgold

    phig = phigold
    rphig = (gami - 1.0d0)/gami*rphigold

    Tg = Tgold
    rTg = (gami - 1.0d0)/gami*rTgold

    !--------------------------------------------
    ! Solve Flow
    !--------------------------------------------
    call solmultiphasethermofluid_stag(istep)
    !--------------------------------------------
    ! Update Old Quantities
    !--------------------------------------------
    dgold = dg
    ugold = ug
    acgold = acg
    ugmold = ugm
    acgmold = acgm
    pgold = pg
    rphigold = rphig
    phigold = phig
    rTgold = rTg
    Tgold = Tg


    if (mod(istep, ifq) == 0) then
      call writeSol(istep)
    end if

  end do

  !--------------------------------------------
  ! Deallocate Matrices and Vectors
  !--------------------------------------------
  if (ismaster) write (*, *) "Deallocating matrices and vectors"
  call deallocMatVec()
  call finalize_config(config)
  ! Finalize MPI
  call MPI_FINALIZE(mpi_err)

end program NURBScode