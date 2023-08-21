!======================================================================
! Main routine to call all the subroutines
!======================================================================
program NURBScode

  use aAdjKeep
  use mpi
  use commonvars
  use class_def
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
  type(MeshData) :: mesh
  type(SparsityPattern) :: sp
  type(RHSData) :: vec
  type(LHSData) :: mat
  type(FieldData) :: solution
  type(DirichletBCData) :: bc

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
  ! Read mesh and MPI-communication Data
  if (ismaster) write (*, *) "Read mesh and communication data"
  ! call input(myid + 1)
  call input(myid + 1, mesh)
  if (numnodes > 1) call ctypes()

  ! Generate Sparse Structures
  if (ismaster) write (*, *) "Generating sparse structure"
  call genSparsityPattern(&
    mesh%NNODE, mesh%maxNSHL, mesh%NELEM, mesh%ELMNSHL, mesh%IEN, &
    sp%index, sp%indptr, sp%nnz)


  ! Allocate Matrices and Vectors
  if (ismaster) write (*, *) "Allocating matrices and vectors"
  ! call allocMatVec()
  call allocField(mesh, solution)
  call allocRHS(mesh, vec)
  call allocLHS(sp, mesh, mat)

  ! Read in restart files
  call readStep(Rstep)
  ! Get initial condition
  if (Rstep == 0) then
    call generateIC()
    call writeSol(Rstep)
  else
    call readSol(Rstep)
  end if
  call init_config(config)

  call setBCs_NSVOF(bc)
  call setBCs_Tem(bc)

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

    solution%ug(:, :) = solution%ugold(:, :)
    solution%acg(:, :) = ((gami - 1.0d0)/gami)*solution%acgold(:, :) 
    pg = pgold

    solution%phig(:) = solution%phigold(:)
    solution%rphig(:) = (gami - 1.0d0)/gami*solution%rphigold(:)

    solution%Tg = solution%Tgold
    solution%rTg = (gami - 1.0d0)/gami*solution%rTgold

    !--------------------------------------------
    ! Solve Flow
    !--------------------------------------------
    call solmultiphasethermofluid_stag(istep, config, mesh, sp, bc, solution, vec, mat)
    !--------------------------------------------
    ! Update Old Quantities
    !--------------------------------------------
    solution%dgold = solution%dg
    solution%ugold = solution%ug
    solution%acgold = solution%acg
    solution%ugmold = solution%ugm
    solution%acgmold = solution%acgm
    solution%pgold = solution%pg
    solution%rphigold = solution%rphig
    solution%phigold = solution%phig
    solution%rTgold = solution%rTg
    solution%Tgold = solution%Tg


    if (mod(istep, ifq) == 0) then
      call writeSol(istep)
    end if

  end do

  !--------------------------------------------
  ! Deallocate Matrices and Vectors
  !--------------------------------------------
  if (ismaster) write (*, *) "Deallocating matrices and vectors"
  ! call deallocMatVec()
  call finalize_config(config)
  call freeField(solution)
  call freeRHS(vec)
  call freeLHS(mat)
  ! Finalize MPI
  call MPI_FINALIZE(mpi_err)

end program NURBScode
