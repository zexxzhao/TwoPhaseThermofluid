!======================================================================
! Solving the Shell problem using fluid force
!======================================================================
subroutine solveKLShell(SH, inewt, istep)

  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  use defs_shell

  implicit none

  type(shell_bld), intent(inout) :: SH
  integer, intent(in)    :: inewt, istep

  real(8), allocatable :: r2tmp1(:, :), r2tmp2(:, :)
!  real(8) :: dshAlpha(SH%TSP%NNODE,NSD), ushAlpha(SH%TSP%NNODE,NSD), &
  !            ashAlpha(SH%TSP%NNODE,NSD), ltq
  real(8) :: dshAlpha(SH%NRB%NNODE, NSD), ushAlpha(SH%NRB%NNODE, NSD), &
             ashAlpha(SH%NRB%NNODE, NSD), ltq
  integer :: i, j, k, bb, inewt_SH, ibld
  character(len=80) :: fname, iname(2)

  if (ismaster) then
    write (*, *)
    write (*, '(I3,a)') inewt, ") Solving Shell Problem ---"
    write (*, *)
  end if

  !===============================================================
  ! Begin SHELL analysis
  !===============================================================

!!$  ! assign the local fluid force from volume element boundary to
!!$  ! shell mesh
!!$  allocate(r2tmp1(SH%FEM%NNODE,NSD))!
!!$
!!$  r2tmp1 = 0.0d0
!!$
!!$  bb = 7
!!$  do i = 1, bound(bb)%NNODE
!!$    r2tmp1(bound(bb)%L2SNODE(i),:) = RHSGu(bound(bb)%BNODES(i),:)
!!$  end do
!!$
!!$  SH%FEM%FORCE = r2tmp1
!!$
!!$  ! sum up the local force such that the shell force vector is complete
!!$  if (numnodes > 1) then
!!$    call MPI_ALLREDUCE(r2tmp1, SH%FEM%FORCE, SH%FEM%NNODE*NSD, &
!!$                       MPI_DOUBLE_PRECISION, &
!!$                       MPI_SUM, MPI_COMM_WORLD, mpi_err)
!!$  end if
!!$
!!$  deallocate(r2tmp1)

  ! assign the local fluid force from volume element boundary to
  ! shell mesh
  do ibld = 1, NBlade
    allocate (r2tmp1(blade(ibld)%NNODE, NSD))
    r2tmp1 = 0.0d0

    allocate (blade(ibld)%FORCE(blade(ibld)%NNODE, NSD))
    blade(ibld)%FORCE = 0.0d0

    bb = ibld + SH%bmap   ! mapping between blade and boundary...
    do i = 1, bound(bb)%NNODE
      r2tmp1(bound(bb)%L2SNODE(i), :) = RHSGu(bound(bb)%BNODES(i), :)
    end do

    blade(ibld)%FORCE = r2tmp1

    ! sum up the local force such that the shell force vector is complete
    if (numnodes > 1) then
      call MPI_ALLREDUCE(r2tmp1, blade(ibld)%FORCE, &
                         blade(ibld)%NNODE*NSD, &
                         MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, mpi_err)
    end if

    deallocate (r2tmp1)
  end do

!!$  ! output pre-integrated force for static analysis
!!$  if (mod(istep,ifreq_sh)==0 .and. ismaster .and. inewt==1) then
!!$    do ibld = 1, NBlade
!!$      write(iname(1),'(I30)') istep
!!$      write(iname(2),'(I30)') ibld + 6
!!$      fname = 'sh.force.fem.'//trim(adjustl(iname(1)))//'.'&
!!$                             //trim(adjustl(iname(2)))
!!$      open(79, file=fname, status='replace')
!!$      write(79,*) blade(ibld)%NNODE
!!$      do i = 1, blade(ibld)%NNODE
!!$        write(79,*) (blade(ibld)%FORCE(i,j), j = 1, 3)
!!$      end do
!!$    end do
!!$  end if

  !-start the shell solve -------------------------------------------

  if (solshell) then

    !-----------------------------------------------------
    ! Since fluid force is pre-integrated RHS, we need to
    ! Project it to nodal points
    !-----------------------------------------------------
    SH%FEM%FORCE = blade(myid + 1)%FORCE

    write (*, *) "Fluid x-Force", sum(SH%FEM%FORCE(:, 1))

    allocate (r2tmp1(SH%FEM%NNODE, NSD))
    r2tmp1 = 0.0d0
    call f2s_l2project(SH%FEM, NSD, SH%FEM%FORCE, r2tmp1)

!!$  ! for checking
!!$  SH%RHSGNorm = sqrt(sum(SH%FEM%FORCE(:,1)**2 + SH%FEM%FORCE(:,2)**2 + &
!!$                         SH%FEM%FORCE(:,3)**2))
!!$
!!$  if (ismaster) then
!!$    write(*,*) "*Shell Force Norm (FEM):", SH%RHSGNorm
!!$  end if
!!$
!!$  SH%RHSGNorm = sum(SH%FEM%FORCE(:,3))
!!$
!!$  if (ismaster) then
!!$    write(*,*) "*Shell Force in Z:", SH%RHSGNorm
!!$  end if

    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !--------------------------------------------------------------
    ! Get the external force (pre-integrated rhs) for t-spline
    ! Only TSP%FORCE is useful. SH%LHSK, SH%mg are useless here
    !--------------------------------------------------------------
!!!    allocate(r2tmp2(SH%TSP%NNODE,NSD))
    allocate (r2tmp2(SH%NRB%NNODE, NSD))
    r2tmp2 = 0.0d0
    !SH%TSP%FORCE = 0.0d0
    SH%NRB%FORCE = 0.0d0

!!!    call f2t_IntElmAss(SH%TSP, SH%BEZ, SH%FEM, SH%icnt, SH%col, SH%row, &
!!!                       nsd, r2tmp1, SH%TSP%FORCE, SH%LHSK, r2tmp2)

    call f2n_IntElmAss(SH%NRB, SH%FEM, SH%icnt, SH%col, SH%row, &
                       nsd, r2tmp1, SH%NRB%FORCE, SH%LHSK, r2tmp2)

    SH%LHSK = 0.0d0
    deallocate (r2tmp1, r2tmp2)

    ! converge geometrical nonlinearity
    do inewt_SH = 1, SH%Nnewt(inewt)

      ! Get quantities at alpha levels:
!!!      ashAlpha = SH%TSP%ashOld + almi*(SH%TSP%ash-SH%TSP%ashOld)
!!!      ushAlpha = SH%TSP%ushOld + alfi*(SH%TSP%ush-SH%TSP%ushOld)
!!!      dshAlpha = SH%TSP%dshOld + alfi*(SH%TSP%dsh-SH%TSP%dshOld)

      ashAlpha = SH%NRB%ashOld + almi*(SH%NRB%ash - SH%NRB%ashOld)
      ushAlpha = SH%NRB%ushOld + alfi*(SH%NRB%ush - SH%NRB%ushOld)
      dshAlpha = SH%NRB%dshOld + alfi*(SH%NRB%dsh - SH%NRB%dshOld)

      ! initialization
      SH%RHSG = 0.0d0    ! initialize
      SH%RHSG_EXT = 0.0d0
      SH%RHSG_GRA = 0.0d0
      SH%LHSK = 0.0d0

!!!      SH%TSP%B_NET_D(:,1:3) = SH%TSP%B_NET(:,1:3) + dshAlpha(:,:)
!!!      call IntElmAss_tsp_sh(SH%TSP, SH%BEZ, SH%icnt, SH%col, SH%row, &
!!!                            nsd, SH%NMat, SH%matA, SH%matB, SH%matD, &
!!!                            SH%G_fact, SH%RHSG, SH%RHSG_GRA, SH%LHSK,  &
!!!                            SH%T_Flag, ashAlpha, alfi, beti, almi, &
!!!                            Delt, SH%BldRot)

      SH%NRB%B_NET_D(:, 1:3) = SH%NRB%B_NET(:, 1:3) + dshAlpha(:, :)
      call IntElmAss_shell(SH, SH%NRB, SH%icnt, SH%col, SH%row, &
                           nsd, &
                           SH%G_fact, SH%RHSG, SH%RHSG_GRA, SH%LHSK, &
                           ashAlpha, alfi, beti, almi, &
                           Delt)

      ! external force
!      SH%RHSG_EXT = SH%TSP%FORCE
      SH%RHSG_EXT = SH%NRB%FORCE

      write (*, *) "Structure x-Force", sum(SH%NRB%FORCE(:, 1))

      ! zero out the residual again
      !     do i = 1, SH%TSP%NNODE
      !      if (SH%TSP%IBC(i,1) == 1) then
      do i = 1, SH%NRB%NNODE
        if (SH%NRB%IBC(i, 1) == 1) then
          SH%RHSG_EXT(i, :) = 0.0d0
          SH%RHSG_GRA(i, :) = 0.0d0
        end if
      end do

      SH%RHSG = SH%RHSG + SH%RHSG_EXT + SH%RHSG_GRA

      !=========================================
      ! Solve Kx = f
      !=========================================
!!$      SH%RHSGNorm = sum(SH%RHSG_EXT(:,3))
!!$      if (ismaster) then
!!$        write(*,*) "*External Force in Z (TSP):", SH%RHSGNorm
!!$      end if

      SH%RHSGNorm = sqrt(sum(SH%RHSG(:, 1)**2 + &
                             SH%RHSG(:, 2)**2 + &
                             SH%RHSG(:, 3)**2))

      if (ismaster) then
        write (*, '(4X,A,I4,A,ES14.6)') "*Shell Inner Iter:", inewt_SH, &
          ",  Res. Norm =", SH%RHSGNorm
        write (*, *)
      end if

      if (SH%RHSGNorm < SH%RHSGtol) exit

      ! Increment for the geometric nonlinearity
      ashAlpha = 0.0d0
!      call SparseCG_BDIAG_shell(SH%TSP%NNODE, SH%TSP%maxNSHL, NSD, &
      !                               SH%icnt, SH%col, SH%row, &
      !                              SH%LHSK, SH%RHSG, &
      !                             ashAlpha, 1.0d-3, 1)

      call SparseCG_BDIAG_shell(SH%NRB%NNODE, SH%NRB%maxNSHL, NSD, &
                                SH%icnt, SH%col, SH%row, &
                                SH%LHSK, SH%RHSG, &
                                ashAlpha, 1.0d-3, 1)

      ! Update Current Values
!      SH%TSP%ash(:,:) = SH%TSP%ash(:,:) + ashAlpha(:,:)
      !     SH%TSP%ush(:,:) = SH%TSP%ush(:,:) + ashAlpha(:,:)*gami*Delt
      !    SH%TSP%dsh(:,:) = SH%TSP%dsh(:,:) + ashAlpha(:,:)*beti*Delt*Delt

      SH%NRB%ash(:, :) = SH%NRB%ash(:, :) + ashAlpha(:, :)
      SH%NRB%ush(:, :) = SH%NRB%ush(:, :) + ashAlpha(:, :)*gami*Delt
      SH%NRB%dsh(:, :) = SH%NRB%dsh(:, :) + ashAlpha(:, :)*beti*Delt*Delt

      ! project the increment back to NURBS
      allocate (r2tmp1(SH%FEM%NNODE, NSD))
      r2tmp1 = 0.0d0
!      call t2f_l2project(SH%FEM, SH%TSP, SH%BEZ, NSD, ashAlpha, r2tmp1)
      call n2f_l2project(SH%FEM, SH%NRB, NSD, ashAlpha, r2tmp1)
      ! zero out the solution at constrained nodes
      do i = 1, SH%FEM%NNODE
        if (SH%FEM%IBC(i, 1) == 1) then
          r2tmp1(i, :) = 0.0d0
        end if
      end do
      SH%FEM%ash(:, :) = SH%FEM%ash(:, :) + r2tmp1(:, :)
      SH%FEM%ush(:, :) = SH%FEM%ush(:, :) + r2tmp1(:, :)*gami*Delt
      SH%FEM%dsh(:, :) = SH%FEM%dsh(:, :) + r2tmp1(:, :)*beti*Delt*Delt

!!!    if (ismaster) then
!!!!      write(*,'(1X,A,2F13.9)') "Tip disp. of FEM and TSP:", &
!!!      write(*,*) "*Tip displacement and increment (TSP|FEM):"
!!!      write(*,*)  SH%TSP%dsh(maxloc(SH%TSP%B_NET(:,2),dim=1),3), &
!!!                  ashAlpha(maxloc(SH%TSP%B_NET(:,2),dim=1),3)*beti*Delt*Delt
!!!      write(*,*)  SH%FEM%dsh(maxloc(SH%FEM%B_NET(:,2),dim=1),3), &
!!!                  r2tmp1(maxloc(SH%FEM%B_NET(:,2),dim=1),3)*beti*Delt*Delt
!!!    end if

      deallocate (r2tmp1)

    end do ! end of local Newton iteration (Nnew_SH)

    ! aerodynamic torque
    SH%Tq1 = sum((SH%FEM%B_NET(:, 2) - Center_Rot(2) + SH%FEM%dsh(:, 2))*SH%FEM%FORCE(:, 3) - &
                 (SH%FEM%B_NET(:, 3) - Center_Rot(3) + SH%FEM%dsh(:, 3))*SH%FEM%FORCE(:, 2))

!    SH%Tq2 = sum((SH%TSP%B_NET(:,1)+SH%TSP%dsh(:,1))*SH%TSP%FORCE(:,2)- &
    !                (SH%TSP%B_NET(:,2)+SH%TSP%dsh(:,2))*SH%TSP%FORCE(:,1))

    SH%Tq2 = sum((SH%NRB%B_NET(:, 2) - Center_Rot(2) + SH%NRB%dsh(:, 2))*SH%NRB%FORCE(:, 3) - &
                 (SH%NRB%B_NET(:, 3) - Center_Rot(3) + SH%NRB%dsh(:, 3))*SH%NRB%FORCE(:, 2))

!!!    if (ismaster) then
!      write(*,'(5X,A,I2,A)'   ) "Blade", myid+1, &
!                                " Leading Edge Deflection (NRB|FEM) = "
    !     write(*,'(5X,3ES15.6)'  ) SH%NRB%dsh(SH%NRB%TipLoc,:)
    !    write(*,'(5X,3ES15.6)'  ) SH%FEM%dsh(SH%FEM%TipLoc,:)

    !   write(*,'(5X,A,I2,A)'   ) "Blade", myid+1, &
    !                               " Trailing Edge Deflection (NRB|FEM) = "
    !  write(*,'(5X,3ES15.6)'  ) SH%NRB%dsh(SH%NRB%TipLocTr,:)
    ! write(*,'(5X,3ES15.6)'  ) SH%FEM%dsh(SH%FEM%TipLocTr,:)
!!!      write(*,'(5X,A,2ES15.6)') "Aerodynamic Torque = ", SH%Tq1, SH%Tq2
!!!    end if

  end if

  !-end solshell -----------------------------------------------------

!!$  ! map back the displacement of the blade
!!$  bb = 7
!!$  ! for weak BC, do NOT set ug = ugm
!!$  if (BCugType(bound(bb)%FACE_ID,1) == 2 .or. &
!!$      BCugType(bound(bb)%FACE_ID,2) == 2 .or. &
!!$      BCugType(bound(bb)%FACE_ID,3) == 2) then
!!$    do i = 1, bound(bb)%NNODE
!!$      acgm(bound(bb)%BNODES(i),:) = SH%FEM%ash(bound(bb)%L2SNODE(i),:)
!!$      ugm( bound(bb)%BNODES(i),:) = SH%FEM%ush(bound(bb)%L2SNODE(i),:)
!!$      dg(  bound(bb)%BNODES(i),:) = SH%FEM%dsh(bound(bb)%L2SNODE(i),:)
!!$    end do
!!$  ! for strong BC, ug = ugm
!!$  else
!!$    do i = 1, bound(bb)%NNODE
!!$      acgm(bound(bb)%BNODES(i),:) = SH%FEM%ash(bound(bb)%L2SNODE(i),:)
!!$      ugm( bound(bb)%BNODES(i),:) = SH%FEM%ush(bound(bb)%L2SNODE(i),:)
!!$      dg(  bound(bb)%BNODES(i),:) = SH%FEM%dsh(bound(bb)%L2SNODE(i),:)
!!$      acg( bound(bb)%BNODES(i),:) = SH%FEM%ash(bound(bb)%L2SNODE(i),:)
!!$      ug(  bound(bb)%BNODES(i),:) = SH%FEM%ush(bound(bb)%L2SNODE(i),:)
!!$    end do
!!$  end if

  ! Compute total torque
  ltq = 0.0d0; torque_SH = 0.0d0
  if (solshell) ltq = SH%Tq1
  call MPI_ALLREDUCE(ltq, torque_SH, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *)
    write (*, *) "    Torque_SH =", torque_SH
    write (*, *)
  end if

  if (solshell) ltq = SH%Tq2
  call MPI_ALLREDUCE(ltq, torque_SH, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, mpi_err)

  if (ismaster) then
    write (*, *)
    write (*, *) "    Torque_SH =", torque_SH
    write (*, *)
  end if

  ! all processors allocate the boundary solution array
  do ibld = 1, NBlade
    allocate (blade(ibld)%ash(blade(ibld)%NNODE, NSD))
    allocate (blade(ibld)%ush(blade(ibld)%NNODE, NSD))
    allocate (blade(ibld)%dsh(blade(ibld)%NNODE, NSD))
    blade(ibld)%ash = 0.0d0
    blade(ibld)%ush = 0.0d0
    blade(ibld)%dsh = 0.0d0
  end do

  ! only those who solve the shell problem would assign
  ! the solution
  if (solshell) then
!    if (myid == 0) then
    blade(1)%ash = SH%FEM%ash
    blade(1)%ush = SH%FEM%ush
    blade(1)%dsh = SH%FEM%dsh
!    else if (myid == 1) then
!      blade(2)%ash = SH%FEM%ash
!      blade(2)%ush = SH%FEM%ush
!      blade(2)%dsh = SH%FEM%dsh
!    else if (myid == 2) then
!      blade(3)%ash = SH%FEM%ash
!      blade(3)%ush = SH%FEM%ush
!      blade(3)%dsh = SH%FEM%dsh
!    end if
  end if

  ! broadcast such that all processors will have the solutions
  do ibld = 1, NBlade
    call MPI_BCAST(blade(ibld)%ash, blade(ibld)%NNODE*NSD, &
                   MPI_DOUBLE_PRECISION, ibld - 1, &
                   MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(blade(ibld)%ush, blade(ibld)%NNODE*NSD, &
                   MPI_DOUBLE_PRECISION, ibld - 1, &
                   MPI_COMM_WORLD, mpi_err)

    call MPI_BCAST(blade(ibld)%dsh, blade(ibld)%NNODE*NSD, &
                   MPI_DOUBLE_PRECISION, ibld - 1, &
                   MPI_COMM_WORLD, mpi_err)
  end do

  ! Now for each blade map back the solution to volume boundary
  do ibld = 1, NBlade
    bb = ibld + SH%bmap   ! mapping between blade and boundary...

    ! for weak BC, do NOT set ug = ugm
    if (BCugType(bound(bb)%FACE_ID, 1) == 2 .or. &
        BCugType(bound(bb)%FACE_ID, 2) == 2 .or. &
        BCugType(bound(bb)%FACE_ID, 3) == 2) then
      do i = 1, bound(bb)%NNODE
        acgm(bound(bb)%BNODES(i), :) = blade(ibld)%ash(bound(bb)%L2SNODE(i), :)
        ugm(bound(bb)%BNODES(i), :) = blade(ibld)%ush(bound(bb)%L2SNODE(i), :)
        dg(bound(bb)%BNODES(i), :) = blade(ibld)%dsh(bound(bb)%L2SNODE(i), :)
      end do
      ! for strong BC, ug = ugm
    else
      do i = 1, bound(bb)%NNODE
        acgm(bound(bb)%BNODES(i), :) = blade(ibld)%ash(bound(bb)%L2SNODE(i), :)
        ugm(bound(bb)%BNODES(i), :) = blade(ibld)%ush(bound(bb)%L2SNODE(i), :)
        dg(bound(bb)%BNODES(i), :) = blade(ibld)%dsh(bound(bb)%L2SNODE(i), :)
        acg(bound(bb)%BNODES(i), :) = blade(ibld)%ash(bound(bb)%L2SNODE(i), :)
        ug(bound(bb)%BNODES(i), :) = blade(ibld)%ush(bound(bb)%L2SNODE(i), :)
      end do
    end if
  end do

  do ibld = 1, NBlade
    deallocate (blade(ibld)%ash, blade(ibld)%ush, blade(ibld)%dsh, blade(ibld)%FORCE)
  end do

  if (numnodes > 1) call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

!  if (ismaster) then
!    write(*,'(I3,a)') inewt,") End Shell Problem ---"
!  end if

end subroutine solveKLShell
