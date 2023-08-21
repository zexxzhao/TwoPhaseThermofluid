!======================================================================
!
!======================================================================
subroutine input(id, mesh)
  ! use aAdjKeep
  ! use commonvars
  use mpi
  use class_def
  implicit none

  integer, intent(in) :: id
  type(MeshData), intent(out) :: mesh

  character(len=30) :: fname, iname
  character(len=10) :: cname

  integer :: mfid, i, j, k, itmp1, itmp2, itmp3, nshl, counter
  integer :: NSD, NNODE, NSHLmax, NELEM, NBOUND, NPATCH

  counter = 1
  ! open mesh files
  mfid = 11
  fname = 'part'//trim(cname(id))//'.dat'
  open (mfid, file=fname, status='old')

  read (mfid, *) mesh%NSD, NSHL, mesh%NSHLBmax
  read (mfid, *) mesh%NNODE, mesh%NELEM, mesh%NBOUND, NPATCH
  NSD = mesh%NSD
  NSHLmax = mesh%NSHLBmax
  NNODE = mesh%NNODE
  NELEM = mesh%NELEM
  NBOUND = mesh%NBOUND


!  if (myid + 1 == 11) then
  !       write(*,*) counter, NSD, NSHL, NSHLBmax
  !      counter = counter + 1
  !     write(*,*) counter, NNODE, NELEM, NBOUND, NPATCH
  !    counter = counter + 1
!endif

  ! read nodes
  allocate (mesh%xg(NNODE, NSD), mesh%NodeID(NNODE))
  do i = 1, NNODE
    read (mfid, *) (mesh%xg(i, j), j=1, NSD), mesh%NodeID(i)
    ! if (myid + 1 == 1) then
    !   write (*, *) counter, xg(i, :)
    !   counter = counter + 1
    ! end if
  end do
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

  ! read elements
  allocate (mesh%ELMNSHL(NELEM), mesh%IEN(NELEM, NSHL))
  allocate (mesh%ELM_ID(NELEM))
  do i = 1, NELEM
    read (mfid, *) mesh%ELMNSHL(i), (mesh%IEN(i, j), j=1, mesh%ELMNSHL(i)), mesh%ELM_ID(i)
    ! if ((id == 120) .and. (i == NELEM)) then
    !   write (*, *), ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i))
    ! end if
  end do
  ! maxNSHL = maxval(mesh%ELMNSHL)
  ! read faces
  allocate (mesh%bound(NBOUND))
  do i = 1, NBOUND

    read (mfid, *) mesh%bound(i)%FACE_ID, mesh%bound(i)%NFACE, mesh%bound(i)%NNODE

    allocate (mesh%bound(i)%FACE_IEN(mesh%bound(i)%NFACE, mesh%NSHLBmax))
    allocate (mesh%bound(i)%F2E(mesh%bound(i)%NFACE))
    allocate (mesh%bound(i)%FACE_OR(mesh%bound(i)%NFACE))
    allocate (mesh%bound(i)%NSHLB(mesh%bound(i)%NFACE))
    do j = 1, mesh%bound(i)%NFACE
      mesh%bound(i)%NSHLB(j) = mesh%NSHLBmax
      read (mfid, *) (mesh%bound(i)%FACE_IEN(j, k), k=1, mesh%NSHLBmax)
    end do
    do j = 1, mesh%bound(i)%NFACE
      read (mfid, *) mesh%bound(i)%F2E(j), mesh%bound(i)%FACE_OR(j)
    end do

    allocate (mesh%bound(i)%BNODES(mesh%bound(i)%NNODE))
    do j = 1, mesh%bound(i)%NNODE
      read (mfid, *) mesh%bound(i)%BNODES(j)
    end do
  end do
  ! read nurbs data
  ! allocate (wg(NNODE))

  ! if (NPATCH > 0) then
  !   iga = .true.
  !   do i = 1, NNODE
  !     read (mfid, *) wg(i)
  !   end do

  !   allocate (EPID(NELEM), EIJK(NELEM, NSD))

  !   do i = 1, NELEM
  !     read (mfid, *) EPID(i), (EIJK(i, j), j=1, NSD)
  !   end do
  ! else
  !   iga = .false.
  !   if(allocated(wg)) then
  !     wg = 1.0d0
  !   endif
  ! end if

  ! read patches
  ! allocate (patch(NPATCH))
  ! do i = 1, NPATCH

  !   read (mfid, *) patch(i)%P, patch(i)%Q, patch(i)%R
  !   read (mfid, *) patch(i)%MCP, patch(i)%NCP, patch(i)%OCP

  !   allocate (patch(i)%U_KNOT(patch(i)%MCP + patch(i)%P + 1))
  !   allocate (patch(i)%V_KNOT(patch(i)%NCP + patch(i)%Q + 1))
  !   allocate (patch(i)%W_KNOT(patch(i)%OCP + patch(i)%R + 1))

  !   read (mfid, *) (patch(i)%U_KNOT(j), j=1, patch(i)%MCP + patch(i)%P + 1)
  !   read (mfid, *) (patch(i)%V_KNOT(j), j=1, patch(i)%NCP + patch(i)%Q + 1)
  !   read (mfid, *) (patch(i)%W_KNOT(j), j=1, patch(i)%OCP + patch(i)%R + 1)
  ! end do

  close (mfid)

  ! Init Flags
  ! allocate (IPER(NNODE))
  ! do i = 1, NNODE
  !   IPER(i) = i
  ! end do

  ! allocate (IBC(NNODE, 2*NSD + 2))
  ! IBC = 0
  ! allocate (IS_SOLID_NODE(NNODE))
  ! IS_SOLID_NODE_ASSIGNED = .false.

  ! allocate (EL_TYP(NELEM))
  ! EL_TYP = 0

  ! allocate (P_Flag(NNODE))
  ! P_Flag = 1

  ! allocate (D_Flag(NNODE))
  ! D_Flag = 0

  ! allocate (ELMNGAUSS(NELEM))
  ! ELMNGAUSS = ELMNSHL
  do i = 1, NBOUND
    allocate (mesh%bound(i)%NGAUSSB(mesh%bound(i)%NFACE))
    mesh%bound(i)%NGAUSSB = mesh%bound(i)%NSHLB
!     if(ismaster) write(*,*) i, bound(i)%NNODE
  end do

  !----------------------------------------------------------
  ! read the mapping between partitioned volume boundary
  ! and unpartitioned shell mesh (bmesh.*.dat)
  !----------------------------------------------------------
  write (iname, '(I8)') id
  fname = 'l2b.'//trim(adjustl(iname))//'.dat'
  open (mfid, file=fname, status='old')

  do i = 1, NBOUND

    if (mesh%bound(i)%NNODE > 0) then
      allocate (mesh%bound(i)%L2SNODE(mesh%bound(i)%NNODE))
      do j = 1, mesh%bound(i)%NNODE
        read (mfid, *) itmp1, mesh%bound(i)%L2SNODE(j)
      end do
    end if

    if (mesh%bound(i)%NFACE > 0) then
      allocate (mesh%bound(i)%L2SELEM(mesh%bound(i)%NFACE))

      do j = 1, mesh%bound(i)%NFACE
        ! read partitioned boundary element number and the corresponding
        ! shell mesh element number
!        read(mfid,*) itmp1, bound(i)%L2SELEM(j)
!        if (itmp1 /= bound(i)%F2E(j)) then
!          write(*,*) "ERROR in l2b: Boundary element doesn't match"
!          stop
!        end if
      end do
    end if

  end do
  close (mfid)

end subroutine input
