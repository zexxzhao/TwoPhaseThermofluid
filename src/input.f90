!======================================================================
!
!======================================================================
subroutine input(id)
  use aAdjKeep
  use commonvars
  use mpi
  implicit none

  integer, intent(in) :: id

  character(len=30) :: fname, iname
  character(len=10) :: cname

  integer :: mfid, i, j, k, itmp1, itmp2, itmp3, nshl, counter

  counter = 1
  ! open mesh files
  mfid = 11
  fname = 'part'//trim(cname(id))//'.dat'
  open (mfid, file=fname, status='old')

  read (mfid, *) NSD, NSHL, NSHLBmax
  read (mfid, *) NNODE, NELEM, NBOUND, NPATCH

!  if (myid + 1 == 11) then
  !       write(*,*) counter, NSD, NSHL, NSHLBmax
  !      counter = counter + 1
  !     write(*,*) counter, NNODE, NELEM, NBOUND, NPATCH
  !    counter = counter + 1
!endif

  write (*, *) myid
  ! read nodes
  allocate (xg(NNODE, NSD), NodeID(NNODE))
  do i = 1, NNODE
    if (fem_flag == 1) then
      read (mfid, *) (xg(i, j), j=1, NSD), NodeID(i)
      if (myid + 1 == 1) then
        write (*, *) counter, xg(i, :)
        counter = counter + 1
      end if
    else
      read (mfid, *) (xg(i, j), j=1, NSD)
    end if
  end do
  write (*, *) myid
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)

  ! read elements
  allocate (ELMNSHL(NELEM), IEN(NELEM, NSHL))
  do i = 1, NELEM
    if (fem_flag == 1) then
      read (mfid, *) ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i)), itmp1
    else
      read (mfid, *) ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i))
    end if
    if ((id == 120) .and. (i == NELEM)) then
      write (*, *), ELMNSHL(i), (IEN(i, j), j=1, ELMNSHL(i))
    end if
  end do
  maxNSHL = maxval(ELMNSHL)
  ! read faces
  allocate (bound(NBOUND))
  do i = 1, NBOUND

    read (mfid, *) bound(i)%FACE_ID, bound(i)%NFACE, bound(i)%NNODE

    allocate (bound(i)%FACE_IEN(bound(i)%NFACE, NSHLBmax))
    allocate (bound(i)%F2E(bound(i)%NFACE))
    allocate (bound(i)%FACE_OR(bound(i)%NFACE))
    allocate (bound(i)%NSHLB(bound(i)%NFACE))

    do j = 1, bound(i)%NFACE
      bound(i)%NSHLB(j) = NSHLBmax
      if (fem_flag == 1) then
        read (mfid, *) (bound(i)%FACE_IEN(j, k), k=1, NSHLBmax)
      else
        read (mfid, *) bound(i)%NSHLB(j), (bound(i)%FACE_IEN(j, k), k=1, NSHLBmax)
      end if
    end do
    do j = 1, bound(i)%NFACE
      read (mfid, *) bound(i)%F2E(j), bound(i)%FACE_OR(j)
    end do

    allocate (bound(i)%BNODES(bound(i)%NNODE))
    do j = 1, bound(i)%NNODE
      read (mfid, *) bound(i)%BNODES(j)
    end do
  end do
  ! read nurbs data
  allocate (wg(NNODE))

  if (NPATCH > 0) then
    iga = .true.
    do i = 1, NNODE
      read (mfid, *) wg(i)
    end do

    allocate (EPID(NELEM), EIJK(NELEM, NSD))

    do i = 1, NELEM
      read (mfid, *) EPID(i), (EIJK(i, j), j=1, NSD)
    end do
  else
    iga = .false.
    wg = 1.0d0
  end if

  ! read patches
  allocate (patch(NPATCH))
  do i = 1, NPATCH

    read (mfid, *) patch(i)%P, patch(i)%Q, patch(i)%R
    read (mfid, *) patch(i)%MCP, patch(i)%NCP, patch(i)%OCP

    allocate (patch(i)%U_KNOT(patch(i)%MCP + patch(i)%P + 1))
    allocate (patch(i)%V_KNOT(patch(i)%NCP + patch(i)%Q + 1))
    allocate (patch(i)%W_KNOT(patch(i)%OCP + patch(i)%R + 1))

    read (mfid, *) (patch(i)%U_KNOT(j), j=1, patch(i)%MCP + patch(i)%P + 1)
    read (mfid, *) (patch(i)%V_KNOT(j), j=1, patch(i)%NCP + patch(i)%Q + 1)
    read (mfid, *) (patch(i)%W_KNOT(j), j=1, patch(i)%OCP + patch(i)%R + 1)
  end do

  close (mfid)

  ! Init Flags
  allocate (IPER(NNODE))
  do i = 1, NNODE
    IPER(i) = i
  end do

  allocate (IBC(NNODE, 2*NSD + 2))
  IBC = 0

  allocate (EL_TYP(NELEM))
  EL_TYP = 0

  allocate (P_Flag(NNODE))
  P_Flag = 1

  allocate (D_Flag(NNODE))
  D_Flag = 0

  allocate (ELMNGAUSS(NELEM))
  ELMNGAUSS = ELMNSHL
  do i = 1, NBOUND
    allocate (bound(i)%NGAUSSB(bound(i)%NFACE))
    bound(i)%NGAUSSB = bound(i)%NSHLB
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

    if (bound(i)%NNODE > 0) then
      allocate (bound(i)%L2SNODE(bound(i)%NNODE))
      do j = 1, bound(i)%NNODE
        read (mfid, *) itmp1, bound(i)%L2SNODE(j)
      end do
    end if

    if (bound(i)%NFACE > 0) then
      allocate (bound(i)%L2SELEM(bound(i)%NFACE))

      do j = 1, bound(i)%NFACE
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
