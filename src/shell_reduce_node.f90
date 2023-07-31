!======================================================================
! program to read in all patches, remove duplicat nodes and creat
! mapping
! NP1 : total number of patches
! NP0 : number of patches without bending strips
!======================================================================
subroutine reduce_node(NSD, NP1, NNODZ, NEL, B_NET, maxNNODZ, mMap, &
                       redNNODZ, sumNNODZ, sumNEL)

  implicit none

  integer, intent(in) :: NSD, NP1, NNODZ(NP1), NEL(NP1), maxNNODZ
  integer, intent(out) :: sumNNODZ, sumNEL, redNNODZ, mMap(NP1, maxNNODZ)
  real(8), intent(in) :: B_NET(NP1, maxNNODZ, NSD + 1)

  real(8), parameter :: tol = 1.0d-6

  integer :: i, j, k, l, ip, ier
  integer :: DFlag
  real(8) :: dist

  real(8), allocatable :: xg(:, :), xl(:, :)

  character(len=30) :: fname

  ! determine maximum local NNODZ (for allocating arrays)
  ! also compute the total number of NNODZ (including duplicate nodes)
  sumNNODZ = sum(NNODZ)

  ! since elements will not overlap, the total number of elements should
  ! simply be sum(NEL)
  sumNEL = sum(NEL)

  !--------------------------------------------------------------------
  ! Now, loop through all the patches, read in mesh files, compare
  ! coordinates for overlaping nodes. Remove overlaping nodes and
  ! creat the mapping.
  ! Note: Loop through NP to 1. i.e., Higher numbering of patches is
  !   the master to the lower numbering. It has to be this way,
  !   since in our algorithm, master owns the node and all information
  !   on that node. Slave are zeros on overlapping nodes.
  !--------------------------------------------------------------------
  allocate (xg(sumNNODZ, NSD + 1), stat=ier)
  if (ier /= 0) stop 'Allocation Error: xg'
  xg = 0.0d0

  redNNODZ = 0
  do ip = NP1, 1, -1

!!!    write(*,*) '====== Patch', ip, '======'

    allocate (xl(NNODZ(ip), NSD + 1), stat=ier)
    if (ier /= 0) stop 'Allocation Error: U_KNOT'
    xl = 0.0d0

    do i = 1, NNODZ(ip)
      xl(i, :) = B_NET(ip, i, :)
    end do

    if ((i - 1) /= NNODZ(ip)) then
      write (*, *) "issue: NPC*MCP & NNODZ don't match"
    end if

    !--------------------------------------------------------------
    ! Now, check for duplicate nodes, or assign to global node
    ! NOTE: always loop backwards so as to maintain master slave relationship,
    !           in which the master has a higher number than the slave
    !--------------------------------------------------------------
    do i = 1, NNODZ(ip)

      DFlag = 0
      do j = 1, redNNODZ
        ! this should be fine even for the first node
        ! since xg(1:3) could be zero, but xg(4) should not
        dist = sqrt(sum((xg(j, 1:3) - xl(i, 1:3))**2))
        ! duplicate node found
        if (dist < tol) then
          DFlag = 1
          exit
        end if
      end do ! end loop global non-duplicat nodes (redNNODZ)

      ! New node found
      if (DFlag == 0) then
        redNNODZ = redNNODZ + 1
        xg(redNNODZ, :) = xl(i, :)
        mMap(ip, i) = redNNODZ
        ! allocate duplicate node
      else if (DFlag == 1) then
        mMap(ip, i) = j
      else
        stop 'wrong DFlag'
      end if
    end do ! end loop local nodes (NNODZ(ip))

    deallocate (xl)
  end do  ! end loop patches

  deallocate (xg)
end subroutine reduce_node
