!======================================================================
!
!======================================================================
subroutine allocField(mesh, field)
  use class_def
  implicit none

  type(MeshData), intent(in) :: mesh
  type(FieldData), intent(inout) :: field

  integer :: NNODE, NSD
  NNODE = mesh%NNODE
  NSD = mesh%NSD

  ! Allocate solution vectors
  allocate (field%dg(NNODE, NSD), field%dgold(NNODE, NSD))

  allocate (field%ug(NNODE, NSD), field%ugold(NNODE, NSD), &
            field%ugm(NNODE, NSD), field%ugmold(NNODE, NSD))

  ! allocate (lhsgq(NNODE), rhsgq(NNODE))
  allocate (field%acg(NNODE, NSD), field%acgold(NNODE, NSD), &
            field%acgm(NNODE, NSD), field%acgmold(NNODE, NSD))

  allocate (field%phig(NNODE), field%phigold(NNODE), &
            field%rphig(NNODE), field%rphigold(NNODE))

  allocate (field%pg(NNODE), field%pgold(NNODE))

  allocate (field%rTg(NNODE), field%rTgold(NNODE))
  allocate (field%Tg(NNODE), field%Tgold(NNODE))

  
end subroutine allocField
!======================================================================
!
!======================================================================

subroutine freeField(field)
  use class_def
  implicit none

  type(FieldData), intent(inout) :: field

  deallocate (field%dg, field%dgold)
  deallocate (field%ug, field%ugold, field%ugm, field%ugmold)
  deallocate (field%acg, field%acgold, field%acgm, field%acgmold)
  deallocate (field%phig, field%phigold, field%rphig, field%rphigold)
  deallocate (field%pg, field%pgold)
  deallocate (field%rTg, field%rTgold)
  deallocate (field%Tg, field%Tgold)

end subroutine freeField

!======================================================================
!
!======================================================================

subroutine allocRHS(mesh, rhs)
  use class_def
  implicit none

  type(MeshData), intent(in) :: mesh
  type(RHSData), intent(inout) :: rhs

  rhs%N = mesh%NNODE
  rhs%NSD = mesh%NSD
  rhs%NVAR = mesh%NSD+3
  ! Allocate Residuals
  allocate (rhs%x(rhs%N, rhs%NVAR))
  ! allocate (rhs%RHSGU(NNODE, NSD))
  ! allocate (rhs%RHSGP(NNODE))
  ! allocate (rhs%RHSGLS(NNODE))
  ! allocate (rhs%RHSGTEM(NNODE))
end subroutine allocRHS


!======================================================================
!
!======================================================================
subroutine freeRHS(rhs)
  use class_def
  implicit none

  type(RHSData), intent(inout) :: rhs

  deallocate (rhs%x)
end subroutine FreeRHS
!======================================================================
!
!======================================================================
subroutine allocLHS(sp, mesh, lhs)
  use class_def
  implicit none

  type(SparsityPattern), intent(in) :: sp
  type(MeshData), intent(in) :: mesh

  type(LHSData), intent(inout) :: lhs

  integer :: nnz, NNODE, NSD
  nnz = sp%nnz
  NNODE = mesh%NNODE
  NSD = mesh%NSD

  ! Allocate Matrices
  allocate (lhs%LHSK11(NSD*NSD, nnz))
  allocate (lhs%LHSG(NSD, nnz))
  allocate (lhs%LHSD1(NSD, nnz))
  allocate (lhs%LHSM(nnz))

  allocate (lhs%LHSLS(nnz))
  allocate (lhs%LHSLSU(NSD, nnz))
  allocate (lhs%LHSULS(NSD, nnz))
  allocate (lhs%LHSPLS(nnz))

  allocate (lhs%LHSTEM(nnz))
end subroutine allocLHS

!======================================================================
!
!======================================================================

subroutine freeLHS(lhs)
  use class_def
  implicit none

  type(LHSData), intent(inout) :: lhs

  deallocate (lhs%LHSK11)
  deallocate (lhs%LHSG)
  deallocate (lhs%LHSD1)
  deallocate (lhs%LHSM)

  deallocate (lhs%LHSLS)
  deallocate (lhs%LHSLSU)
  deallocate (lhs%LHSULS)
  deallocate (lhs%LHSPLS)
end subroutine freeLHS
!======================================================================
!
!======================================================================

subroutine allocDirichletBC(mesh, bc)
  use class_def
  implicit none

  type(MeshData), intent(in) :: mesh
  type(DirichletBCData), intent(inout) :: bc

  integer :: NNODE, NSD, NBOUND
  NNODE = mesh%NNODE
  NSD = mesh%NSD
  NBOUND = mesh%NBOUND

  bc%NBOUND = NBOUND
  bc%NSD = NSD
  allocate (bc%BCugType(NBOUND, NSD))
  allocate (bc%BCugValu(NBOUND, NSD))

  allocate (bc%BCphigType(NBOUND))
  allocate (bc%BCphigValu(NBOUND))

  allocate (bc%BCTgType(NBOUND))
  allocate (bc%BCTgValu(NBOUND))

  allocate (bc%IBC(NNODE, NSD*2+2))
end subroutine allocDirichletBC
!======================================================================
!
!======================================================================

subroutine freeDirichletBC(bc)
  use class_def
  implicit none

  type(DirichletBCData), intent(inout) :: bc

  deallocate (bc%BCugType)
  deallocate (bc%BCugValu)

  deallocate (bc%BCphigType)
  deallocate (bc%BCphigValu)

  deallocate (bc%BCTgType)
  deallocate (bc%BCTgValu)

  deallocate (bc%IBC)
end subroutine freeDirichletBC
!======================================================================
!
!======================================================================
! subroutine allocMatVec(sp, mesh, field)
! 
!   ! use aAdjKeep
!   ! use commonvars
!   ! use defs_shell
!   use class_def
! 
!   implicit none
! 
!   type(SparsityPattern), intent(in) :: sp
!   type(MeshData), intent(in) :: mesh
! 
!   type(FieldData), intent(out) :: field
! 
!   integer :: NGAUSSMAX, i
!   integer :: icnt, NNODE
!   integer :: NSD
! 
!   NSD = mesh%NSD
!   NNODE = mesh%NNODE
!   icnt = sp%nnz
! 
! 
!   ! Allocate solution vectors
!   allocate (field%dg(NNODE, NSD), field%dgold(NNODE, NSD))
! 
!   allocate (field%ug(NNODE, NSD), field%ugold(NNODE, NSD), &
!             field%ugm(NNODE, NSD), field%ugmold(NNODE, NSD))
! 
!   ! allocate (lhsgq(NNODE), rhsgq(NNODE))
!   allocate (field%acg(NNODE, NSD), field%acgold(NNODE, NSD), &
!             field%acgm(NNODE, NSD), field%acgmold(NNODE, NSD))
! 
!   allocate (field%phig(NNODE), field%phigold(NNODE), &
!             field%rphig(NNODE), field%rphigold(NNODE))
! 
!   allocate (field%pg(NNODE), field%pgold(NNODE))
! 
!   allocate (field%rTg(NNODE), field%rTgold(NNODE))
!   allocate (field%Tg(NNODE), field%Tgold(NNODE))
! 
!   ! allocate (FPKS(NELEM, NSD, NSD))
! 
!   ! allocate averaged solution vectors
!   ! allocate (uavg(NNODE, NSD), pavg(NNODE))
! 
!   ! Allocate Residuals
!   allocate (RHSGu(NNODE, NSD), RHSGm(NNODE, NSD), &
!             RHSGp(NNODE), RHSGls(NNODE), RHSGTem(NNODE))
! 
!   ! Allocate Matrices
!   allocate (lhsK11(NSD*NSD, icnt))
!   allocate (lhsK12(NSD*NSD, icnt), lhsK22(NSD*NSD, icnt))
! 
!   allocate (lhsG(NSD, icnt), lhsD1(NSD, icnt), &
!             lhsD2(NSD, icnt))
! 
!   allocate (LHSM(icnt), LHSLS(icnt), LHSmass(icnt))
! 
!   allocate (LHSlsu(NSD, icnt), LHSUls(NSD, icnt), LHSPLS(icnt))
! 
!   allocate (LHSTem(icnt))
! 
!   ! NGAUSSMAX = maxval(ELMNGAUSS)
! 
! end subroutine allocMatVec

!======================================================================
!
!======================================================================
! subroutine deallocMatVec
! 
!   use aAdjKeep
!   use mpi
!   use commonvars
! 
!   implicit none
! 
!   integer :: i
! 
!   ! deallocate mesh and flags
!   deallocate (xg, wg, IEN)
!   deallocate (ELM_ID)
!   deallocate (IPER, IBC)
!   ! deallocate (IS_SOLID_NODE)
! 
!   deallocate (EL_TYP, P_Flag, D_Flag)
! 
!   !if (NPATCH .gt. 0) deallocate (EPID, EIJK)
!   !do i = 1, NPATCH
!   !  deallocate (patch(i)%U_KNOT)
!   !  deallocate (patch(i)%V_KNOT)
!   !  deallocate (patch(i)%W_KNOT)
!   !end do
!   !deallocate (patch)
! 
!   do i = 1, NBOUND
!     deallocate (bound(i)%FACE_IEN)
!     deallocate (bound(i)%F2E)
!     deallocate (bound(i)%FACE_OR)
!     deallocate (bound(i)%BNODES)
!   end do
!   deallocate (bound)
! 
!   ! deallocate sparseStruct
!   ! deallocate (col, row)
! 
!   ! dellocate solution vectors
!   deallocate (dg, dgold)
!   deallocate (ug, ugold, ugm, ugmold)
!   deallocate (acg, acgold, acgm, acgmold)
!   deallocate (phig, phigold, rphig, rphigold)
!   deallocate (pg, pgold)
!   deallocate (rTg, rTgold)
!   deallocate (Tg, Tgold)
!   ! deallocate (uavg, pavg)
!   ! deallocate (FPKS)
! 
!   ! deallocate Residuals
!   deallocate (RHSGu, RHSGm, RHSGp, RHSGls, RHSGtem)
! 
!   ! deallocate Matrices
!   deallocate (lhsK11, lhsK12, lhsK22)
!   deallocate (lhsG, lhsD1, lhsD2)
!   deallocate (LHSM, LHSLS, LHSmass)
!   deallocate (LHSlsu, LHSUls, lhsPls)
!   deallocate (LHStem)
! 
!   ! Conditional deallocations
!   if (numnodes .gt. 1) deallocate (ilworku)
! 
! end subroutine deallocMatVec
