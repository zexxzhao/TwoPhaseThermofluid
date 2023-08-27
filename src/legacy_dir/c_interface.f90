!======================================================================
! C
!======================================================================

subroutine allocField_C(mesh_c, field_c) bind(C, name="allocField_C")
  use iso_c_binding

  use class_def_c
  implicit none

  type(FieldData_C), pointer :: field_c
  type(FieldData), allocatable, target :: field
  type(MeshData_C), pointer :: mesh_c
  type(MeshData), pointer :: mesh

  ! Assciate the pointer
  call c_f_pointer(mesh_c%data_f, mesh)

  ! Allocate the instance
  call allocField(mesh, field)

  ! Associate the pointer
  field_c%data_f = c_loc(field)

end subroutine allocField_C 

!======================================================================
subroutine freeField_C(field_c) bind(C, name="freeField_C")
  use iso_c_binding

  use class_def_c
  implicit none

  type(FieldData_C), pointer :: field_c
  type(FieldData), pointer :: field

  ! Assciate the pointer
  call c_f_pointer(field_c%data_f, field)

  ! Free the instance
  call freeField(field)

end subroutine freeField_C

!======================================================================
! C
!======================================================================

subroutine allocRHS_C(mesh_c, rhs_c) bind(C, name="allocRHS_C")
  use iso_c_binding

  use class_def_c
  implicit none

  type(RHSData_C), pointer :: rhs_c
  type(RHSData), allocatable, target :: rhs
  type(MeshData_C), pointer :: mesh_c
  type(MeshData), pointer :: mesh

  ! Assciate the pointer
  call c_f_pointer(mesh_c%data_f, mesh)

  ! Allocate the instance
  call allocRHS(mesh, rhs)

  ! Associate the pointer
  rhs_C%data_f = c_loc(rhs)

end subroutine allocRHS_C

!======================================================================

subroutine freeRHS_C(rhs_c) bind(C, name="freeRHS_C")
  use iso_c_binding

  use class_def_c
  implicit none

  type(RHSData_C), pointer :: rhs_c
  type(RHSData), pointer :: rhs

  ! Assciate the pointer
  call c_f_pointer(rhs_c%data_f, rhs)

  ! Free the instance
  call freeRHS(rhs)

end subroutine freeRHS_C

!======================================================================
! C
!======================================================================

subroutine allocLHS_C(sp_c, mesh_c, lhs_c) bind(C, name="allocLHS_C")
    use iso_c_binding
    use class_def_c
    implicit none

    type(LHSData_C), pointer :: lhs_c
    type(LHSData), allocatable, target :: lhs
    type(MeshData_C), pointer :: mesh_c
    type(MeshData), pointer :: mesh
    type(SparsityPattern_C), pointer :: sp_c
    type(SparsityPattern), pointer :: sp

    ! Assciate the pointer
    call c_f_pointer(mesh_c%data_f, mesh)
    call c_f_pointer(sp_c%data_f, sp)

    ! Allocate the instance
    call allocLHS(sp, mesh, lhs)

    ! Associate the pointer
    lhs_c%data_f = c_loc(lhs)


end subroutine allocLHS_C

!======================================================================

subroutine freeLHS_C(lhs_c) bind(C, name="freeLHS_C")
    use iso_c_binding
    use class_def_c
    implicit none

    type(LHSData_C), pointer :: lhs_c
    type(LHSData), pointer :: lhs

    ! Assciate the pointer
    call c_f_pointer(lhs_c%data_f, lhs)

    ! Free the instance
    call freeLHS(lhs)

end subroutine freeLHS_C

!======================================================================
! C
!======================================================================

subroutine allocDirichletBC_C(mesh_c, bc_c) bind(C, name='allocDirichletBC_C')
    use iso_c_binding
    use class_def_c

    implicit none

    type(DirichletBCData_C), pointer :: bc_c
    type(DirichletBCData), allocatable, target :: bc
    type(MeshData_C), pointer :: mesh_c
    type(MeshData), pointer :: mesh

    ! Assciate the pointer
    call c_f_pointer(mesh_c%data_f, mesh)

    ! Allocate the instance
    call allocDirichletBC(mesh, bc)

    ! Associate the pointer
    bc_c%data_f = c_loc(bc)

end subroutine allocDirichletBC_C

!======================================================================

subroutine freeDirichletBC_C(bc_c) bind(C, name='freeDirichletBC_C')
    use iso_c_binding
    use class_def_c

    implicit none

    type(DirichletBCData_C), pointer :: bc_c
    type(DirichletBCData), pointer :: bc

    ! Assciate the pointer
    call c_f_pointer(bc_c%data_f, bc)

    ! Free the instance
    call freeDirichletBC(bc)

end subroutine freeDirichletBC_C