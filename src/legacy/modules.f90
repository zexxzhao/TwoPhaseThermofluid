!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
module class_def
  implicit none

  ! type NURBSpatch

  !   integer :: P, Q, R
  !   integer :: MCP, NCP, OCP
  !   real(8), allocatable :: U_KNOT(:), V_KNOT(:), W_KNOT(:)

  ! end type NURBSpatch
  type bnd_class

    integer :: FACE_ID

    integer :: NFACE

    integer, allocatable :: FACE_IEN(:, :)

    integer, allocatable :: F2E(:)
    integer, allocatable :: FACE_OR(:)
    integer, allocatable :: NSHLB(:)
    integer, allocatable :: NGAUSSB(:)

    integer :: NNODE
    integer, allocatable :: BNODES(:)

    ! mapping between partitioned (local) boundary node/element
    ! and unpartitioned boundary (shell) node/element
    integer, allocatable :: L2SNODE(:), L2SELEM(:)

  end type bnd_class

  type MeshData
    integer :: NSD, NNODE, NELEM, NBOUND
    integer :: NPATCH, NSHLBmax, maxNSHL
    real(8), allocatable :: xg(:, :)
  
    integer, allocatable :: IEN(:, :), NodeID(:)
    integer, allocatable :: ELM_ID(:)

    type(bnd_class), allocatable :: bound(:)
    integer, allocatable :: ELMNSHL(:), ELMNGAUSS(:)
  end type MeshData

  type SparsityPattern
    integer :: NNODE
    integer :: nnz ! number of non-zero entries
    integer, allocatable :: index(:) ! size(index) = nnz
    integer, allocatable :: indptr(:) ! size(indptr) = NNODE+1
  end type SparsityPattern

  type FieldData
    real(8), allocatable :: dg(:, :), dgold(:, :)
    real(8), allocatable :: ugm(:, :), ugmold(:, :)
    real(8), allocatable :: acgm(:, :), acgmold(:, :)

    real(8), allocatable :: ug(:, :), ugold(:, :)
    real(8), allocatable :: acg(:, :), acgold(:, :)
    real(8), allocatable :: pg(:), pgold(:)

    real(8), allocatable :: phig(:), phigold(:)
    real(8), allocatable :: rphig(:), rphigold(:)

    real(8), allocatable :: rTg(:), rTgold(:)
    real(8), allocatable :: Tg(:), Tgold(:)

  end type FieldData

  type RHSData
    real(8), allocatable :: RHSGU(:, :)
    real(8), allocatable :: RHSGP(:)
    real(8), allocatable :: RHSGLS(:)
    real(8), allocatable :: RHSGTEM(:)
  end type RHSData

  type LHSData
    real(8), allocatable :: LHSK11(:, :)
    real(8), allocatable :: LHSG(:, :)
    real(8), allocatable :: LHSD1(:, :)
    real(8), allocatable :: LHSM(:)
    real(8), allocatable :: LHSLS(:)
    real(8), allocatable :: LHSLSU(:, :)
    real(8), allocatable :: LHSPLS(:)
    real(8), allocatable :: LHSULS(:, :)
    real(8), allocatable :: LHSTEM(:)
  end type LHSData

  type DirichletBCData
    integer :: NBOUND, NSD
    integer, allocatable :: BCugType(:,:)
    real(8), allocatable :: BCugValu(:,:)
    integer, allocatable :: BCphigType(:)
    real(8), allocatable :: BCphigValu(:)
    integer, allocatable :: BCTgType(:)
    real(8), allocatable :: BCTgValu(:)

    integer, allocatable :: IBC(:, :)
  end type DirichletBCData

end module class_def
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
module class_def_c
  use iso_c_binding
  use class_def
  implicit none

  type, bind(C) :: bnd_class_C
    integer(C_INT) :: FACE_ID

    integer(C_INT) :: NFACE

    type(C_PTR) :: FACE_IEN

    type(C_PTR) :: F2E
    type(C_PTR) :: FACE_OR
    type(C_PTR) :: NSHLB
    type(C_PTR) :: NGAUSSB

    integer(C_INT) :: NNODE
    type(C_PTR) :: BNODES

    ! mapping between partitioned (local) boundary node/element
    ! and unpartitioned boundary (shell) node/element
    type(C_PTR) :: L2SNODE, L2SELEM
  end type bnd_class_C

  type, bind(C) :: MeshData_C
    integer(C_INT) :: NSD, NNODE, NELEM, NBOUND
    integer(C_INT) :: NPATCH, NSHLBmax, maxNSHL
    type(C_PTR) :: xg
  
    type(C_PTR) :: IEN, NodeID
    type(C_PTR) :: ELM_ID

    type(C_PTR) :: bound
    type(C_PTR) :: ELMNSHL, ELMNGAUSS

    type(C_PTR) :: mesh_data_f
  end type MeshData_C

  type, bind(C) :: SparsityPattern_C
    integer(C_INT) :: NNODE
    integer(C_INT) :: nnz ! number of non-zero entries
    type(C_PTR) :: index ! size(index) = nnz
    type(C_PTR) :: indptr ! size(indptr) = NNODE+1
  end type SparsityPattern_C

  type, bind(C) :: FieldData_C
    type(C_PTR) :: dg, dgold
    type(C_PTR) :: ugm, ugmold
    type(C_PTR) :: acgm, acgmold

    type(C_PTR) :: ug, ugold
    type(C_PTR) :: acg, acgold
    type(C_PTR) :: pg, pgold

    type(C_PTR) :: phig, phigold
    type(C_PTR) :: rphig, rphigold

    type(C_PTR) :: rTg, rTgold
    type(C_PTR) :: Tg, Tgold
  end type FieldData_C

  type, bind(C) :: RHSData_C
    type(C_PTR) :: RHSGU
    type(C_PTR) :: RHSGP
    type(C_PTR) :: RHSGLS
    type(C_PTR) :: RHSGTEM
  end type RHSData_C

  type, bind(C) :: LHSData_C
    type(C_PTR) :: LHSK11
    type(C_PTR) :: LHSG
    type(C_PTR) :: LHSD1
    type(C_PTR) :: LHSM
    type(C_PTR) :: LHSLS
    type(C_PTR) :: LHSLSU
    type(C_PTR) :: LHSPLS
    type(C_PTR) :: LHSULS
    type(C_PTR) :: LHSTEM
  end type LHSData_C

  type, bind(C) :: DirichletBCData_C
    integer(C_INT) :: NBOUND, NSD
    type(C_PTR) :: BCugType
    type(C_PTR) :: BCugValu
    type(C_PTR) :: BCphigType
    type(C_PTR) :: BCphigValu
    type(C_PTR) :: BCTgType
    type(C_PTR) :: BCTgValu

    type(C_PTR) :: IBC
  end type DirichletBCData_C

end module class_def_c

!------------------------------------------------------------------------
! Module for defining shell types and variables
!------------------------------------------------------------------------
! module defs_shell
! 
!   implicit none
! 
!   ! Declare type mesh
!   type :: mesh
!     ! degree in U and V for each patch
!     integer, allocatable :: P(:), Q(:)
! 
!     ! patch type
!     ! 1-blade; 0-bending strips; 2-shear web; ...
!     integer, allocatable :: PTYPE(:)
! 
!     ! Knot vectors in u and v directions for each element
!     real(8), allocatable :: U_KNOT(:, :), V_KNOT(:, :)
! 
!     ! Size of the knot vector for each elements (e.g. NUK=P+MCP+1)
!     integer, allocatable :: NUK(:), NVK(:)
! 
!     ! Control Net
!     ! B_NET is reference config, B_NET_D is current config.
!     ! For the pre-bend case, reference config. could change, so
!     ! we created B_NET_U for undeformed, original config.
!     real(8), allocatable :: B_NET(:, :), B_NET_U(:, :), B_NET_D(:, :)
! 
!     ! Boundary condition indicator for global nodes and edges
!     ! respectively
!     integer, allocatable :: IBC(:, :)
! 
!     ! array to store force vectors on the wind turbine blade
!     real(8), allocatable :: FORCE(:, :)
! 
!     ! Global connectivity array
!     integer, allocatable :: IEN(:, :), INN(:, :)
! 
!     ! number of shape functions for every element
!     integer, allocatable :: NSHL(:)
! 
!     ! Bezier extraction operator
!     real(8), allocatable :: Ext(:, :, :)
! 
!     ! Array for closest points (and the corresponding element)
!     integer, allocatable :: CLE(:, :)
!     real(8), allocatable :: CLP(:, :, :)
! 
!     ! Array for element list
!     integer, allocatable :: Elm_Close(:, :, :)
!     integer :: NEL_Close, NEL_Close_t
! 
!     ! Array for element gauss point location and for buiding the
!     ! element list based on radial location
!     real(8) :: Elm_Size
!     real(8), allocatable :: Elm_Loc(:, :, :)
!     integer, allocatable :: RAD_ELM_LIST(:, :, :), RAD_ELM_NUM(:, :)
! 
!     integer :: NGAUSS, NNODE, NEL, maxNSHL
! 
!     ! Store node number of the tip
!     integer :: TipLoc, TipLocTr
! 
!     ! array for Solution vectors
!     real(8), allocatable :: dsh(:, :), dshold(:, :), &
!                             ush(:, :), ushold(:, :), &
!                             ash(:, :), ashold(:, :)
! 
!     ! Surface ID and Boundary number
!     integer :: FaceID, iBound
!   end type mesh
! 
!   ! Declare type mesh for multi-patch
!   type :: mesh_mp
!     ! degree in U and V for each patch
!     integer, allocatable :: P(:), Q(:)
! 
!     ! number of control points in U and V for each patch
!     ! (no need for T-spline)
!     integer, allocatable :: MCP(:), NCP(:)
! 
!     ! Total number of control points and elements for each patch
!     integer, allocatable :: NNODE(:), NEL(:)
! 
!     ! patch type
!     ! 1-blade surface; 0-bending strips; 2-shear web; ...
!     integer, allocatable :: PTYPE(:)
! 
!     ! Knot vectors in u and v directions for each patch
!     real(8), allocatable :: U_KNOT(:, :), V_KNOT(:, :)
! 
!     ! Control Net
!     real(8), allocatable :: B_NET(:, :, :)
! 
!     ! Boundary condition indicator for global nodes and edges
!     ! respectively
!     integer, allocatable :: IBC(:, :, :)
! 
!     ! array to store force vectors on the wind turbine blade
!     real(8), allocatable :: FORCE(:, :, :)
! 
!     ! Mapping between patches and global reduced numbering
!     ! e.g., MAP(1,2) = 3 means 1st patch, 2nd node points to
!     !   global reduced node number 3
!     integer, allocatable :: MAP(:, :)
!   end type mesh_mp
! 
!   ! Declare type shell (for wind turbine blade)
!   type :: shell_bld
! 
!     type(mesh_mp) :: mNRB
!     type(mesh)    :: NRB
! 
!     type(mesh)    :: TSP, BEZ
! 
!     type(mesh)    :: FEM
! 
!     ! number of patches for Blade Surface (S). Matches the fluid mesh
!     ! number of patches for blade structure (B). May include shear webs
!     ! number of total patches including bending strips (T)
!     integer :: NPS, NPB, NPT
! 
!     integer :: M_Flag, T_Flag
!     real(8) :: RHSGtol, G_fact(3), RHSGNorm
! 
!     ! row, col, and total of nonzero entries for sparse structure
!     integer, allocatable :: row(:), col(:)
!     integer :: icnt
! 
!     ! The right hand side load vector G and left hand stiffness matrix K
!     real(8), allocatable :: RHSG(:, :), LHSK(:, :), &
!                             RHSG_EXT(:, :), RHSG_GRA(:, :)
! 
! !    ! Solution vectors
! !    real(8), allocatable :: yg(:,:), dg(:,:), tg(:), &
! !                            mg(:), dl(:,:)
! 
!     ! material matrix for composite
!     integer :: NMat
! !    real(8), allocatable :: matA(:,:,:), matB(:,:,:), matD(:,:,:)
!     real(8), allocatable :: matA(:, :, :, :), matB(:, :, :, :), matD(:, :, :, :), Density(:, :), Thickness(:, :)
! 
!     ! number of newton iterations for shell
!     integer, allocatable :: Nnewt(:)
! 
!     ! Torque computed on shell mesh
!     real(8) :: Tq1, Tq2
! 
!     ! Blade rotation. 0 degree is the straight-up position
!     real(8) :: BldRot
! 
!     integer :: bmap
!   end type shell_bld
! 
!   ! Declare type shell (for non-matching boundaries)
!   type :: shell_nmb
! 
!     type(mesh), allocatable :: FEM(:)
! 
!   end type shell_nmb
! end module defs_shell

!------------------------------------------------------------------------
!     Module for storing arrays and allocation routines
!------------------------------------------------------------------------
module aAdjKeep

  use class_def
  ! use defs_shell

  implicit none
  save

  ! Mesh
  real(8), allocatable :: xg(:, :), wg(:)

  integer, allocatable :: IEN(:, :), NodeID(:)
  integer, allocatable :: EPID(:), EIJK(:, :)
  integer, allocatable :: ELM_ID(:)

  type(bnd_class), allocatable :: bound(:)
  !type(NURBSpatch), allocatable :: patch(:)
  ! Array for Prism
  integer, allocatable :: ELMNSHL(:), ELMNGAUSS(:)

  ! Contraint flags
  integer, allocatable :: IPER(:)
  integer, allocatable :: IBC(:, :)
  ! logical :: IS_SOLID_NODE_ASSIGNED
  ! integer, allocatable :: IS_SOLID_NODE(:)

  ! Type flags
  integer, allocatable :: EL_TYP(:), D_FLAG(:), P_FLAG(:)

  ! Spars Struc
  ! integer, allocatable :: row(:), col(:)

  integer :: InNNODE1, InCELL
  integer, allocatable :: ien_inlet(:, :)
  real(8), allocatable :: xg_inlet(:, :)

  ! The right hand side load vector G and left hand stiffness matrix K
  real(8), allocatable :: RHSGU(:, :), RHSGM(:, :), RHSGP(:), &
                          RHSGLS(:), RHSGTEM(:)

  real(8), allocatable :: lhsK11(:, :), lhsK12(:, :), lhsK22(:, :), &
                          lhsG(:, :), lhsD1(:, :), lhsD2(:, :), &
                          lhsM(:), lhsLS(:), LHSPi(:, :), &
                          LHSPti(:, :), LHSKi(:, :), lhsMi(:), &
                          LHSLSi(:), LHSmass(:), &
                          LHSlsu(:, :), LHSPls(:), LHSlsP(:), LHSUls(:, :), &
                          LHSTem(:)

  real(8), allocatable :: lhsgq(:), rhsgq(:)
  ! Solution vectors
  real(8), allocatable :: dg(:, :), dgold(:, :), &
                          ug(:, :), ugold(:, :), &
                          ugm(:, :), ugmold(:, :), &
                          acg(:, :), acgold(:, :), &
                          acgm(:, :), acgmold(:, :), &
                          pg(:), pgold(:), &
                          phig(:), phigold(:), &
                          rphig(:), rphigold(:), &
                          rTg(:), rTgold(:), &
                          Tg(:), Tgold(:)

  ! real(8), allocatable :: uavg(:, :), pavg(:)

  ! Rigid body
  ! real(8) :: vbn0(3), vbn1(3)
  ! real(8) :: dbn0(3), dbn1(3)
  ! real(8) :: wbn0(3), wbn1(3)
  ! real(8) :: Rn0(3, 3), Rn1(3, 3)

  ! First P-K Stress
  ! real(8), allocatable :: FPKS(:, :, :)


  ! global information for individual blades
  ! type(mesh) :: blade(3)

  ! real(8) :: Center_Rot(3)

end module aAdjKeep

!----------------------------------------------------------------------
! Module for common variables
!----------------------------------------------------------------------
module commonvars

  implicit none
  save

  integer::check_num, InNNODE
  integer::index_plane_p, index_plane_n, max_planes
  real(8)::Delt_x, x_inflow, remain_x, position_init
  integer::tmp_num1, tmp_num2
  real(8), dimension(1:3)::C_tmp1, C_tmp2
  real(8)::V_tmp1, V_tmp2
  real(8), allocatable :: inflow_velocity(:, :), inflow_pressure(:), inflow_rho(:)
  real(8), allocatable :: inflow_x(:), inflow_y(:), inflow_z(:)
  integer, allocatable :: inflow_bnodes(:)

  real(8), allocatable :: inflow_velocity_check(:, :)
  character(len=1000) :: filename_p, filename_n, filename_check1, filename_check2
  real(8), dimension(1:3)::Inflow_velocity_p, Inflow_velocity_n
  logical::turbulence_flag
  real(8)::interpolation_factor, inflow_pressure_p, inflow_pressure_n, inflow_rho_p, inflow_rho_n
  real(8)::top_flag
  ! Mesh
  integer :: NSD, NNODE, NELEM, NBOUND, NPATCH, NSHLBmax, &
             NBlade, maxNSHL
  logical :: iga
  integer :: use_vms
  real(8) :: DetJ, DetJb, DetJinv, hglob

  ! Time step
  real(8) :: Delt, Dtgl, rhoinf, beti, gami, alfi, almi, &
             mgam, ogam, time, &
             conv_time, mono_time, move_time, shel_time
  ! real(8) :: lambda
  integer :: Nstep, ifq, ifq_sh, ifq_tq, mono_iter

  ! Navier-Stokes solver
  real(8) :: mua, rhoa, muw, rhow
  real(8) :: cpa, cpw, kappaa, kappaw
  real(8) :: rhos, mus, cps, kappas

  real(8) :: NS_kdc_w, NS_kdc_a, fine_tau

  real(8) :: NS_GMRES_atol, NS_NL_Uatol, NS_NL_Patol
  real(8) :: NS_GMRES_tol, NS_NL_Utol, NS_NL_Ptol
  integer :: NS_GMRES_itermin, NS_GMRES_itermax, NS_NL_itermax, &
             NS_hess_flag

  ! LevelSet Convection solver
  real(8) :: LSC_kdc
  real(8) :: LSC_GMRES_atol, LSC_NL_atol
  real(8) :: LSC_GMRES_tol, LSC_NL_tol
  integer :: LSC_GMRES_itermin, LSC_GMRES_itermax, LSC_NL_itermax, &
             LSC_pred_step

  ! Temperature Advection-diffusion solver
  real(8) :: TEM_kdc
  real(8) :: TEM_GMRES_atol, TEM_NL_atol
  real(8) :: TEM_GMRES_tol, TEM_NL_tol
  integer :: TEM_GMRES_itermin, TEM_GMRES_itermax, TEM_NL_itermax

  ! Rigid body solver
  real(8) :: RB_NL_Ftol, RB_NL_Mtol

  ! Mesh solver
  real(8) :: Mesh_GMRES_tol, Mesh_NL_tol
  integer :: Mesh_GMRES_itermin, Mesh_GMRES_itermax
  integer :: fem_flag
  ! Level Set Redistance solver
  real(8) :: LSRD_kdc, LSRD_penalty_fac
  real(8) :: LSRD_pseudoDTGL_fac
  real(8) :: LSRD_GMRES_tol, LSRD_NL_tol
  integer :: LSRD_GMRES_itermin, LSRD_GMRES_itermax, LSRD_NL_itermax

  ! Level Set massfix
  real(8) :: Mass_init
  real(8) :: Mass_NL_tol
  integer :: Mass_NL_itermax

  ! Rigid body
  real(8) :: massb, Ibhat(3, 3), xcg(3), Fb(3), Mb(3)
  integer :: VBC(3), MBC(3)

  ! Domain and Hull
  real(8) :: domain_top, domain_bottom, &
             domain_left, domain_right, &
             domain_front, domain_back, &
             water_level, water_depth, air_height, hull_length

  ! Wave generating wall
  real(8) :: wave_periode, wave_length, wave_amp, wave_angle

  ! Gravity
  real(8) :: gravity, gravvec(3)

  ! Interface
  real(8) :: mp_eps, beta_t, phi_inf, &
             phi_t, C_t, Re, Pe, Gr, Fr, &
             Ra, Sr, Pr, dphi_bg(3), cross_flag

  !real(8) :: mu, rho, dmudphi, drhodphi
  real(8) :: kappa
  !real(8) :: cp, hk
  real(8) :: Ts, c_cond, c_evap, lh
  real(8) :: C_DMDOT
  ! Setup
  real(8) :: Froude, Uin
  integer :: BCtype(99)

  integer, allocatable :: BCugType(:, :)
  real(8), allocatable :: BCugValu(:, :), phi_bg(:)
  integer, allocatable :: BCphigType(:), BCTgType(:)
  real(8), allocatable :: BCphigValu(:), BCTgValu(:)

  ! real(8) :: usettle
  ! Flags
  ! logical :: move, mono, conv, shel

  ! rotation
  real(8) :: theta, thetd, thedd, thetaOld, thetdOld, theddOld, &
             maxthetd, torque1, torque2, torque3, torque4, torque_SH, &
             moment1, moment2, force_trac(3), force_cons(3)

  logical :: solshell, nonmatch
end module commonvars

!----------------------------------------------------------------------
! Module for common parameters (e.g. pi)
!----------------------------------------------------------------------
module commonpars
  implicit none
  save

  real(8), parameter :: pi = 3.14159265358979323846264338328d0
  real(8), parameter :: gravvec(3) = (/0.d0, 0.d0, -9.81d0/)

  ! Assemble field
  integer, parameter :: ASSEMBLE_FIELD_NONE = 0
  integer, parameter :: ASSEMBLE_FIELD_NS = 1
  integer, parameter :: ASSEMBLE_FIELD_LS = 2
  integer, parameter :: ASSEMBLE_FIELD_VOF = ASSEMBLE_FIELD_LS
  integer, parameter :: ASSEMBLE_FIELD_TEM = 4

  ! Assemble tensor
  integer, parameter :: ASSEMBLE_TENSOR_NONE = 0
  integer, parameter :: ASSEMBLE_TENSOR_SCALAR = 1
  integer, parameter :: ASSEMBLE_TENSOR_VEC = 2
  integer, parameter :: ASSEMBLE_TENSOR_MAT = 4
end module commonpars

!----------------------------------------------------------------------
!     Module for mpi
!----------------------------------------------------------------------
module mpi
  implicit none
  save

  include "mpif.h"

  integer, parameter :: mpi_master = 0
  integer, parameter :: maxtask = 50
  integer, parameter :: maxseg = 15000

  integer :: numnodes, myid, mpi_err
  integer :: status(MPI_STATUS_SIZE)
  integer :: lstrideu, lstridep, numtask
  integer :: lfrontu, maxfrontu
  integer :: lfrontp, maxfrontp
  integer :: nlworku, nlworkp
  integer :: itag, iacc, iother, numseg, isgbeg, itkbeg, isgend
  integer :: sevsegtypeu(maxtask, 15)
  integer :: sevsegtypep(maxtask, 15)
  logical :: ismaster

  integer, allocatable :: ilworku(:), ilworkp(:)
end module mpi

module configuration
  use iso_c_binding

  type VMSConfigType
    logical :: use_vms
    logical :: use_taubar
    logical :: use_sliding_velocity
    real(8) :: NS_kdc_w, NS_kdc_a
    real(8) :: LSC_kdc
    real(8) :: Tem_kdc
    real(8) :: c_dmdot

  end type VMSConfigType

  type KSPConfigType

    integer :: NRES
    integer, allocatable :: max_iter(:), min_iter(:)
    real(8), allocatable :: atol(:), rtol(:)
  end type KSPConfigType

  type NewtonRaphsonConfigType

    integer :: NRES
    integer :: max_iter, min_iter
    real(8), allocatable :: atol(:), rtol(:)
  end type NewtonRaphsonConfigType

  type TimeIntegralConfigType
    integer :: Nstep
    integer :: ifq
    real(8) :: Delt
    real(8) :: rhoinf
  end type TimeIntegralConfigType

  type PropertyType
    real(8) :: rhoa, rhow, rhos
    real(8) :: mua, muw, mus
    real(8) :: cpa, cpw, cps
    real(8) :: kappaa, kappaw, kappas
    real(8) :: Ts, lh, c_cond, c_evap
  end type PropertyType

  type BCConfigType
    integer :: NBOUND, NSD
    ! ug
    integer, allocatable :: BCugType(:,:)
    real(8), allocatable :: BCugValu(:,:)
    ! phig
    integer, allocatable :: BCphigType(:)
    real(8), allocatable :: BCphigValu(:)
    ! Tg
    integer, allocatable :: BCTgType(:)
    real(8), allocatable :: BCTgValu(:)
  end type BCConfigType

  type MPIConfigType
    integer :: numnodes, myid
    logical :: ismaster
  end type MPIConfigType

  type ConfigType
    logical :: iga
    logical :: fem_flag
    logical :: use_hessian
    logical :: calc_cfl
    type(TimeIntegralConfigType) :: time_integral
    type(PropertyType) :: property
    type(BCConfigType) :: bc
    type(VMSConfigType) :: vms
    type(KSPConfigType) :: ksp
    type(NewtonRaphsonConfigType) :: newton_raphson
    type(MPIConfigType) :: mpi
  end type ConfigType

  contains
  subroutine init_config(config)
    use aAdjKeep
    use commonvars
  
    implicit none
  
    include "mpif.h"
    type(ConfigType), intent(out) :: config

    integer :: NRES = 4
    integer :: mpi_err

    config%iga = .false.
    config%fem_flag = fem_flag /= 0
    config%use_hessian = .false.
    config%calc_cfl = .true.

    config%time_integral%Nstep = Nstep
    config%time_integral%ifq = ifq
    config%time_integral%Delt = Delt
    config%time_integral%rhoinf = rhoinf

    config%property%rhoa = rhoa
    config%property%rhow = rhow
    config%property%mua = mua
    config%property%muw = muw
    config%property%cpa = cpa
    config%property%cpw = cpw
    config%property%kappaa = kappaa
    config%property%kappaw = kappaw
    config%property%rhos = rhos
    config%property%cps = cps
    config%property%kappas = kappas
    config%property%Ts = Ts
    config%property%lh = lh
    config%property%c_cond = c_cond
    config%property%c_evap = c_evap

    config%bc%NBOUND = NBOUND
    config%bc%NSD = NSD
    allocate (config%bc%BCugType(config%bc%NBOUND, config%bc%NSD), &
              config%bc%BCugValu(config%bc%NBOUND, config%bc%NSD), &
              config%bc%BCphigType(config%bc%NBOUND), &
              config%bc%BCphigValu(config%bc%NBOUND), &
              config%bc%BCTgType(config%bc%NBOUND), &
              config%bc%BCTgValu(config%bc%NBOUND))
    config%bc%BCugType = BCugType
    config%bc%BCugValu = BCugValu
    config%bc%BCphigType = BCphigType
    config%bc%BCphigValu = BCphigValu
    config%bc%BCTgType = BCTgType
    config%bc%BCTgValu = BCTgValu

    config%vms%use_vms = USE_VMS /= 0
    config%vms%use_taubar = .false.
    config%vms%use_sliding_velocity = .false.
    config%vms%NS_kdc_w = NS_kdc_w
    config%vms%NS_kdc_a = NS_kdc_a
    config%vms%LSC_kdc = LSC_kdc
    config%vms%Tem_kdc = TEM_kdc
    config%vms%c_dmdot = C_DMDOT
  
    config%ksp%NRES = NRES
    allocate(config%ksp%max_iter(NRES))
    allocate(config%ksp%min_iter(NRES))
    allocate(config%ksp%atol(NRES))
    allocate(config%ksp%rtol(NRES))
  
    config%ksp%max_iter(:) = (/NS_GMRES_itermax, NS_GMRES_itermax, LSC_GMRES_itermax, TEM_GMRES_itermax/)
    config%ksp%min_iter(:) = (/NS_GMRES_itermin, NS_GMRES_itermin, LSC_GMRES_itermin, TEM_GMRES_itermin/)
    config%ksp%atol(:) = (/NS_GMRES_atol, NS_GMRES_atol, LSC_GMRES_atol, TEM_GMRES_atol/)
    config%ksp%rtol(:) = (/NS_GMRES_tol, NS_GMRES_tol, LSC_GMRES_tol, TEM_GMRES_tol/)

    config%newton_raphson%NRES = NRES
    allocate(config%newton_raphson%atol(NRES))
    allocate(config%newton_raphson%rtol(NRES))

    config%newton_raphson%max_iter = NS_NL_itermax
    config%newton_raphson%min_iter = 1
    config%newton_raphson%atol(:) = (/NS_NL_Uatol, NS_NL_Patol, LSC_NL_atol, LSC_NL_atol/)
    config%newton_raphson%rtol(:) = (/NS_NL_Utol, NS_NL_Ptol, LSC_NL_tol, LSC_NL_tol/)

    call MPI_Comm_size(MPI_COMM_WORLD, config%mpi%numnodes, mpi_err)
    call MPI_Comm_rank(MPI_COMM_WORLD, config%mpi%myid, mpi_err)
    config%mpi%ismaster = config%mpi%myid == 0
  end subroutine init_config

  subroutine finalize_config(config)
    use aAdjKeep
    use commonvars
  
    implicit none
    type(ConfigType), intent(inout) :: config

    deallocate(config%bc%BCugType)
    deallocate(config%bc%BCugValu)
    deallocate(config%bc%BCphigType)
    deallocate(config%bc%BCphigValu)
    deallocate(config%bc%BCTgType)
    deallocate(config%bc%BCTgValu)

    deallocate(config%ksp%max_iter)
    deallocate(config%ksp%min_iter)
    deallocate(config%ksp%atol)
    deallocate(config%ksp%rtol)
    deallocate(config%newton_raphson%atol)
    deallocate(config%newton_raphson%rtol)

  end subroutine finalize_config
end module configuration
