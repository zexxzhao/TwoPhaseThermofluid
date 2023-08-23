!======================================================================
!
!======================================================================
subroutine assembleQuenching(assemble_tensor_flag, assemble_field_flag, &
                             config, mesh, sp, bc, field, &
                             vec, mat)
  ! use aAdjKeep
  use mpi
  ! use commonvars
  use commonpars
  use class_def
  use configuration

  implicit none
  ! include 'mpif.h'

  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: assemble_field_flag ! assemble NS + LS/VOF + Tem

  type(ConfigType), intent(in) :: config
  type(MeshData), intent(in) :: mesh
  type(SparsityPattern), intent(in) :: sp
  type(DirichletBCData), intent(in) :: bc
  type(FieldData), intent(in) :: field

  type(RHSData), intent(inout) :: vec
  type(LHSData), intent(inout) :: mat


  real(8) :: dgAlpha(mesh%NNODE, mesh%NSD), &
             ugAlpha(mesh%NNODE, mesh%NSD), ugmAlpha(mesh%NNODE, mesh%NSD), &
             acgAlpha(mesh%NNODE, mesh%NSD), acgmAlpha(mesh%NNODE, mesh%NSD), &
             pgAlpha(mesh%NNODE), phigAlpha(mesh%NNODE), rphigAlpha(mesh%NNODE), &
             rTgAlpha(mesh%NNODE), TgAlpha(mesh%NNODE)

  real(8) :: t1, t2
  real(8) :: gami, Delt, alfi, almi, beti, rhoinf

  Delt = config%time_integral%delt
  rhoinf = config%time_integral%rhoinf

  almi = (3d0 - rhoinf) / (1d0 + rhoinf) * 0.5d0
  alfi = 1.0d0 / (1d0 + rhoinf)

  gami = 0.5d0 + almi - alfi
  beti = 0.25d0 * (1d0 + almi - alfi) * (1d0 + almi - alfi)
 
  !---------------------------------
  ! Alpha stage
  !---------------------------------
  acgAlpha = field%acgold + almi*(field%acg - field%acgold)
  acgmAlpha = field%acgmold + almi*(field%acgm - field%acgmold)
  ugAlpha = field%ugold + alfi*(field%ug - field%ugold)
  ugmAlpha = field%ugmold + alfi*(field%ugm - field%ugmold)
  dgAlpha = field%dgold + alfi*(field%dg - field%dgold)
  pgAlpha = field%pg

  phigAlpha = field%phigold + alfi*(field%phig - field%phigold)
  rphigAlpha = field%rphigold + almi*(field%rphig - field%rphigold)

  TgAlpha = field%Tgold + alfi*(field%Tg - field%Tgold)
  rTgAlpha = field%rTgold + almi*(field%rTg - field%rTgold)

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    vec%RHSGu = 0.0d0
    vec%RHSGp = 0.0d0
    vec%RHSGls = 0.0d0
    vec%RHSGtem = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    mat%LHSK11 = 0.0d0
    mat%LHSG = 0.0d0
    mat%LHSD1 = 0.0d0
    mat%LHSM = 0.0d0

    mat%LHSLS = 0.0d0
    mat%LHSULS = 0.0d0
    mat%LHSLSU = 0.0d0
    mat%LHSPLS = 0.0d0

    mat%LHSTem = 0.0d0
  endif

  if (ismaster) then
    call CPU_TIME(t1)
  endif
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF) > 0) then
    if(config%vms%use_vms) then
    else
      call IntElmAss_NSVOF_Quenching_STAB( &
        config, mesh, sp, bc, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag, &
        vec, mat)
    endif
    call FaceAssembly_NS_weak_CF(&
      config, mesh, sp, bc, &
      dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
      acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
      assemble_tensor_flag, &
      vec, mat)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    if(config%vms%use_vms) then
    else
      call IntElmAss_Tem_Quenching_STAB( &
        config, mesh, sp, bc, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag, &
        vec, mat)
    endif
  end if

  if (ismaster) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
      call commu(vec%RHSGp, mesh%NNODE, 1, 'in ')
      call commu(vec%RHSGu, mesh%NNODE, mesh%NSD, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0) then
      call commu(vec%RHSGls, mesh%NNODE, 1, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
      call commu(vec%RHSGTem, mesh%NNODE, 1, 'in ')
    endif
  end if

end subroutine assembleQuenching
!======================================================================
!
!======================================================================

subroutine IntElmAss_NSVOF_Quenching_STAB(&
  config, mesh, sp, bc, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag, &
  vec, mat)

  ! use aAdjKeep
  ! use commonvars
  use commonpars
  use mpi
  use class_def
  use configuration
  implicit none
  ! include 'mpif.h'

  type(ConfigType), intent(in) :: config
  type(MeshData), intent(in) :: mesh
  type(SparsityPattern), intent(in) :: sp
  type(DirichletBCData), intent(in) :: bc


  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(mesh%NNODE, mesh%NSD), ugAlpha(mesh%NNODE, mesh%NSD), &
                         acgAlpha(mesh%NNODE, mesh%NSD), ugmAlpha(mesh%NNODE, mesh%NSD), &
                         acgmAlpha(mesh%NNODE, mesh%NSD), pgAlpha(mesh%NNODE), &
                         phigAlpha(mesh%NNODE), rphigAlpha(mesh%NNODE), &
                         TgAlpha(mesh%NNODE), rTgAlpha(mesh%NNODE)
                         ! Local variables

  type(RHSData), intent(inout) :: vec
  type(LHSData), intent(inout) :: mat

  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(mesh%NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, k_dc, k_dc_phi, tauP, tauLS

  !real(8), allocatable :: xKebe11(:, :, :), xGebe(:, :, :), xDebe1(:, :, :), &
  !                        xMebe(:, :), Rhsu(:, :), Rhsm(:, :), Rhsp(:), Rhsq(:), Rhsl(:), &
  !                        xLSebe(:, :), xLSUebe(:, :, :), xPLSebe(:, :), &
  !                        xULSebe(:, :, :), Rhsphi(:)

  real(8), allocatable :: xKebe11(:, :, :), xGebe(:, :, :)
  real(8), allocatable :: xDebe1(:, :, :), xMebe(:, :)
  real(8), allocatable :: xLSebe(:, :), xLSUebe(:, :, :)
  real(8), allocatable :: xULSebe(:, :, :), xPLSebe(:, :)

  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsphi(:)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: xl(:, :)
  real(8), allocatable :: dl(:, :), wl(:)
  real(8), allocatable :: acml(:, :), uml(:, :)
  real(8), allocatable :: fl(:, :), ul(:, :), acl(:, :), pl(:)
  real(8), allocatable :: phil(:), rphil(:)
  real(8), allocatable :: Tl(:), rTl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:)

  ! integer, allocatable :: ibc_loc(:, :)

  real(8) :: di(mesh%NSD), ui(mesh%NSD), aci(mesh%NSD), umi(mesh%NSD), acmi(mesh%NSD), rphii, &
             fi(mesh%NSD), xi(mesh%NSD), ddidxi(mesh%NSD, mesh%NSD), phii, &
             dphidxi(mesh%NSD), duidxi(mesh%NSD, mesh%NSD), dphidxidxj(mesh%NSD, mesh%NSD), &
             duidxixj(mesh%NSD, mesh%NSD, mesh%NSD), duiolddxi(mesh%NSD, mesh%NSD), &
             dxidx(mesh%NSD, mesh%NSD), dpridxi(mesh%NSD), rLi(mesh%NSD)
  real(8) :: uadvi(mesh%NSD), uadvi_full(mesh%NSD)
  real(8) :: Ti, rTi, dTdxi(mesh%NSD)
  real(8) :: ns_kdc

  real(8) :: Gij(mesh%NSD, mesh%NSD), Ginv(mesh%NSD, mesh%NSD)
  real(8) :: cfl, gcfl, cfl_loc
  real(8) :: Ftens(mesh%NSD, mesh%NSD), Stens(mesh%NSD, mesh%NSD), Ctens(mesh%NSD, mesh%NSD, mesh%NSD, mesh%NSD)
  real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)
  logical :: is_fluid

  real(8) :: rhoi, mui, cpi, hki
  real(8) :: mdot, vdot, phic, dmdphii, drhodphi, dmudphi
  real(8) :: divu, lambda, rm(mesh%NSD), tmp1(mesh%NSD)
  real(8) :: fact1, fact2

  integer :: NNODE, NSD
  real(8) :: gami, Delt, alfi, almi, beti, rhoinf

  real(8) :: c_evap, c_cond, Ts, lh
  real(8) :: rhoa, rhow, rhos
  real(8) :: muw, mua, mus
  real(8) :: c_dmdot

  real(8) :: DetJ, DetJw
  ! integer :: mpi_err

  NNODE = mesh%NNODE
  NSD = mesh%NSD

  Delt = config%time_integral%delt
  rhoinf = config%time_integral%rhoinf

  almi = (3d0 - rhoinf) / (1d0 + rhoinf) * 0.5d0
  alfi = 1.0d0 / (1d0 + rhoinf)

  gami = 0.5d0 + almi - alfi
  beti = 0.25d0 * (1d0 + almi - alfi) * (1d0 + almi - alfi)
 
  !---------------------------------
  c_evap = config%property%c_evap
  c_cond = config%property%c_cond
  Ts = config%property%Ts
  lh = config%property%lh

  rhoa = config%property%rhoa
  rhow = config%property%rhow
  rhos = config%property%rhos

  mua = config%property%mua
  muw = config%property%muw
  mus = config%property%mus

  c_dmdot = config%vms%c_dmdot

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  fact1 = almi
  fact2 = alfi * delt * gami


  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

  ! loop over elements
  do iel = 1, mesh%NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= mesh%ELMNSHL(iel)) then

      if (NSHL >= 0) then
        deallocate (xKebe11, xGebe)                      ! 2
        deallocate (xDebe1, xMebe)                       ! 4
        deallocate (xLSebe, xLSUebe)                     ! 6
        deallocate (xULSebe, xPLSebe)                    ! 8
        deallocate (Rhsu, Rhsp, Rhsphi)                  ! 11
        deallocate (shlu, shgradgu, shhessgu)            ! 14
        deallocate (xl)                                  ! 15
        deallocate (dl, wl)                              ! 17
        deallocate (acml, uml)                           ! 19
        deallocate (fl, ul, acl, pl)                     ! 23
        deallocate (phil, rphil)                         ! 25
        deallocate (Tl, rTl)                             ! 27
        deallocate (gp, gw)                              ! 29
        deallocate (shconv, shconv_full)                 ! 31
      end if

      NSHL = mesh%ELMNSHL(iel)
      NGAUSS = mesh%ELMNGAUSS(iel)
      allocate (xKebe11(NSD * NSD, NSHL, NSHL), xGebe(NSD, NSHL, NSHL))
      allocate (xDebe1(NSD, NSHL, NSHL), xMebe(NSHL, NSHL))
      allocate (xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL))
      allocate (xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL))
      allocate (Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsphi(NSHL))
      allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD))
      allocate (xl(NSHL, NSD))
      allocate (dl(NSHL, NSD), wl(NSHL))
      allocate (acml(NSHL, NSD), uml(NSHL, NSD))
      allocate (fl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), pl(NSHL))
      allocate (phil(NSHL), rphil(NSHL))
      allocate (Tl(NSHL), rTl(NSHL))
      allocate (gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (shconv(NSHL), shconv_full(NSHL))
      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if

    is_fluid = mesh%ELM_ID(iel) == 102

    fl = 0.0d0
!    do i = 1, NSHL
!      fl(i,3) = -1.0d0/(Fr**2.0d0)
!    end do

    ! Get local solution arrays
    do i = 1, NSHL
      idx = mesh%IEN(iel, i)
      xl(i, :) = mesh%xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      ul(i, :) = ugAlpha(idx, :)
      acl(i, :) = acgAlpha(idx, :)
      uml(i, :) = ugmAlpha(idx, :)
      acml(i, :) = acgmAlpha(idx, :)
      pl(i) = pgAlpha(idx)
      phil(i) = phigAlpha(idx)
      rphil(i) = rphigAlpha(idx)
      Tl(i) = TgAlpha(idx)
      rTl(i) = rTgAlpha(idx)
    end do


    ! initialize local stiffness matrix
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      xLSebe = 0.0d0
      xLSUebe = 0.0d0
      xULSebe = 0.0d0
      xPLSebe = 0.0d0

    end if
    ! initialize local load vector
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      Rhsu = 0.0d0
      Rhsp = 0.0d0
      Rhsphi = 0.0d0
    end if
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      ! initialize
      shlu(:) = 0.0d0
      shgradgu(:, :) = 0.0d0
      shhessgu(:, :, :) = 0.0d0
      ! hess_flag = .false.
      call eval_shape(NSD, nshl, iel, gp(igauss, :), xl, dl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, DetJ, 0)

      DetJw = DetJ * gw(igauss)
      !call e3int_fluid(nshl, xl, dl, ul, acl, uml, acml, &
      !                 pl, fl, phil, Tl, rTl, shlu, shgradgu, &
      !                 shgradgu, shhessgu, dxidx, Ginv, &
      !                 di, ui, aci, umi, acmi, pri, fi, &
      !                 ddidxi, duidxi, duidxixj, dpridxi, &
      !                 phi, dphidxi, dphidxidxj, rphil, dphidti, xi, &
      !                 rTi, Ti, dTdxi, rhoi, mui, cpi, hki, &
      !                 rLi, res_phic_tmp1)
      call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi)
      call e3int_qr(NSHL, NSD, NSD, shlu, dl, di)
      call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui)
      call e3int_qr(NSHL, NSD, NSD, shlu, acl, aci)
      call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi)
      call e3int_qr(NSHL, NSD, NSD, shlu, acml, acmi)
      call e3int_qr(NSHL, NSD, NSD, shlu, fl, fi)
      call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
      call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
      call e3int_qr(NSHL, NSD, 1, shlu, rphil, rphii)
      call e3int_qr(NSHL, NSD, 1, shlu, Tl, Ti)
      call e3int_qr(NSHL, NSD, 1, shlu, rTl, rTi)

      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, dl, ddidxi)
      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, pl, dpridxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, phil, dphidxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, Tl, dTdxi)

      call e3int_qr_hess(NSHL, NSD, NSD, shhessgu, ul, duidxixj)
      call e3int_qr_hess(NSHL, NSD, 1, shhessgu, phil, dphidxidxj)

        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate on qr", igauss, NGAUSS
      phic = max(min(phii, 1.0d0), 0.0d0)
      ! phic = phii
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdphii = -c_evap * rhow * (Ti - Ts) / Ts * c_dmdot
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdphii = c_cond * rhoa * (Ti - Ts) / Ts * c_dmdot
      endif
      vdot = mdot / rhoa - mdot / rhow

      if(is_fluid) then
        call prop_interp(rhow, rhoa, phii, rhoi)
        call prop_interp(muw, mua, phii, mui)
      else ! solid
        rhoi = rhos
        mui = mus
        if(.not. config%vms%use_sliding_velocity) then
          ui(:) = umi(:)
        endif
      endif
      lambda = -2d0 / 3d0 * mui
      drhodphi = rhoa - rhow
      dmudphi = mua - muw
      if(phii < 0d0 .or. phii > 1d0) then
        drhodphi = 0d0
        dmudphi = 0d0
        dmdphii = 0d0
      endif

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      divu = duidxi(1, 1) + duidxi(2, 2) + duidxi(3, 3)

      fi(:) = gravvec(:) * rhoi
      rm(1) = rhoi * (aci(1) + sum(duidxi(1, :) * uadvi(:))) - fi(1)
      rm(2) = rhoi * (aci(2) + sum(duidxi(2, :) * uadvi(:))) - fi(2)
      rm(3) = rhoi * (aci(3) + sum(duidxi(3, :) * uadvi(:))) - fi(3)

      ! call a function to get residual for NS and VOF
      ! call e3int_rLi(NSD, rhoi, mui, aci, duidxi, uadvi, dpridxi, fi, duidxixj, rLi)
      ! call e3int_resphi(NSD, rphii, phii, dphidxi, uadvi, mdot, divu, rhoa, res_phic_tmp1)
      rLi(:) = rhoi * aci(:) + &
        rhoi * (duidxi(:, 1) * uadvi(1) + duidxi(:, 2) * uadvi(2) + duidxi(:, 3) * uadvi(3)) + &
        dpridxi(:) - fi(:) !- &
        !mui * (duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) - &
        !(mui+lambda) * (duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))

      res_phic_tmp1 = rphii + sum(uadvi(:) * dphidxi(:)) + phii * vdot - mdot / rhoa
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate on RESIDUAL", igauss, NGAUSS
      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:) * shgradgu(aa, :))
      enddo

      tauBar = 0.0d0
      if(config%vms%use_taubar) then
        call e3STAB_3D_NSVOF_TAUBAR(NSD, Gij, uadvi, uprime, tauBar)
      endif

      k_dc = 0.0d0
      if(abs(config%vms%NS_kdc_w) + abs(config%vms%NS_kdc_a)  > 0.0d0) then
        call prop_interp(config%vms% NS_kdc_w, &
                         config%vms%NS_kdc_a, phii, ns_kdc)
        call e3DC_beta2(NSD, ns_kdc, Gij, rLi, k_dc)
      endif

      k_dc_phi = 0.0d0
      if(abs(config%vms%LSC_kdc) > 0.0d0) then
        call e3DC_scalar(NSD, config%vms%LSC_kdc, Gij, res_phic_tmp1, dphidxi, k_dc_phi)
      endif
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate KDC", igauss, NGAUSS

      if(config%calc_cfl) then
        call e3CFL(NSD, uadvi, Gij, Delt, cfl_loc)
        cfl = max(cfl, cfl_loc)
      endif
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate CFL", igauss, NGAUSS

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0 .and. is_fluid) then
        do aa = 1, NSHL
          do bb = 1, NSHL
            do i = 1, NSD
              do j = 1, NSD
                xKebe11((i-1)*NSD+j, aa, bb) = xKebe11((i-1)*NSD+j, aa, bb) + &
                  fact2 * shgradgu(aa, j) * shgradgu(bb, i) * (mui + k_dc) * DetJw + &
                  fact2 * shgradgu(aa, i) * shgradgu(bb, j) * (tauC+lambda) * DetJw
              enddo
              xKebe11((i-1)*NSD+i, aa, bb) = xKebe11((i-1)*NSD+i, aa, bb) + &
                fact1 * shlu(aa) * rhoi * shlu(bb) * DetJw + &
                fact2 * shlu(aa) * rhoi * shconv(bb) * DetJw + &
                fact2 * (mui + k_dc) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJw + &
                fact1 * rhoi * shconv(aa) * tauM * rhoi * shlu(bb) * DetJw + &
                fact2 * rhoi * shconv(aa) * tauM * rhoi * shconv(bb) * DetJw
            enddo
            xGebe(:, aa, bb) = xGebe(:, aa, bb) - &
              fact2 * shgradgu(aa, :) * shlu(bb) * DetJw + &
              fact2 * rhoi * shconv(aa) * tauM * shgradgu(bb, :) * DetJw

            xDebe1(:, aa, bb) = xDebe1(:, aa, bb) + &
              fact2 * shlu(aa) * shgradgu(bb, :) * DetJw + &
              fact1 * shgradgu(aa, :) * tauM * rhoi * shlu(bb) * DetJw + &
              fact2 * shgradgu(aa, :) * tauM * rhoi * shconv(bb) * DetJw

            xMebe(aa, bb) = xMebe(aa, bb) + &
              fact2 * tauM * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJw

            xLSebe(aa, bb) = xLSebe(aa, bb)  + &
              (shlu(aa) + tauls * shconv(aa)) * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJw
            xLSebe(aa, bb) = xLSebe(aa, bb)  + &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * shlu(bb) * vdot * DetJw + &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * phii * (1/rhoa - 1/rhow) * dmdphii * shlu(bb) * DetJw - &
              (shlu(aa) + tauls * shconv(aa)) * fact2 * dmdphii * shlu(bb) / rhoa * DetJw
            xLSebe(aa, bb) = xLSebe(aa, bb) + &
              fact2 * k_dc_phi * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJw

            ! if(.false.) then
            xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
              fact2 * shlu(aa) * drhodphi * rm(:) / rhoi * shlu(bb) * DetJw + &
              fact2 * rhoi * shconv(aa) * tauM * drhodphi * rm(:) / rhoi * shlu(bb) * DetJw + &
              fact2 * shgradgu(aa, :) * divu * (-2d0/3d0) * dmudphi * shlu(bb) * DetJw + &
              fact2 * shgradgu(aa, :) * divu * drhodphi * tauC / rhoi * shlu(bb) * DetJw

            xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(1, :) + duidxi(:, 1))) * dmudphi * shlu(bb) * DetJw
            xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(2, :) + duidxi(:, 2))) * dmudphi * shlu(bb) * DetJw
            xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(3, :) + duidxi(:, 3))) * dmudphi * shlu(bb) * DetJw

            xPLSebe(aa, bb) = xPLSebe(aa, bb) - &
              fact2 * shlu(aa) * dmdphii * shlu(bb) * (1/rhoa - 1/rhow) * DetJw !- &
              !fact2 * sum(shgradgu(aa, :) * dpridxi(:)) * tauM / rhoi * drhodphi * shlu(bb) * DetJw
            
            xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) &
              + fact2 * (shlu(aa) + tauls * shconv(aa)) * shlu(bb) * dphidxi(:) * DetJw !&
              ! + fact2 * tauls * shgradgu(aa, :) * shlu(bb) * res_phic_tmp1 * DetJw !&
              ! + fact2 * (shlu(aa) + tauls * shconv(aa)) * phii * shgradgu(bb, :) * DetJw
            ! endif
          end do
        end do

        ! call e3LHS_3D_fluid_quenching(&
        !   nshl, ui, umi, aci, pri, duidxi, dpridxi, &
        !   dphidxi, dphidxidxj, rphii, phii, &
        !   rLi, fi, rhoi, mui, Ti, &
        !   tauM, tauP, tauLS, tauC, tauBar, tauBar1, &
        !   k_dc, k_dc_phi, gw(igauss), shlu, shgradgu, &
        !   shhessgu, xKebe11, xGebe, xDebe1, xMebe, &
        !   xLSebe, xLSUebe, xULSebe, xPLSebe)
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_fluid) then
        tmp1(:) = rhoi *aci(:) + &
          rhoi * uadvi(1) * duidxi(:, 1) + &
          rhoi * uadvi(2) * duidxi(:, 2) + &
          rhoi * uadvi(3) * duidxi(:, 3) - fi(:)
        do aa = 1, NSHL
          RHSu(:, aa) = RHSu(:, aa) - &
            shlu(aa) * tmp1(:) * DetJw - &
            shgradgu(aa, :) * (-pri + (tauC + lambda) * divu) * DetJw - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 1) + duidxi(1, :))) * DetJw - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 2) + duidxi(2, :))) * DetJw - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 3) + duidxi(3, :))) * DetJw + &
            rhoi * shconv(aa) * uprime(:) * DetJw
          RHSp(aa) = RHSp(aa) - &
            shlu(aa) * (divu - vdot) * DetJw + &
            sum(shgradgu(aa, :) * uprime(:)) * DetJw
          RHSphi(aa) = RHSphi(aa) - &
            (shlu(aa) + tauls * shconv(aa)) * res_phic_tmp1 * DetJw - &
            k_dc_phi * sum(shgradgu(aa, :) * dphidxi(:)) * DetJw
        enddo
      end if

    end do

    ! Apply Dirichlet BCs
    call apply_bc(&
      mesh%nsd, nshl, mesh%ien(iel, :), bc, &
      xKebe11, xGebe, xDebe1, xMebe, &
      Rhsu, Rhsp, &
      xLSebe, xLSUebe, xULSebe, xPLSebe, &
      Rhsphi)
    ! call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
    !               xMebe, Rhsu, Rhsp, &
    !               xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
    !    call MPI_barrier(MPI_COMM_WORLD, mpi_err)
    !if(ismaster) write(*,*) "Evaluate BC", NSHL
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(mesh%NNODE, mesh%NSD, &
                            sp%nnz, sp%indices, sp%indptr, &
                            nshl, mesh%ien(iel, :), &
                            xKebe11, xGebe, xDebe1, xMebe, &
                            xLSebe, xLSUebe, xULSebe, xPLSebe, &
                            mat%LHSK11, mat%LHSG, mat%LHSD1, mat%LHSM, &
                            mat%LHSLS, mat%LHSULS, mat%LHSLSU, mat%LHSPLS)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      call LocalToGlobalNSVOF_3D(vec%RHSGu, vec%RHSGp, vec%RHSGls, &
                                 mesh%NNODE, mesh%NSD, mesh%ELMNSHL(iel), mesh%maxNSHL, &
                                 mesh%IEN(iel, :), rhsu, rhsp, rhsphi)
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate add_local"
    end if

  end do

  deallocate (xKebe11, xGebe)                      ! 2
  deallocate (xDebe1, xMebe)                       ! 4
  deallocate (xLSebe, xLSUebe)                     ! 6
  deallocate (xULSebe, xPLSebe)                    ! 8
  deallocate (Rhsu, Rhsp, Rhsphi)                  ! 11
  deallocate (shlu, shgradgu, shhessgu)            ! 14
  deallocate (xl)                                  ! 15
  deallocate (dl, wl)                              ! 17
  deallocate (acml, uml)                           ! 19
  deallocate (fl, ul, acl, pl)                     ! 23
  deallocate (phil, rphil)                         ! 25
  deallocate (Tl, rTl)                             ! 27
  deallocate (gp, gw)                              ! 29
  deallocate (shconv, shconv_full)                 ! 31
  ! write(*,*) "DEBUG, myid=", myid, config%calc_cfl
  if (config%calc_cfl) then
    ! Find largest CFL-number and output to screen
    if (numnodes > 1) then
      call MPI_ALLReduce(cfl, gcfl, 1, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD, mpi_err)
    else
      gcfl = cfl
    end if
    if (ismaster) then
      write (*, '(40("-"))')
      write (*, *) "    CFL = ", gcfl
    end if
  end if

end subroutine IntElmAss_NSVOF_Quenching_STAB


!======================================================================
!
!======================================================================
subroutine IntElmAss_Tem_Quenching_STAB(&
  config, mesh, sp, bc, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag, &
  vec, mat)
  ! use aAdjKeep
  !use commonvars
  use class_def
  use commonpars
  ! use mpi
  use configuration

  implicit none
  include 'mpif.h'

  type(ConfigType), intent(in) :: config
  type(MeshData), intent(in) :: mesh
  type(SparsityPattern), intent(in) :: sp
  type(DirichletBCData), intent(in) :: bc


  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(mesh%NNODE, mesh%NSD), ugAlpha(mesh%NNODE, mesh%NSD), &
                         acgAlpha(mesh%NNODE, mesh%NSD), ugmAlpha(mesh%NNODE, mesh%NSD), &
                         acgmAlpha(mesh%NNODE, mesh%NSD), pgAlpha(mesh%NNODE), &
                         phigAlpha(mesh%NNODE), dphidtgAlpha(mesh%NNODE), &
                         TgAlpha(mesh%NNODE), rTgAlpha(mesh%NNODE)
  
  type(RHSData), intent(inout) :: vec
  type(LHSData), intent(inout) :: mat
  ! Local variables
  real(8), parameter :: damp = 0.5d0

  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(mesh%NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, kappa_str, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: RHStem(:), xTebe(:, :)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)

  real(8), allocatable :: rTl(:), Tl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:), tmp(:), drdTi(:)

  real(8) :: di(mesh%NSD), ui(mesh%NSD), aci(mesh%NSD), umi(mesh%NSD), acmi(mesh%NSD), rphii, &
             fi(mesh%NSD), uadvi(mesh%NSD), xi(mesh%NSD), ddidxi(mesh%NSD, mesh%NSD), phii, &
             dphidxi(mesh%NSD), duidxi(mesh%NSD, mesh%NSD), dphidxidxj(mesh%NSD, mesh%NSD), &
             duidxixj(mesh%NSD, mesh%NSD, mesh%NSD), duiolddxi(mesh%NSD, mesh%NSD), &
             dxidx(mesh%NSD, mesh%NSD), dpridxi(mesh%NSD), rLi(mesh%NSD)
  real(8) :: uadvi_full(mesh%NSD)
  real(8) :: Gij(mesh%NSD, mesh%NSD), Ginv(mesh%NSD, mesh%NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(mesh%NSD, mesh%NSD), Stens(mesh%NSD, mesh%NSD), Ctens(mesh%NSD, mesh%NSD, mesh%NSD, mesh%NSD)
  !real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)

  !real(8) :: phi

  real(8) :: rTi, Ti, dTdxi(mesh%NSD), dTdxixj(mesh%NSD, mesh%NSD)
  real(8) :: res_tem1, res_tem2, tau_tem
  logical :: is_fluid
  real(8) :: fact1, fact2
  real(8) :: Se, mdot, vdot
  real(8) :: rhoi, mui, cpi, hki, rhocpi
  real(8) :: k_dc_tem
  real(8) :: dmdTi, phic, tmp1, dSedTi

  integer :: NNODE, NSD
  real(8) :: gami, Delt, alfi, almi, beti, rhoinf

  real(8) :: c_evap, c_cond, Ts, lh
  real(8) :: rhoa, rhow, rhos
  real(8) :: mua, muw, mus
  real(8) :: cpa, cpw, cps
  real(8) :: kappaa, kappaw, kappas
  real(8) :: c_dmdot

  real(8) :: DetJ, DetJw

  NNODE = mesh%NNODE
  NSD = mesh%NSD

  Delt = config%time_integral%delt
  rhoinf = config%time_integral%rhoinf

  almi = (3d0 - rhoinf) / (1d0 + rhoinf) * 0.5d0
  alfi = 1.0d0 / (1d0 + rhoinf)

  gami = 0.5d0 + almi - alfi
  beti = 0.25d0 * (1d0 + almi - alfi) * (1d0 + almi - alfi)
 
  !---------------------------------
  c_evap = config%property%c_evap
  c_cond = config%property%c_cond
  Ts = config%property%Ts
  lh = config%property%lh

  rhoa = config%property%rhoa
  rhow = config%property%rhow
  rhos = config%property%rhos

  mua = config%property%mua
  muw = config%property%muw
  mus = config%property%mus

  cpa = config%property%cpa
  cpw = config%property%cpw
  cps = config%property%cps

  kappaa = config%property%kappaa
  kappaw = config%property%kappaw
  kappas = config%property%kappas
  c_dmdot = config%vms%c_dmdot


  fact1 = almi
  fact2 = alfi * gami * Delt

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

  ! loop over elements
  do iel = 1, mesh%NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= mesh%ELMNSHL(iel)) then

      if (NSHL >= 0) then
        deallocate (RHStem, xTebe)
        deallocate (shlu, shgradgu, shhessgu)
        deallocate (fl, dl, ul, wl, acl)
        deallocate (uml, acml, pl, dlold)
        deallocate (xl, dumb, phil, ulold, dphidtl)
        deallocate (rTl, Tl)
        deallocate (gp, gw)
        deallocate (shconv, shconv_full, tmp, drdTi)
      end if

      NSHL = mesh%ELMNSHL(iel)
      NGAUSS = mesh%ELMNGAUSS(iel)
      allocate (RHStem(NSHL), xTebe(NSHL, NSHL))
      allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD))
      allocate (fl(NSHL, NSD), dl(NSHL, NSD), ul(NSHL, NSD), wl(NSHL), acl(NSHL, NSD))
      allocate (uml(NSHL, NSD), acml(NSHL, NSD), pl(NSHL), dlold(NSHL, NSD))
      allocate (xl(NSHL, NSD), dumb(NSHL, NSD), phil(NSHL), ulold(NSHL, NSD), dphidtl(NSHL))
      allocate (rTl(NSHL), Tl(NSHL))
      allocate (gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (shconv(NSHL), shconv_full(NSHL), tmp(NSHL), drdTi(NSHL))

      ! get Gaussian points and weights
      call genGPandGW(gp, gw, NGAUSS)
    end if
    is_fluid = mesh%ELM_ID(iel) == 102
    fl = 0.0d0

    ! Get local solution arrays
    do i = 1, NSHL
      idx = mesh%IEN(iel, i)
      xl(i, :) = mesh%xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      ! dlold(i, :) = field%dgold(idx, :)
      ! wl(i) = wg(idx)
      ul(i, :) = ugAlpha(idx, :)
      ! ulold(i, :) = field%ugold(idx, :)
      acl(i, :) = acgAlpha(idx, :)
      uml(i, :) = ugmAlpha(idx, :)
      acml(i, :) = acgmAlpha(idx, :)
      pl(i) = pgAlpha(idx)
      phil(i) = phigAlpha(idx)
      dphidtl(i) = dphidtgAlpha(idx)
      Tl(i) = TgAlpha(idx)
      rTl(i) = rTgAlpha(idx)
    end do

    ! initialize local stiffness matrix
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      xTebe(:, :) = 0.0d0
    end if
    ! initialize local load vector
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      RHStem(:) = 0.0d0
    end if
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      ! initialize
      shlu = 0.0d0
      shgradgu = 0.0d0
      shhessgu = 0.0d0
  
      call eval_shape(NSD, nshl, iel, gp(igauss, :), xl, dl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, DetJ, 0)

      DetJw = DetJ * gw(igauss)
      ! Fluid
      call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
      call e3int_qr(NSHL, NSD, 1, shlu, rTl, rTi)
      call e3int_qr(NSHL, NSD, 1, shlu, Tl, Ti)
      call e3int_qr(NSHL, NSD, NSD, shlu, acl, aci)
      call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui)
      call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi)

      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, pl, dpridxi)
      call e3int_qr_grad(NSHL, NSD, 1, shgradgu, Tl, dTdxi)
      call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

      call e3int_qr_hess(NSHL, NSD, 1, shhessgu, Tl, dTdxixj)
      call e3int_qr_hess(NSHL, NSD, NSD, shhessgu, ul, duidxixj)

      phic = max(min(phii, 1.0d0), 0.0d0)
      ! phic = phii
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdTi = c_evap * (1-phic) * rhow / Ts * c_dmdot
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdTi = c_cond * (phic) * rhoa / Ts * c_dmdot
      endif
      vdot = mdot / rhoa - mdot / rhow

      if(is_fluid) then
        call prop_interp(rhow, rhoa, phii, rhoi)
        call prop_interp(muw, mua, phii, mui)
        call prop_interp(cpw, cpa, phii, cpi)
        call prop_interp(kappaw, kappaa, phii, hki)
        call prop_interp(rhow * cpw, rhoa * cpa, phii, rhocpi)
      else ! solid
        rhoi = rhos
        mui = mus
        cpi = cps
        hki = kappas
        rhocpi = rhos * cps
        if(.not. config%vms%use_sliding_velocity) then
          ui(:) = umi(:)
        endif
      endif

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      ! call a function to get residual for Tem
      ! call e3int_rLi(NSD, rhoi, mui, aci, duidxi, uadvi, dpridxi, fi, duidxixj, rLi)


      ! call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      ! uprime(:) = -tauM * rLi(:)
  
      call e3int_restem(NSD, rTi, Ti, dTdxi, dTdxixj, uadvi, Se, rhocpi, hki, res_tem1)
      call e3STAB_3D_TEM(NSD, Gij, uadvi, Delt, rhoi, cpi, hki, tau_tem)

      k_dc_tem = 0.0d0
      if(abs(config%vms%Tem_kdc) > 0d0) then
        call e3DC_scalar(NSD, config%vms%Tem_kdc, Gij, res_tem1, dTdxi, k_dc_tem)
      endif

      ! uadvi_full(:) = uadvi(:) + uprime(:)
      ! do aa = 1, NSHL
      !   shconv(aa) = sum(uadvi(:)*shgradgu(aa, :))
      !   ! shconv_full(aa) = sum(uadvi_full(:)*shgradgu(aa, :))
      ! enddo
      shconv(:) = matmul(shgradgu(:, :), uadvi(:))
      ! tmp1 = rhocpi * vdot + (rhow*cpw-rhoa*cpa) * sum(uadvi_full(:) * dphidxi(:))
      Se = ((cpw - cpa) * (Ti - Ts) - Lh) * mdot 
      dSedTi = (cpw - cpa) * mdot + ((cpw - cpa) * (Ti - Ts) - lh)* dmdTi
      ! drdTi(:) = rhocpi * (fact1 * shlu(:) + fact2 * shconv(:)) - fact2 * dSedTi * shlu(:)
      ! shconv_full(:) = shconv(:)
      tmp(:)  = shlu(:) + tau_tem * rhocpi * shconv(:)
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        ! e3LHS_3D_tem()
        ! xTebe(:, :) = 0.0d0
        if(is_fluid) then
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) &
                + shlu(aa) * rhocpi * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJw &
                - shlu(aa) * fact2 * dSedTi * shlu(bb) * DetJw &
                + tau_tem * rhocpi * shconv(aa) * rhocpi &
                    * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJw &
                - tau_tem * rhocpi * shconv(aa) * fact2 * dSedTi * shlu(bb) * DetJw &
                + fact2 * (hki + k_dc_tem) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJw
                
            enddo
          enddo
        else
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) + &
                shlu(aa) * rhocpi * fact1 * shlu(bb) * DetJw + &
                fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJw
            enddo
          enddo
        endif
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        if(is_fluid) then
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * (rTi + dot_product(uadvi(:), dTdxi(:))) * DetJw
          RHSTem(:) = RHSTem(:) + shlu(:) * Se * DetJw
          RHSTem(:) = RHSTem(:) - rhocpi * tau_tem * shconv(:) * res_tem1 * DetJw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 1) * dTdxi(1) * DetJw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 2) * dTdxi(2) * DetJw
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 3) * dTdxi(3) * DetJw
        else
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * rTi * DetJw
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 1) * dTdxi(1) * DetJw
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 2) * dTdxi(2) * DetJw
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 3) * dTdxi(3) * DetJw
        endif
      end if
    end do

    ! Apply Dirichlet BCs
    call BCLHS_tem(nshl, mesh%ien(iel, :), bc, xTebe, RHSTem)
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      ! Assemble global LHS Matrix
      call FillSparseMat_tem(mesh%NNODE, mesh%NSD, &
                             sp%nnz, sp%indices, sp%indptr, &
                             nshl, mesh%ien(iel, :), xTebe, mat%LHStem)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! Assemble load RHS vector
      do aa = 1, NSHL
        vec%RHSGTEM(mesh%IEN(iel, aa)) = vec%RHSGTEM(mesh%IEN(iel, aa)) + Rhstem(aa)
      end do
    end if
  end do

  deallocate (RHStem, xTebe)
  deallocate (shlu, shgradgu, shhessgu)
  deallocate (fl, dl, ul, wl, acl)
  deallocate (uml, acml, pl, dlold)
  deallocate (xl, dumb, phil, ulold, dphidtl)
  deallocate (rTl, Tl)
  deallocate (gp, gw)
  deallocate (shconv, shconv_full, tmp, drdTi)


end subroutine IntElmAss_Tem_Quenching_STAB
