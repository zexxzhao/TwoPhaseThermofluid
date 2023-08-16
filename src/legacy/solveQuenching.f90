!======================================================================
!
!======================================================================
subroutine assembleQuenching(config, assemble_tensor_flag, assemble_field_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  use configuration

  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: assemble_tensor_flag ! assembe Jacobian mat or vec
  integer, intent(in) :: assemble_field_flag ! assemble NS + LS/VOF + Tem

  real(8) :: dgAlpha(NNODE, NSD), &
             ugAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
             acgAlpha(NNODE, NSD), acgmAlpha(NNODE, NSD), &
             pgAlpha(NNODE), phigAlpha(NNODE), rphigAlpha(NNODE), &
             rTgAlpha(NNODE), TgAlpha(NNODE)

  real(8) :: t1, t2

  !---------------------------------
  ! Alpha stage
  !---------------------------------
  acgAlpha = acgold + almi*(acg - acgold)
  acgmAlpha = acgmold + almi*(acgm - acgmold)
  ugAlpha = ugold + alfi*(ug - ugold)
  ugmAlpha = ugmold + alfi*(ugm - ugmold)
  dgAlpha = dgold + alfi*(dg - dgold)
  pgAlpha = pg

  phigAlpha = phigold + alfi*(phig - phigold)
  rphigAlpha = rphigold + almi*(rphig - rphigold)

  TgAlpha = Tgold + alfi*(Tg - Tgold)
  rTgAlpha = rTgold + almi*(rTg - rTgold)

  !---------------------------------
  ! zero out RHS and LHS
  !---------------------------------

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    RHSGu = 0.0d0
    RHSGp = 0.0d0
    RHSGls = 0.0d0
    RHSGtem = 0.0d0
  end if

  if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    LHSK11 = 0.0d0
    LHSG = 0.0d0
    LHSD1 = 0.0d0
    LHSM = 0.0d0

    LHSLS = 0.0d0
    LHSULS = 0.0d0
    LHSLSU = 0.0d0
    LHSPLS = 0.0d0

    LHSTem = 0.0d0
  endif

  if (myid .eq. 0) then
    call CPU_TIME(t1)
  endif
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_NS + ASSEMBLE_FIELD_VOF) > 0) then
    if(config%vms%use_vms) then
      call IntElmAss_NSVOF_Quenching_VMS( &
        config, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag)
    else
      call IntElmAss_NSVOF_Quenching_STAB( &
        config, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag)
    endif
    call FaceAssembly_NS_weak(config, &
                              dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                              acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                              assemble_tensor_flag)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    if(config%vms%use_vms) then
      call IntElmAss_Tem_Quenching_VMS( &
        config, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag)
    else
      call IntElmAss_Tem_Quenching_STAB( &
        config, &
        dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
        acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
        TgAlpha, rTgAlpha, &
        assemble_tensor_flag)
    endif
  end if

  if (myid .eq. 0) then
    call CPU_TIME(t2)
    write (*, *) "Total time assemble:", t2 - t1, "seconds"
  end if

  if (numnodes > 1 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_NS) > 0) then
      call commu(RHSGp, 1, 'in ')
      call commu(RHSGu, NSD, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_VOF) > 0) then
      call commu(RHSGls, 1, 'in ')
    endif
    if(iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
      call commu(RHSGTem, 1, 'in ')
    endif
  end if

end subroutine assembleQuenching
!======================================================================
!
!======================================================================

subroutine IntElmAss_NSVOF_Quenching_VMS(&
  config, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  use configuration
  implicit none

  type(ConfigType), intent(in) :: config

  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), rphigAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
                         ! Local variables
  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
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

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi(NSD), uadvi_full(NSD)
  real(8) :: Ti, rTi, dTdxi(NSD)
  real(8) :: ns_kdc

  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl, gcfl, cfl_loc
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)
  logical :: is_fluid

  real(8) :: rhoi, mui, cpi, hki
  real(8) :: mdot, vdot, phic, dmdphii, drhodphi, dmudphi
  real(8) :: divu, lambda, rm(NSD), tmp1(NSD)
  real(8) :: fact1, fact2

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  fact1 = almi
  fact2 = alfi * delt * gami

  drhodphi = rhow - rhoa
  dmudphi = muw - mua

  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

  ! loop over elements
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

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

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
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

    is_fluid = ELM_ID(iel) == 102

    fl = 0.0d0
!    do i = 1, NSHL
!      fl(i,3) = -1.0d0/(Fr**2.0d0)
!    end do

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      wl(i) = wg(idx)
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
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)

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
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdphii = -c_evap * rhow * (Ti - Ts) / Ts
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdphii = c_cond * rhoa * (Ti - Ts) / Ts
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
      lambda = -2d0 / 3 * mui

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
        dpridxi(:) - fi(:) - &
        mui * (duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) + &
        (mui+lambda) * (duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))

      res_phic_tmp1 = rphii + sum(uadvi(:) * dphidxi(:)) + phii * vdot - mdot / rhoa
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate on RESIDUAL", igauss, NGAUSS
      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)
      uadvi_full(:) = uadvi(:) + uprime(:)
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:) * shgradgu(aa, :))
        shconv_full(aa) = sum(uadvi_full(:) * shgradgu(aa, :))
      enddo
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate STAB", igauss, NGAUSS, &
        !     config%vms%use_taubar

      tauBar = 0.0d0
      if(config%vms%use_taubar) then
        call e3STAB_3D_NSVOF_TAUBAR(NSD, Gij, uadvi, uprime, tauBar)
      endif

      k_dc = 0.0d0
      if(abs(config%vms%NS_kdc_w) + abs(config%vms%NS_kdc_a)  > 0.0d0) then
        call prop_interp(NS_kdc_w, NS_kdc_a, phii, ns_kdc)
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

        !call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(.false.) then
        !   do aa = 1, NSHL
        !     do bb = 1, NSHL
        !       do i = 1, NSD
        !         do j = 1, NSD
        !           xKebe11((i-1)*NSD+j, aa, bb) = xKebe11((i-1)*NSD+j, aa, bb) + &
        !             fact2 * shgradgu(aa, j) * shgradgu(bb, i) * (mui + k_dc) * DetJ * gw(igauss) + &
        !             fact2 * shgradgu(aa, i) * shgradgu(bb, j) * (tauC+lambda) * DetJ * gw(igauss)
        !         enddo
        !         xKebe11((i-1)*NSD+i, aa, bb) = xKebe11((i-1)*NSD+i, aa, bb) + &
        !           fact1 * shlu(aa) * rhoi * shlu(bb) * DetJ * gw(igauss) + &
        !           fact2 * shlu(aa) * rhoi * shconv(bb) * DetJ * gw(igauss) + &
        !           fact2 * (mui + k_dc) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJ * gw(igauss) + &
        !           fact1 * rhoi * shconv(aa) * tauM * rhoi * shlu(bb) * DetJ * gw(igauss) + &
        !           fact2 * rhoi * shconv(aa) * tauM * rhoi * shconv(bb) * DetJ * gw(igauss)
        !       enddo
        !       xGebe(:, aa, bb) = xGebe(:, aa, bb) - &
        !         fact2 * shgradgu(aa, :) * shlu(bb) * DetJ * gw(igauss) + &
        !         fact2 * rhoi * shconv(aa) * tauM * shgradgu(bb, :) * DetJ * gw(igauss)
        !       xDebe1(:, aa, bb) = xDebe1(:, aa, bb) + &
        !         fact2 * shlu(aa) * shgradgu(bb, :) * DetJ * gw(igauss) + &
        !         fact1 * shgradgu(aa, :) * tauM * rhoi * shlu(bb) * DetJ * gw(igauss) + &
        !         fact2 * shgradgu(aa, :) * tauM * rhoi * shconv(bb) * DetJ * gw(igauss)
        !       xMebe(aa, bb) = xMebe(aa, bb) + &
        !         fact2 * tauM * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJ * gw(igauss)

        !       xLSebe(aa, bb) = xLSebe(aa, bb)  + &
        !         (shlu(aa) + tauls * shconv(aa)) * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJ * gw(igauss)
        !       xLSebe(aa, bb) = xLSebe(aa, bb)  + (shlu(aa) + tauls * shconv(aa)) * &
        !         fact2 * (shlu(bb) * divu - dmdphii * shlu(bb) / rhoa) * DetJ * gw(igauss)
        !         ! (fact1 * shlu(bb) + fact2 * shconv(bb) + fact2 * (divu - dmdphii) * shlu(bb)) * DetJ * gw(igauss)
        !       xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        !         fact2 * shlu(aa) * drhodphi * rm(:) / rhoi * DetJ * gw(igauss) + &
        !         fact2 * shconv(aa) * tauM * tauM * rLi(:) * drhodphi * DetJ * gw(igauss) + &
        !         fact2 * rhoi * shconv(aa) * tauM * drhodphi * rm(:) / rhoi * DetJ * gw(igauss) + &
        !         fact2 * shgradgu(aa, :) * divu * (-2d0/3) * dmudphi * DetJ * gw(igauss) + &
        !         fact2 * shgradgu(aa, :) * divu * drhodphi * tauC / rhoi * DetJ * gw(igauss)

        !       xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
        !         fact2 * sum(shgradgu(aa, :) * (duidxi(1, :) + duidxi(:, 1))) * dmudphi * DetJ * gw(igauss)
        !       xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
        !         fact2 * sum(shgradgu(aa, :) * (duidxi(2, :) + duidxi(:, 2))) * dmudphi * DetJ * gw(igauss)
        !       xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
        !         fact2 * sum(shgradgu(aa, :) * (duidxi(3, :) + duidxi(:, 3))) * dmudphi * DetJ * gw(igauss)

        !       xPLSebe(aa, bb) = xPLSebe(aa, bb) + &
        !         fact2 * shlu(aa) * drhodphi * shlu(bb) * (1/rhoa - 1/rhow) * DetJ * gw(igauss) + &
        !         fact2 * sum(shgradgu(aa, :) * rm(:)) / rhoi * tauM * drhodphi * shlu(bb) * DetJ * gw(igauss)
        !       
        !       xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
        !         fact2 * tauls * shgradgu(aa, :) * shlu(bb) * res_phic_tmp1 * DetJ * gw(igauss) + &
        !         fact2 * (shlu(aa) + tauls * shconv(aa)) * shlu(bb) * dphidxi(:) * DetJ * gw(igauss) + &
        !         fact2 * (shlu(aa) + tauls * shconv(aa)) * phii * shgradgu(bb, :) * DetJ * gw(igauss)
        !     end do
        !   end do
        ! endif
        ! if(ismaster) write(*,*) "Assemble NSLS Jac"
        call e3LHS_3D_fluid_quenching(&
          nshl, ui, umi, aci, pri, duidxi, dpridxi, &
          dphidxi, dphidxidxj, rphii, phii, &
          rLi, fi, rhoi, mui, Ti, &
          tauM, tauP, tauLS, tauC, tauBar, tauBar1, &
          k_dc, k_dc_phi, gw(igauss), shlu, shgradgu, &
          shhessgu, xKebe11, xGebe, xDebe1, xMebe, &
          xLSebe, xLSUebe, xULSebe, xPLSebe)
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_fluid) then
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) myid, "Before assemble NSLS RHS", iel, NELEM
        tmp1(:) = rhoi *aci(:) + &
          rhoi * uadvi_full(1) * duidxi(:, 1) + &
          rhoi * uadvi_full(2) * duidxi(:, 2) + &
          rhoi * uadvi_full(3) * duidxi(:, 3) - fi(:) -&
          uprime(:) * (vdot * rhoi + drhodphi * sum(uadvi_full(:) * dphidxi(:)))
        do aa = 1, NSHL
          RHSu(:, aa) = RHSu(:, aa) - &
            shlu(aa) * tmp1(:) * DetJ * gw(igauss) - &
            shgradgu(aa, :) * (-pri + (tauC + lambda) * divu) * DetJ * gw(igauss) - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 1) + duidxi(:, 1))) * DetJ * gw(igauss) + &
            rhoi * shconv_full(aa) * uprime(:) * DetJ * gw(igauss)
          RHSp(aa) = RHSp(aa) - &
            shlu(aa) * (divu - vdot) * DetJ * gw(igauss) + &
            sum(shgradgu(aa, :) * uprime(:)) * DetJ * gw(igauss)
          RHSphi(aa) = RHSphi(aa) - &
            shlu(aa) * (rphii + sum(uadvi_full(:) * dphidxi(:)) + phii * vdot - mdot / rhoa) * DetJ * gw(igauss) - &
            shconv_full(aa) * tauls * res_phic_tmp1 * DetJ * gw(igauss) - &
            k_dc_phi * sum(shgradgu(aa, :) * dphidxi(:)) * DetJ * gw(igauss)
        enddo
      end if

    end do

    ! Apply Dirichlet BCs
    call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                  xMebe, Rhsu, Rhsp, &
                  xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
    !    call MPI_barrier(MPI_COMM_WORLD, mpi_err)
    !if(ismaster) write(*,*) "Evaluate BC", NSHL
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                            xLSebe, xLSUebe, xULSebe, xPLSebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      call LocalToGlobalNSVOF_3D(RHSGu, RHSGp, RHSGls, &
                                 NNODE, NSD, NSHL, maxNSHL, &
                                 IEN(iel, :), rhsu, rhsp, rhsphi)
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

end subroutine IntElmAss_NSVOF_Quenching_VMS
!======================================================================
!
!======================================================================

subroutine IntElmAss_NSVOF_Quenching_STAB(&
  config, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  use configuration
  implicit none

  type(ConfigType), intent(in) :: config

  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), rphigAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
                         ! Local variables
  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
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

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi(NSD), uadvi_full(NSD)
  real(8) :: Ti, rTi, dTdxi(NSD)
  real(8) :: ns_kdc

  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl, gcfl, cfl_loc
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)
  logical :: is_fluid

  real(8) :: rhoi, mui, cpi, hki
  real(8) :: mdot, vdot, phic, dmdphii, drhodphi, dmudphi
  real(8) :: divu, lambda, rm(NSD), tmp1(NSD)
  real(8) :: fact1, fact2

  volm = 0.0d0
  vol_ex = 0.0d0
  cfl = 0.0d0

  fact1 = almi
  fact2 = alfi * delt * gami

  drhodphi = rhow - rhoa
  dmudphi = muw - mua

  NGAUSS = -1
  NSHL = -1
!  rhsgq = 0d0
!  lhsgq = 0d0

  ! loop over elements
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

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

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
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

    is_fluid = ELM_ID(iel) == 102

    fl = 0.0d0
!    do i = 1, NSHL
!      fl(i,3) = -1.0d0/(Fr**2.0d0)
!    end do

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      wl(i) = wg(idx)
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
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)

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
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdphii = -c_evap * rhow * (Ti - Ts) / Ts
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdphii = c_cond * rhoa * (Ti - Ts) / Ts
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
      lambda = -2d0 / 3 * mui

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
        dpridxi(:) - fi(:) - &
        mui * (duidxixj(:, 1, 1) + duidxixj(:, 2, 2) + duidxixj(:, 3, 3)) + &
        (mui+lambda) * (duidxixj(1, 1, :) + duidxixj(2, 2, :) + duidxixj(3, 3, :))

      res_phic_tmp1 = rphii + sum(uadvi(:) * dphidxi(:)) + phii * vdot - mdot / rhoa
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate on RESIDUAL", igauss, NGAUSS
      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)
      uadvi_full(:) = uadvi(:) + uprime(:)
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:) * shgradgu(aa, :))
        shconv_full(aa) = sum(uadvi_full(:) * shgradgu(aa, :))
      enddo
        ! call MPI_barrier(MPI_COMM_WORLD, mpi_err)
        ! if(ismaster) write(*,*) "Evaluate STAB", igauss, NGAUSS, &
        !     config%vms%use_taubar

      tauBar = 0.0d0
      if(config%vms%use_taubar) then
        call e3STAB_3D_NSVOF_TAUBAR(NSD, Gij, uadvi, uprime, tauBar)
      endif

      k_dc = 0.0d0
      if(abs(config%vms%NS_kdc_w) + abs(config%vms%NS_kdc_a)  > 0.0d0) then
        call prop_interp(NS_kdc_w, NS_kdc_a, phii, ns_kdc)
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
                  fact2 * shgradgu(aa, j) * shgradgu(bb, i) * (mui + k_dc) * DetJ * gw(igauss) + &
                  fact2 * shgradgu(aa, i) * shgradgu(bb, j) * (tauC+lambda) * DetJ * gw(igauss)
              enddo
              xKebe11((i-1)*NSD+i, aa, bb) = xKebe11((i-1)*NSD+i, aa, bb) + &
                fact1 * shlu(aa) * rhoi * shlu(bb) * DetJ * gw(igauss) + &
                fact2 * shlu(aa) * rhoi * shconv(bb) * DetJ * gw(igauss) + &
                fact2 * (mui + k_dc) * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJ * gw(igauss) + &
                fact1 * rhoi * shconv(aa) * tauM * rhoi * shlu(bb) * DetJ * gw(igauss) + &
                fact2 * rhoi * shconv(aa) * tauM * rhoi * shconv(bb) * DetJ * gw(igauss)
            enddo
            xGebe(:, aa, bb) = xGebe(:, aa, bb) - &
              fact2 * shgradgu(aa, :) * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * rhoi * shconv(aa) * tauM * shgradgu(bb, :) * DetJ * gw(igauss)

            xDebe1(:, aa, bb) = xDebe1(:, aa, bb) + &
              fact2 * shlu(aa) * shgradgu(bb, :) * DetJ * gw(igauss) + &
              fact1 * shgradgu(aa, :) * tauM * rhoi * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * shgradgu(aa, :) * tauM * rhoi * shconv(bb) * DetJ * gw(igauss)

            xMebe(aa, bb) = xMebe(aa, bb) + &
              fact2 * tauM * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJ * gw(igauss)

            xLSebe(aa, bb) = xLSebe(aa, bb)  + &
              (shlu(aa) + tauls * shconv(aa)) * (fact1 * shlu(bb) + fact2 * shconv(bb)) * DetJ * gw(igauss)
            xLSebe(aa, bb) = xLSebe(aa, bb)  + (shlu(aa) + tauls * shconv(aa)) * &
              fact2 * (shlu(bb) * divu - dmdphii * shlu(bb) / rhoa) * DetJ * gw(igauss)
            xLSebe(aa, bb) = xLSebe(aa, bb) + &
              fact2 * k_dc_phi * sum(shgradgu(aa, :) * shgradgu(bb, :)) * DetJ * gw(igauss)
              ! (fact1 * shlu(bb) + fact2 * shconv(bb) + fact2 * (divu - dmdphii) * shlu(bb)) * DetJ * gw(igauss)

            xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
              fact2 * shlu(aa) * drhodphi * rm(:) / rhoi * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * shconv(aa) * tauM * rLi(:) * drhodphi * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * rhoi * shconv(aa) * tauM * drhodphi * rm(:) / rhoi * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * shgradgu(aa, :) * divu * (-2d0/3) * dmudphi * shlu(bb) * DetJ * gw(igauss) + &
              fact2 * shgradgu(aa, :) * divu * drhodphi * tauC / rhoi * shlu(bb) * DetJ * gw(igauss)

            xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(1, :) + duidxi(:, 1))) * dmudphi * DetJ * gw(igauss)
            xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(2, :) + duidxi(:, 2))) * dmudphi * DetJ * gw(igauss)
            xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
              fact2 * sum(shgradgu(aa, :) * (duidxi(3, :) + duidxi(:, 3))) * dmudphi * DetJ * gw(igauss)

            xPLSebe(aa, bb) = xPLSebe(aa, bb) + &
              fact2 * shlu(aa) * dmdphii * shlu(bb) * (1/rhoa - 1/rhow) * DetJ * gw(igauss) + &
              fact2 * sum(shgradgu(aa, :) * rm(:)) / rhoi * tauM * drhodphi * shlu(bb) * DetJ * gw(igauss)
            
            xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
              fact2 * tauls * shgradgu(aa, :) * shlu(bb) * res_phic_tmp1 * DetJ * gw(igauss) + &
              fact2 * (shlu(aa) + tauls * shconv(aa)) * shlu(bb) * dphidxi(:) * DetJ * gw(igauss) + &
              fact2 * (shlu(aa) + tauls * shconv(aa)) * phii * shgradgu(bb, :) * DetJ * gw(igauss)
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
            shlu(aa) * tmp1(:) * DetJ * gw(igauss) - &
            shgradgu(aa, :) * (-pri + (tauC + lambda) * divu) * DetJ * gw(igauss) - &
            (mui + k_dc) * sum(shgradgu(aa, :) * (duidxi(:, 1) + duidxi(:, 1))) * DetJ * gw(igauss) + &
            rhoi * shconv(aa) * uprime(:) * DetJ * gw(igauss)
          RHSp(aa) = RHSp(aa) - &
            shlu(aa) * (divu - vdot) * DetJ * gw(igauss) + &
            sum(shgradgu(aa, :) * uprime(:)) * DetJ * gw(igauss)
          RHSphi(aa) = RHSphi(aa) - &
            (shlu(aa) + tauls * shconv(aa)) * res_phic_tmp1 * DetJ * gw(igauss) - &
            k_dc_phi * sum(shgradgu(aa, :) * dphidxi(:)) * DetJ * gw(igauss)
        enddo
      end if

    end do

    ! Apply Dirichlet BCs
    call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                  xMebe, Rhsu, Rhsp, &
                  xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
    !    call MPI_barrier(MPI_COMM_WORLD, mpi_err)
    !if(ismaster) write(*,*) "Evaluate BC", NSHL
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                            xLSebe, xLSUebe, xULSebe, xPLSebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      call LocalToGlobalNSVOF_3D(RHSGu, RHSGp, RHSGls, &
                                 NNODE, NSD, NSHL, maxNSHL, &
                                 IEN(iel, :), rhsu, rhsp, rhsphi)
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
subroutine IntElmAss_Tem_Quenching_VMS(&
  config, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  use configuration

  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), dphidtgAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
  ! Local variables
  real(8), parameter :: damp = 0.5d0

  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, kappa_str, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: RHStem(:), xTebe(:, :)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)

  real(8), allocatable :: rTl(:), Tl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:), tmp(:), drdTi(:)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), uadvi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi_full(NSD)
  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  !real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)

  !real(8) :: phi

  real(8) :: rTi, Ti, dTdxi(NSD), dTdxixj(NSD, NSD)
  real(8) :: res_tem1, res_tem2, tau_tem
  logical :: is_fluid
  real(8) :: fact1, fact2
  real(8) :: Se, mdot, vdot
  real(8) :: rhoi, mui, cpi, hki, rhocpi
  real(8) :: k_dc_tem
  real(8) :: dmdTi, phic, tmp1, dSedTi

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
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

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

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
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
    is_fluid = ELM_ID(iel) == 102
    fl = 0.0d0

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      dlold(i, :) = dgold(idx, :)
      wl(i) = wg(idx)
      ul(i, :) = ugAlpha(idx, :)
      ulold(i, :) = ugold(idx, :)
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
  
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)

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
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdTi = c_evap * (1-phic) * rhow / Ts
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdTi = c_cond * (phic) * rhoa / Ts
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
      fi(:) = rhoi * gravvec(:)

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      ! call a function to get residual for Tem
      Se = ((cpw - cpa) * (Ti - Ts) - Lh) * mdot 
      call e3int_rLi(NSD, rhoi, mui, aci, duidxi, uadvi, dpridxi, fi, duidxixj, rLi)


      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)
  
      call e3int_restem(NSD, rTi, Ti, dTdxi, dTdxixj, uadvi, Se, rhocpi, hki, res_tem1)
      call e3STAB_3D_TEM(NSD, Gij, uadvi, Delt, rhoi, cpi, hki, tau_tem)

      k_dc_tem = 0.0d0
      if(abs(config%vms%Tem_kdc) > 0d0) then
        call e3DC_scalar(NSD, config%vms%Tem_kdc, Gij, res_tem1, dTdxi, k_dc_tem)
      endif

      uadvi_full(:) = uadvi(:) + uprime(:)
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:)*shgradgu(aa, :))
        shconv_full(aa) = sum(uadvi_full(:)*shgradgu(aa, :))
      enddo
      tmp1 = rhocpi * vdot + (rhow*cpw-rhoa*cpa) * sum(uadvi_full(:) * dphidxi(:))
      dSedTi = (cpw - cpa) * mdot + ((cpw - cpa) * (Ti - Ts) - lh)* dmdTi
      drdTi(:) = rhocpi * (fact1 * shlu(:) + fact2 * shconv(:)) - fact2 * dSedTi * shlu(:)
      ! shconv_full(:) = shconv(:)
      ! tmp(:)  = shlu(:) + tau_tem * rhocpi * shconv(:)
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        ! e3LHS_3D_tem()
        ! xTebe(:, :) = 0.0d0
        if(is_fluid) then
          do aa = 1, NSHL
            do bb = 1,NSHL
              !xTebe(aa, bb) = xTebe(aa, bb) + shlu(aa) * rhocpi * &
              !                (fact1 * shlu(bb) + fact2 * shconv_full(bb)) * gw(igauss) * DetJ
              !xTebe(aa, bb) = xTebe(aa, bb) + fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
              !xTebe(aa, bb) = xTebe(aa, bb) + rhocpi * shconv_full(aa) * tau_tem * &
              !                rhocpi * (fact1 * shlu(bb) + fact2 * shconv(bb)) * gw(igauss) * DetJ
              xTebe(aa, bb) = xTebe(aa, bb) + &
                shlu(aa) * rhocpi * fact1 * shlu(bb) * gw(igauss) * DetJ + & 
                shlu(aa) * rhocpi * fact2 * shconv_full(bb) * gw(igauss) * DetJ -&
                !shlu(aa) * fact2 * shlu(bb) * (cpw - cpa) * mdot * gw(igauss) * DetJ - &
                !shlu(aa) * fact2 * shlu(bb) * ((cpw-cpa)*(Ti-Ts)-lh) * dmdTi * gw(igauss) * DetJ + &
                shlu(aa) * fact2 * shlu(bb) * dSedTi * gw(igauss) * DetJ + &
                fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
              xTebe(aa, bb) = xTebe(aa, bb) + &
                rhocpi * shconv_full(aa) * tau_tem * drdTi(bb) * gw(igauss) * DetJ + &
                shlu(aa) * tau_tem * drdTi(bb) * tmp1 * gw(igauss) * DetJ + &
                shlu(aa) * tau_tem * res_tem1 * rhocpi * fact2 * dmdTi * shlu(bb) *(1/rhoa - 1/rhow) * gw(igauss) * DetJ
                
            enddo
          enddo
        else
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) + &
                shlu(aa) * rhocpi * fact1 * shlu(bb) * gw(igauss) * DetJ + & 
                fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
            enddo
          enddo
        endif
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        ! e3RHS_3D_tem()
        ! Stabilized FEM
        if(is_fluid) then
          RHSTem(:) = RHSTem(:) - shlu(:) * (rhocpi * (rTi + sum(uadvi_full(:) * dTdxi(:))) - Se) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - rhocpi * shconv_full(:) * tau_tem * res_tem1 * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki+k_dc_tem) * shgradgu(:, 1) * dTdxi(1) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki+k_dc_tem) * shgradgu(:, 2) * dTdxi(2) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki+k_dc_tem) * shgradgu(:, 3) * dTdxi(3) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - shlu(:) * tau_tem * res_tem1 * tmp1 * gw(igauss) * DetJ
        else
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * rTi * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 1) * dTdxi(1) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 2) * dTdxi(2) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 3) * dTdxi(3) * gw(igauss) * DetJ
        endif
        ! RHSTem(:) = RHSTem(:) + shlu(:) * sum(uadvi(:) * dphidxi(:)) * (rhow * cpw - rhoa * cpa) * gw(igauss) * DetJ
        ! if(isnan(sum(RHSTem))) write(*,*) "DEBUG", myid, tau_tem, rhoi, cpi, hki, res_tem1, dTdxi
      end if
    end do

    ! Apply Dirichlet BCs
    call BCLHS_tem(nshl, iel, xTebe, RHSTem)
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      ! Assemble global LHS Matrix
      call FillSparseMat_tem(nshl, iel, xTebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! Assemble load RHS vector
      do aa = 1, NSHL
        RHSGTEM(IEN(iel, aa)) = RHSGTEM(IEN(iel, aa)) + Rhstem(aa)
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


end subroutine IntElmAss_Tem_Quenching_VMS
!======================================================================
!
!======================================================================
subroutine IntElmAss_Tem_Quenching_STAB(&
  config, &
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  use configuration

  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), dphidtgAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
  ! Local variables
  real(8), parameter :: damp = 0.5d0

  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, kappa_str, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: RHStem(:), xTebe(:, :)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)

  real(8), allocatable :: rTl(:), Tl(:)

  real(8), allocatable :: gp(:, :), gw(:)
  real(8), allocatable :: shconv(:), shconv_full(:), tmp(:), drdTi(:)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), rphii, &
             fi(NSD), uadvi(NSD), xi(NSD), ddidxi(NSD, NSD), phii, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: uadvi_full(NSD)
  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  !real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)

  !real(8) :: phi

  real(8) :: rTi, Ti, dTdxi(NSD), dTdxixj(NSD, NSD)
  real(8) :: res_tem1, res_tem2, tau_tem
  logical :: is_fluid
  real(8) :: fact1, fact2
  real(8) :: Se, mdot, vdot
  real(8) :: rhoi, mui, cpi, hki, rhocpi
  real(8) :: k_dc_tem
  real(8) :: dmdTi, phic, tmp1, dSedTi

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
  do iel = 1, NELEM
!  write(*,*) "iel:",iel, NELEM
    if (NSHL /= ELMNSHL(iel)) then

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

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
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
    is_fluid = ELM_ID(iel) == 102
    fl = 0.0d0

    ! Get local solution arrays
    do i = 1, NSHL
      idx = IEN(iel, i)
      xl(i, :) = xg(idx, :)
      dl(i, :) = dgAlpha(idx, :)
      dlold(i, :) = dgold(idx, :)
      wl(i) = wg(idx)
      ul(i, :) = ugAlpha(idx, :)
      ulold(i, :) = ugold(idx, :)
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
  
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, 0)

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
      if(Ti > Ts) then
        mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
        dmdTi = c_evap * (1-phic) * rhow / Ts
      else
        mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
        dmdTi = c_cond * (phic) * rhoa / Ts
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
      fi(:) = rhoi * gravvec(:)

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      ! call a function to get residual for Tem
      Se = ((cpw - cpa) * (Ti - Ts) - Lh) * mdot 
      call e3int_rLi(NSD, rhoi, mui, aci, duidxi, uadvi, dpridxi, fi, duidxixj, rLi)


      call e3STAB_3D_NSVOF(NSD, Gij, Delt, uadvi, rhoi, mui, tauM, tauP, tauC, tauLS)
      uprime(:) = -tauM * rLi(:)
  
      call e3int_restem(NSD, rTi, Ti, dTdxi, dTdxixj, uadvi, Se, rhocpi, hki, res_tem1)
      call e3STAB_3D_TEM(NSD, Gij, uadvi, Delt, rhoi, cpi, hki, tau_tem)

      k_dc_tem = 0.0d0
      if(abs(config%vms%Tem_kdc) > 0d0) then
        call e3DC_scalar(NSD, config%vms%Tem_kdc, Gij, res_tem1, dTdxi, k_dc_tem)
      endif

      ! uadvi_full(:) = uadvi(:) + uprime(:)
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:)*shgradgu(aa, :))
        ! shconv_full(aa) = sum(uadvi_full(:)*shgradgu(aa, :))
      enddo
      ! tmp1 = rhocpi * vdot + (rhow*cpw-rhoa*cpa) * sum(uadvi_full(:) * dphidxi(:))
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
              xTebe(aa, bb) = xTebe(aa, bb) + &
                tmp(aa) * rhocpi * (fact1 * shlu(bb) + fact2 * shconv(bb)) * gw(igauss) * DetJ - &
                tmp(aa) * fact2 * dSedTi * shlu(bb) * gw(igauss) * DetJ + &
                fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
                
            enddo
          enddo
        else
          do aa = 1, NSHL
            do bb = 1,NSHL
              xTebe(aa, bb) = xTebe(aa, bb) + &
                tmp(aa) * rhocpi * fact1 * shlu(bb) * gw(igauss) * DetJ - &
                fact2 * hki * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
            enddo
          enddo
        endif
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        if(is_fluid) then
          RHSTem(:) = RHSTem(:) - tmp(:) * res_tem1 * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 1) * dTdxi(1) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 2) * dTdxi(2) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - (hki + k_dc_tem) * shgradgu(:, 3) * dTdxi(3) * gw(igauss) * DetJ
        else
          RHSTem(:) = RHSTem(:) - shlu(:) * rhocpi * rTi * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 1) * dTdxi(1) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 2) * dTdxi(2) * gw(igauss) * DetJ
          RHSTem(:) = RHSTem(:) - hki * shgradgu(:, 3) * dTdxi(3) * gw(igauss) * DetJ
        endif
      end if
    end do

    ! Apply Dirichlet BCs
    call BCLHS_tem(nshl, iel, xTebe, RHSTem)
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      ! Assemble global LHS Matrix
      call FillSparseMat_tem(nshl, iel, xTebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! Assemble load RHS vector
      do aa = 1, NSHL
        RHSGTEM(IEN(iel, aa)) = RHSGTEM(IEN(iel, aa)) + Rhstem(aa)
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
