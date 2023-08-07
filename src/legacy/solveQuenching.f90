!======================================================================
!
!======================================================================
subroutine assembleQuenching(assemble_tensor_flag, assemble_field_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars

  implicit none

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
    ! write(*,*) myid, "ug-alpha", sum(ugAlpha(:, :) ** 2), assemble_tensor_flag 
    call IntElmAss_NSVOF_Quenching(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                         acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                         TgAlpha, rTgAlpha, &
                         assemble_tensor_flag)
    ! write(*,*) myid, "RHSgu1", sum(RHSGu(:, :) ** 2), &
    !       sum(LHSK11**2), sum(lhsG**2), sum(lhsD1**2), sum(lhsM**2)
    call FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                              acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                              assemble_tensor_flag)
    ! write(*,*) myid, "RHSgu2", sum(RHSGu(:, :) ** 2), &
    !      sum(LHSK11**2), sum(lhsG**2), sum(lhsD1**2), sum(lhsM**2)
  end if
  if (iand(assemble_field_flag, ASSEMBLE_FIELD_TEM) > 0) then
    call IntElmAss_Tem_Quenching(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                       acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                       TgAlpha, rTgAlpha, &
                       assemble_tensor_flag)

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

subroutine IntElmAss_NSVOF_Quenching(&
  dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
  acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
  TgAlpha, rTgAlpha, &
  assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none

  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), dphidtgAlpha(NNODE), &
                         TgAlpha(NNODE), rTgAlpha(NNODE)
                         ! Local variables
  integer :: iel, igauss, aa, bb, i, j, hess_flag, idx, nshl, NGAUSS

  real(8) :: volm, vol_ex, uprime(NSD)
  real(8) :: tauM, tauC, tauBar, tauBar1, pri, k_dc, k_dc_phi, tauP, tauLS

  real(8), allocatable :: xKebe11(:, :, :), xGebe(:, :, :), xDebe1(:, :, :), &
                          xMebe(:, :), Rhsu(:, :), Rhsm(:, :), Rhsp(:), Rhsq(:), Rhsl(:), &
                          xLSebe(:, :), xLSUebe(:, :, :), xPLSebe(:, :), &
                          xULSebe(:, :, :), Rhsphi(:)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)
  real(8), allocatable :: Tl(:), rTl(:)

  real(8), allocatable :: gp(:, :), gw(:)

  integer, allocatable :: ibc_loc(:, :)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), dphidti, &
             fi(NSD), uadvi(NSD), uadvi_ls(NSD), xi(NSD), ddidxi(NSD, NSD), phi, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)
  real(8) :: Ti, rTi, dTdxi(NSD)

  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)
  logical :: is_fluid
  real(8) :: rhoi, mui

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
        deallocate (shlu, shgradgu, shhessgu, gp, gw, &
                    fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
                    xl, dumb, phil, ulold, dphidtl, &
                    xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsm, Rhsp, Rhsq, Rhsl, &
                    xLSebe, xLSUebe, xULSebe, xPLsebe, Rhsphi)
        deallocate (Tl, rTl)
        deallocate (ibc_loc)
      end if

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
      allocate (xKebe11(NSD*NSD, NSHL, NSHL), xGebe(NSD, NSHL, NSHL), &
                xDebe1(NSD, NSHL, NSHL), xMebe(NSHL, NSHL), &
                xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL), &
                xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL), &
                Rhsu(NSD, NSHL), Rhsm(NSD, NSHL), Rhsp(NSHL), Rhsq(NSHL), Rhsl(NSHL), &
                Rhsphi(NSHL), &
                shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), &
                fl(NSHL, NSD), dl(NSHL, NSD), ul(NSHL, NSD), wl(NSHL), &
                acl(NSHL, NSD), uml(NSHL, NSD), acml(NSHL, NSD), &
                pl(NSHL), dlold(NSHL, NSD), xl(NSHL, NSD), &
                dumb(NSHL, NSD), phil(NSHL), ulold(NSHL, NSD), dphidtl(NSHL), &
                gp(NGAUSS, NSD), gw(NGAUSS))
      allocate (Tl(NSHL), rTl(NSHL))
      allocate (ibc_loc(NSHL, NSD))
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

    dumb = 0.0d0

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
      Rhsm = 0.0d0
      Rhsp = 0.0d0
      Rhsq = 0.0d0
      Rhsl = 0.0d0
      Rhsphi = 0.0d0
    end if
    ! Loop over integration points (NGAUSS in each direction)
    do igauss = 1, NGAUSS

      ! Get Element Shape functions and their gradients
      ! initialize
      shlu = 0.0d0
      shgradgu = 0.0d0
      shhessgu = 0.0d0
      hess_flag = NS_hess_flag
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, hess_flag)

      call e3int_fluid(nshl, xl, dl, ul, acl, uml, acml, &
                       pl, fl, phil, shlu, shgradgu, &
                       shgradgu, shhessgu, dxidx, Ginv, &
                       di, ui, aci, umi, acmi, pri, fi, &
                       ddidxi, duidxi, duidxixj, dpridxi, &
                       phi, dphidxi, dphidxidxj, dphidtl, dphidti, xi, rLi)
      rTi = sum(rTl(:)*shlu(:))
      Ti = sum(Tl(:)*shlu(:))
      dTdxi(1) = sum(Tl(:)*shgradgu(:, 1))
      dTdxi(2) = sum(Tl(:)*shgradgu(:, 2))
      dTdxi(3) = sum(Tl(:)*shgradgu(:, 3))
      rho = rhow * phi + rhoa * (1d0 - phi)
      mu = muw * phi + mua * (1d0 - phi)
      cp = cpw * phi + cpa * (1d0 - phi)
      hk = kappaw * phi + kappaa * (1d0 - phi)
      if (.not. is_fluid) then
        rho = 2.7d3
        mu = 1d3
        cp = 921.0d0
        hk = 205.0d0
        uadvi(:) = 0d0
        ui(:) = 0d0
        umi(:) = 0d0
      end if
      rhoi = rho
      mui = mu
      fi(:) = rhoi * gravvec(:)
  
      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      uadvi_ls(:) = uadvi(:) + gravvec(:)*usettle ! usettle = 0d0
  
      tauM = 0.0d0; tauP = 0.0d0; tauC = 0.0d0; tauBar = 0.0d0; tauLS = 0.0d0

      call e3STAB_3D(Gij, Ginv, uadvi, uadvi_ls, rLi, &
                     tauM, tauP, tauLS, tauC, tauBar, tauBar1, uprime, cfl_loc, &
                     rho, mu)

      cfl(1) = max(cfl(1), cfl_loc(1))
!!!      cfl(2)  = max(cfl(2), cfl_loc(2))
      k_dc = 0.0d0
      k_dc_phi = 0.0d0
      call e3DC_beta2(uadvi, duidxi, Gij, rLi, dxidx, k_dc)

      res_phic_tmp1 = dphidti + sum(uadvi_ls(:)*dphidxi(:)) - (kappa*dphidxidxj(1, 1) &
                                                               + kappa*dphidxidxj(2, 2) &
                                                               + kappa*dphidxidxj(3, 3))
      res_phic_tmp2 = 0d0

      do i = 1, NSD
        do j = 1, NSD
          res_phic_tmp2 = res_phic_tmp2 + dphidxi(i)*Gij(i, j)*dphidxi(j)
        end do
      end do

      k_dc_phi = 1d0*abs(res_phic_tmp1)/(sqrt(res_phic_tmp2) + 0.0000000000001d0)
      k_dc_phi = 0d0

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0 .and. is_fluid) then

        call e3LHS_3D_fluid_quenching(&
          nshl, ui, umi, aci, pri, duidxi, dpridxi, &
          dphidxi, dphidxidxj, dphidti, phi, &
          rLi, rhoi, mui, &
          tauM, tauP, tauLS, tauC, tauBar, tauBar1, &
          k_dc, k_dc_phi, gw(igauss), shlu, shgradgu, &
          shhessgu, xKebe11, xGebe, xDebe1, xMebe, &
          xLSebe, xLSUebe, xULSebe, xPLSebe)

      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0 .and. is_fluid) then
        call e3RHS_3D_fluid_quenching(&
          nshl, ui, aci, umi, acmi, uadvi, &
          pri, rLi, fi, duidxi, ddidxi, &
          tauM, tauP, tauLS, tauC, tauBar, tauBar1, k_dc, k_dc_phi, &
          gw(igauss), shlu, shgradgu, uprime, &
          Rhsu, Rhsp, phi, &
          dpridxi, dphidxi, dphidxidxj, dphidti, &
          Ti, rTi, rhoi, mui, &
          RHSphi)

      end if
      ! if(ismaster) write(*,*) "flag = ", assemble_tensor_flag

    end do
    do aa = 1, NSHL
      do bb = 1,NSD
        ibc_loc(aa, bb) = ibc(ien(iel, aa), bb)
      enddo
    enddo

    ! Apply Dirichlet BCs
    ! if(myid == 1 .and. iel == 1794 .and. iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
    !   write(*,*) "xKebe11", xKebe11
    ! endif
    call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, &
                  xMebe, Rhsu, Rhsp, &
                  xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
      call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                            xLSebe, xLSUebe, xULSebe, xPLSebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)
      ! do aa = 1, NSHL
      !   RHSGLS(IEN(iel, aa)) = RHSGLS(IEN(iel, aa)) + Rhsphi(aa)
      ! end do
      call LocalToGlobalNSVOF_3D(RHSGu, RHSGp, RHSGls, &
                                 NNODE, NSD, NSHL,&
                                 IEN(iel, :), rhsu, rhsp, rhsphi)
    end if

  end do

  deallocate (xKebe11, xGebe, shlu, shgradgu, shhessgu, &
              fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
              xl, dumb, phil, ulold, dphidtl, &
              Rhsu, Rhsm, Rhsp, Rhsq, Rhsl, xDebe1, xMebe, gp, gw, &
              xLSebe, xLSUebe, xULSebe, Rhsphi)
  deallocate (Tl, rTl)
  deallocate (ibc_loc)
  if (.true.) then
    ! Find largest CFL-number and output to screen
    if (numnodes > 1) then
      gcfl = cfl
      call MPI_ALLReduce(gcfl, cfl, 2, MPI_DOUBLE_PRECISION, &
                         MPI_MAX, MPI_COMM_WORLD, mpi_err)
    end if
    if (ismaster) then
      write (*, '(40("-"))')
      write (*, *) "    CFL = ", cfl(1)
    end if
  end if

end subroutine IntElmAss_NSVOF_Quenching
!======================================================================
!
!======================================================================
subroutine IntElmAss_Tem_Quenching(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                         acgmAlpha, pgAlpha, phigAlpha, dphidtgAlpha, &
                         TgAlpha, rTgAlpha, &
                         assemble_tensor_flag)
  use aAdjKeep
  use commonvars
  use commonpars
  use mpi
  implicit none

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

  real(8), allocatable :: RHStem(:)
  real(8), allocatable :: xTebe(:, :)

  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: fl(:, :), dl(:, :), ul(:, :), wl(:), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), dlold(:, :), &
                          xl(:, :), dumb(:, :), phil(:), ulold(:, :), dphidtl(:)

  real(8), allocatable :: rTl(:), Tl(:)

  real(8), allocatable :: gp(:, :), gw(:)

  real(8) :: di(NSD), ui(NSD), aci(NSD), umi(NSD), acmi(NSD), dphidti, &
             fi(NSD), uadvi(NSD), uadvi_ls(NSD), xi(NSD), ddidxi(NSD, NSD), phi, &
             dphidxi(NSD), duidxi(NSD, NSD), dphidxidxj(NSD, NSD), &
             duidxixj(NSD, NSD, NSD), duiolddxi(NSD, NSD), &
             dxidx(NSD, NSD), dpridxi(NSD), rLi(NSD)

  real(8) :: Gij(NSD, NSD), Ginv(NSD, NSD)
  real(8) :: cfl(2), gcfl(2), cfl_loc(2)
  real(8) :: Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
  !real(8) :: res_phic_tmp1, res_phic_tmp2
  real(8) :: Qi, Si(3, 3), Omegai(3, 3)

  !real(8) :: phi

  real(8) :: rTi, Ti, dTdxi(NSD)
  real(8) :: res_tem1, res_tem2, tau_tem
  logical :: is_fluid
  real(8) :: fact1, fact2
  real(8), allocatable :: shconv(:), tmp(:)
  real(8) :: Se, mdot, vdot

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
        deallocate (shlu, shgradgu, shhessgu, gp, gw, &
                    fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
                    xl, dumb, phil, ulold, dphidtl)
        deallocate (rTl, Tl, RHStem, xTebe)
        deallocate (shconv, tmp)
      end if

      NSHL = ELMNSHL(iel)
      NGAUSS = ELMNGAUSS(iel)
      allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), &
                gp(NGAUSS, NSD), gw(NGAUSS), &
                fl(NSHL, NSD), dl(NSHL, NSD), ul(NSHL, NSD), wl(NSHL), &
                acl(NSHL, NSD), uml(NSHL, NSD), acml(NSHL, NSD), &
                pl(NSHL), dlold(NSHL, NSD), xl(NSHL, NSD), &
                dumb(NSHL, NSD), phil(NSHL), ulold(NSHL, NSD), dphidtl(NSHL))
      allocate (rTl(NSHL), Tl(NSHL), RHStem(NSHL), xTebe(NSHL, NSHL))
      allocate (shconv(NSHL), tmp(NSHL))

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
      hess_flag = NS_hess_flag
      call eval_shape(nshl, iel, gp(igauss, :), xl, dl, wl, shlu, &
                      shgradgu, shhessgu, dxidx, Gij, Ginv, hess_flag)

      ! Fluid
      phi = sum(phil(:)*shlu(:))
      rTi = sum(rTl(:)*shlu(:))
      Ti = sum(Tl(:)*shlu(:))
      dTdxi(1) = sum(Tl(:)*shgradgu(:, 1))
      dTdxi(2) = sum(Tl(:)*shgradgu(:, 2))
      dTdxi(3) = sum(Tl(:)*shgradgu(:, 3))
      ui(1) = sum(ul(:, 1)*shlu(:))
      ui(2) = sum(ul(:, 2)*shlu(:))
      ui(3) = sum(ul(:, 3)*shlu(:))
      umi(1) = sum(uml(:, 1)*shlu(:))
      umi(2) = sum(uml(:, 2)*shlu(:))
      umi(3) = sum(uml(:, 3)*shlu(:))
      rho = phi*rhow + (1 - phi)*rhoa
      mu = phi*muw + (1 - phi)*mua
      cp = phi*cpw + (1 - phi)*cpa
      hk = phi*kappaw + (1 - phi)*kappaa
      if (.not. is_fluid) then
        rho = 2.7d3
        mu = 1d3
        cp = 921.0d0
        hk = 205.0d0
        uadvi(:) = 0d0
        ui(:) = 0d0
        umi(:) = 0d0
      end if

      ! ALE Advective Velocity
      uadvi(:) = ui(:) - umi(:)
      uadvi_ls(:) = uadvi(:) + gravvec(:)*usettle
      tauM = 0.0d0; tauP = 0.0d0; tauC = 0.0d0; tauBar = 0.0d0; tauLS = 0.0d0
      do aa = 1, NSHL
        shconv(aa) = sum(uadvi(:)*shgradgu(aa, :))
      enddo
      tau_tem = 0.0
      do aa = 1, NSD
        do bb = 1, NSD
          tau_tem = tau_tem + uadvi(bb) *  Gij(bb, aa) * uadvi(aa)
        enddo
      enddo
      tau_tem = tau_tem + 4.0d0 / Delt ** 2
      tau_tem = tau_tem + 3.0d0 * (hk / rho / cp) ** 2 * sum(Gij ** 2)
      tau_tem = 1d0/rho/cp/sqrt(tau_tem)

      if(Ti > 1d2 + 2.73d2) then
        mdot = c_evap * (1 - phi) * rhow * (Ti - Ts) / Ts
      else
        mdot = c_cond * phi * rhoa * (Ti - Ts) / Ts
      endif
      vdot = mdot / rhoa - mdot / rhow

      Se = ((cpw - cpa) * (Ti - Ts) - Lh) * mdot 
      res_tem1 = rho*cp*(rTi + sum(uadvi(:)*dTdxi(:))) - Se
      ! resume here
      res_tem2 = 0d0

      do i = 1, NSD
        do j = 1, NSD
          res_tem2 = res_tem2 + dTdxi(i)*Gij(i, j)*dTdxi(j)
        end do
      end do

      hk = hk + 1.0d0 * abs(res_tem1) / sqrt(res_tem2 + 1d-10)
      !k_dc_phi = 1d0*abs(res_phic_tmp1)/(sqrt(res_phic_tmp2) + 0.0000000000001d0)
!        write(*,*) "kc:", k_dc, k_dc_phi
!        k_dc = 0d0
      !k_dc_phi = 0d0
      ! k_dc = 1d0 * abs(res_tem1) / (sqrt(res_tem2) + 1d-10)

      tmp(:)  = shlu(:) + tau_tem * rho * cp * shconv(:)
      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        ! e3LHS_3D_tem()
        ! xTebe(:, :) = 0.0d0
        do aa = 1, NSHL
          do bb = 1,NSHL
            xTebe(aa, bb) = xTebe(aa, bb) + tmp(aa) * fact1 * rho * cp * shlu(bb) * gw(igauss) * DetJ
            xTebe(aa, bb) = xTebe(aa, bb) + tmp(aa) * fact2 * rho * cp * shconv(bb) * gw(igauss) * DetJ
            xTebe(aa, bb) = xTebe(aa, bb) + fact2 * hk * sum(shgradgu(aa, :) * shgradgu(bb, :)) * gw(igauss) * DetJ
          enddo
        enddo
      end if

      if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        ! e3RHS_3D_tem()
        RHSTem(:) = RHStem(:) -  tmp(:) * res_tem1 * gw(igauss) * DetJ
        RHSTem(:) = RHSTem(:) - hk * shgradgu(:, 1) * dTdxi(1) * gw(igauss) * DetJ
        RHSTem(:) = RHSTem(:) - hk * shgradgu(:, 2) * dTdxi(2) * gw(igauss) * DetJ
        RHSTem(:) = RHSTem(:) - hk * shgradgu(:, 3) * dTdxi(3) * gw(igauss) * DetJ
      end if
    end do
    !if(sum(dTdxi(:) ** 2) > 0d0) then
    !  write(*,*) myid, iel, "dTdxi:", dTdxi(:), "R1=", RHSTem(:)
    !endif
    ! Apply Dirichlet BCs
    call BCLHS_tem(nshl, iel, xTebe, RHSTem)
    !if(sum(dTdxi(:) ** 2) > 0d0) then
    !  write(*,*) myid, iel, "dTdxi:", dTdxi(:), "R2=", RHSTem(:)
    !endif
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then

      ! Assemble global LHS Matrix

      !call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
      !                      xLSebe, xLSUebe, xULSebe, xPLSebe)
      call FillSparseMat_tem(nshl, iel, xTebe)
    end if
    if (iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
      ! Assemble load RHS vector
      do aa = 1, NSHL
        RHSGTEM(IEN(iel, aa)) = RHSGTEM(IEN(iel, aa)) + Rhstem(aa)
      end do
    end if

!    write(*,*) sum(Rhsphi)

  end do

  deallocate (shlu, shgradgu, shhessgu, &
              fl, dl, ul, wl, acl, uml, acml, pl, dlold, &
              xl, dumb, phil, ulold, dphidtl, &
              gp, gw)
  deallocate (rTl, Tl, RHStem, xTebe)
  deallocate(shconv, tmp)

end subroutine IntElmAss_Tem_Quenching
