!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_weak(config, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                                assemble_tensor_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  use configuration
  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), rphigAlpha(NNODE)

  ! Local variables
  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: dl(:, :), ul(:, :), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), &
                          xl(:, :), phil(:), rphil(:), wl(:), dphidtl(:), phi_bgl(:)
  real(8), allocatable :: xKebe11(:, :, :), xMebe(:, :), xGebe(:, :, :), &
                          xDebe1(:, :, :), xDebe2(:, :, :), &
                          xLSebe(:, :), xLSUebe(:, :, :), &
                          xULSebe(:, :, :), xPLSebe(:, :)
  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsm(:, :), Rhsphi(:)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, outflow, traction

  real(8), allocatable :: gp(:, :), gw(:), mgp(:, :)
  real(8) :: ui(NSD), pri, duidxi(NSD, NSD), phii, rphii, dphidxi(NSD), &
             Gij(NSD, NSD), Ginv(NSD, NSD), ti(NSD), nor(NSD), &
             bui(NSD), umi(NSD), dxidx(NSD, NSD), dphidx(NSD), xi(3)
  real(8) :: tmp1(NSD), tmp2(NSD, NSD)
  real(8) :: fact1, fact2, tauB, tauNor, gi(NSD), gphi, tauBLS
  real(8) :: gneg, wave_u, wave_phi
  real(8) :: e3(3)
  real(8) :: rhoi, mui
  integer :: NGAUSSf
  e3 = 0d0
  e3(3) = -1d0
  ! "shb" will be the shape function array while "shbg" will hold the
  ! gradients of the shape functions
  fact1 = almi
  fact2 = alfi*gami*Delt

  outflow = 0 !not outflow boundary
  traction = 0 ! no traction BC

  gi = 0.0d0

  NSHL = -1

  ! loop over faces
  do b = 1, NBOUND

    if (BCugType(b, 1) /= 2 .and. &
        BCugType(b, 2) /= 2 .and. &
        BCugType(b, 3) /= 2) then
      cycle
    endif


    do ifac = 1, bound(b)%NFACE

      iel = bound(b)%F2E(ifac)

      if (NSHL /= ELMNSHL(iel)) then
        if (NSHL >= 0) then
          deallocate (shlu, shgradgu, shhessgu, &
                      dl, ul, acl, uml, acml, pl, xl, phil, rphil, wl, phi_bgl, &
                      xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                      Rhsu, Rhsp, Rhsm, xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
        end if

        NSHL = ELMNSHL(iel)
        NGAUSSf = bound(b)%NGAUSSB(ifac)

        allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD), &
                  dl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), uml(NSHL, NSD), &
                  acml(NSHL, NSD), pl(NSHL), xl(NSHL, NSD), phil(NSHL), rphil(NSHL), &
                  wl(NSHL), &
                  xKebe11(NSD*NSD, NSHL, NSHL), xMebe(NSHL, NSHL), &
                  xGebe(NSD, NSHL, NSHL), xDebe1(NSD, NSHL, NSHL), &
                  xDebe2(NSD, NSHL, NSHL), &
                  Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsm(NSD, NSHL), &
                  xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL), &
                  xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL), Rhsphi(NSHL), phi_bgl(NSHL))
        allocate (gp(NGAUSSf, 2), gw(NGAUSSf), mgp(NGAUSSf, 3))
      end if

      call genGPandGWb(gp, gw, NGAUSSf)
      call genGPMap(nshl, NGAUSSf, bound(b)%FACE_OR(ifac), iga, mgp)

      ! Get local solution arrays
      do i = 1, NSHL
        j = IEN(iel, i)
        xl(i, :) = xg(j, :)
        dl(i, :) = dgAlpha(j, :)
        ul(i, :) = ugAlpha(j, :)
        acl(i, :) = acgAlpha(j, :)
        uml(i, :) = ugmAlpha(j, :)
        acml(i, :) = acgmAlpha(j, :)
        pl(i) = pgAlpha(j)
        phil(i) = phigAlpha(j)
        rphil(i) = rphigAlpha(j)
        wl(i) = wg(j)
        ! phi_bgl(i) = phi_bg(j)
      end do

      ! initialize local resistance matrix
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      xLSUebe = 0d0
      xPLSebe = 0d0
      xULSebe = 0d0

      xLSebe = 0.0d0

      ! initialize local load vector
      Rhsu = 0.0d0
      Rhsp = 0.0d0
      Rhsm = 0.0d0

      Rhsphi = 0.0d0

      do igauss = 1, NGAUSSf

        call eval_faceshape(nshl, iel, gp(igauss, :), mgp(igauss, :), &
                            bound(b)%FACE_OR(ifac), &
                            xl, dl, wl, shlu, shgradgu, dxidx, &
                            Gij, Ginv, nor)

        ! Interpolate
        call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
        call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
        call e3int_qr(NSHL, NSD, 1, shlu, rphil, rphii)
        call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui)
        call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi)
        call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi)

        call e3int_qr_grad(NSHL, NSD, 1, shgradgu, phil, dphidxi)
        call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

        call prop_interp(rhow, rhoa, phii, rhoi)
        call prop_interp(muw, mua, phii, mui)
        ! velocity on the rotor is mesh velocity, not zero...
        gi = umi
        gphi = 0d0


        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
          call e3bSTAB_weak(tauB, tauNor, ui - umi, nor, dxidx, rhoi, mui)

          call e3bRHS_weak( &
            nshl, nor, tauB, tauNor, gw(igauss), &
            shlu, shgradgu, ui, umi, pri, duidxi, &
            gi, rhoi, mui, &
            Rhsu, Rhsp, ti, tmp1, tmp2)

        end if

        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
          call e3bLHS_weak( &
            nshl, ui, umi, duidxi, rhoi, mui, &
            tauB, tauNor, &
            gw(igauss), shlu, shgradgu, xKebe11, &
            xGebe, xDebe1, nor)
        end if

      end do
      call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp, &
                    xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                              xLSebe, xLSUebe, xULSebe, xPLSebe)
      endif
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)
      endif

    end do


  end do

  if (NSHL >= 0) then
    deallocate (shlu, shgradgu, shhessgu, &
                dl, ul, acl, uml, acml, xl, phil, rphil, wl, &
                xKebe11, xMebe, xGebe, xDebe1, xDebe2, &
                Rhsu, Rhsp, Rhsm, Rhsphi, phi_bgl, xLSebe, xULSebe, xLSUebe, xPLSebe, &
                gp, gw, mgp)
  end if

end subroutine FaceAssembly_NS_weak
!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_weak_CF(config, dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                   acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                                   assemble_tensor_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  use configuration
  implicit none

  type(ConfigType), intent(in) :: config
  integer, intent(in) :: assemble_tensor_flag
  real(8), intent(in) :: dgAlpha(NNODE, NSD), ugAlpha(NNODE, NSD), &
                         acgAlpha(NNODE, NSD), ugmAlpha(NNODE, NSD), &
                         acgmAlpha(NNODE, NSD), pgAlpha(NNODE), &
                         phigAlpha(NNODE), rphigAlpha(NNODE)

  ! Local variables
  real(8), allocatable :: shlu(:), shgradgu(:, :), shhessgu(:, :, :)

  real(8), allocatable :: dl(:, :), ul(:, :), acl(:, :), &
                          uml(:, :), acml(:, :), pl(:), &
                          xl(:, :), phil(:), rphil(:), wl(:), dphidtl(:), phi_bgl(:)
  real(8), allocatable :: xKebe11(:, :, :), xMebe(:, :), xGebe(:, :, :), &
                          xDebe1(:, :, :), xDebe2(:, :, :), &
                          xLSebe(:, :), xLSUebe(:, :, :), &
                          xULSebe(:, :, :), xPLSebe(:, :)
  real(8), allocatable :: Rhsu(:, :), Rhsp(:), Rhsm(:, :), Rhsphi(:)

  integer :: b, ifac, iel, igauss, i, j, k, hess_flag, aa, bb, nshl, outflow, traction
  integer :: ii, jj

  real(8), allocatable :: gp(:, :), gw(:), mgp(:, :)
  real(8) :: ui(NSD), pri, duidxi(NSD, NSD), phii, rphii, dphidxi(NSD), &
             Gij(NSD, NSD), Ginv(NSD, NSD), ti(NSD), nor(NSD), &
             bui(NSD), umi(NSD), dxidx(NSD, NSD), dphidx(NSD), xi(3)
  real(8) :: tmp1(NSD), tmp2(NSD, NSD)
  real(8) :: fact1, fact2, tauB, tauNor, gi(NSD), gphi, tauBLS
  real(8) :: gneg, wave_u, wave_phi
  real(8) :: e3(3)
  real(8) :: rhoi, mui
  real(8) :: drhodphii, dmudphii, lambda
  real(8) :: unor, upos, uneg, gnor, divu
  real(8) :: DetJgw
  real(8), allocatable :: tmp(:, :)
  integer :: NGAUSSf
  e3 = 0d0
  e3(3) = -1d0
  ! "shb" will be the shape function array while "shbg" will hold the
  ! gradients of the shape functions
  fact1 = almi
  fact2 = alfi*gami*Delt

  outflow = 0 !not outflow boundary
  traction = 0 ! no traction BC


  gi = 0.0d0

  NSHL = -1

  ! loop over faces
  do b = 1, NBOUND

    if (BCugType(b, 1) /= 2 .and. &
        BCugType(b, 2) /= 2 .and. &
        BCugType(b, 3) /= 2) then
      cycle
    endif


    do ifac = 1, bound(b)%NFACE

      iel = bound(b)%F2E(ifac)

      if (NSHL /= ELMNSHL(iel)) then
        if (NSHL >= 0) then
          deallocate (shlu, shgradgu, shhessgu)
          deallocate (dl, ul, acl, uml)
          deallocate (acml, pl, xl, phil, rphil)
          deallocate (wl)
          deallocate (xKebe11, xMebe)
          deallocate (xGebe, xDebe1)
          deallocate (Rhsu, Rhsp, Rhsphi)
          deallocate (xLSebe, xLSUebe)
          deallocate (xULSebe, xPLSebe)
          deallocate (gp, gw, mgp)
          deallocate (tmp)
        end if

        NSHL = ELMNSHL(iel)
        NGAUSSf = bound(b)%NGAUSSB(ifac)

        allocate (shlu(NSHL), shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD))
        allocate (dl(NSHL, NSD), ul(NSHL, NSD), acl(NSHL, NSD), uml(NSHL, NSD))
        allocate (acml(NSHL, NSD), pl(NSHL), xl(NSHL, NSD), phil(NSHL), rphil(NSHL))
        allocate (wl(NSHL))
        allocate (xKebe11(NSD*NSD, NSHL, NSHL), xMebe(NSHL, NSHL))
        allocate (xGebe(NSD, NSHL, NSHL), xDebe1(NSD, NSHL, NSHL))
        allocate (Rhsu(NSD, NSHL), Rhsp(NSHL), Rhsphi(NSHL))
        allocate (xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHL, NSHL))
        allocate (xULSebe(NSD, NSHL, NSHL), xPLSebe(NSHL, NSHL))
        allocate (gp(NGAUSSf, 2), gw(NGAUSSf), mgp(NGAUSSf, 3))
        allocate (tmp(NSHL, NSHL))
      end if

      call genGPandGWb(gp, gw, NGAUSSf)
      call genGPMap(nshl, NGAUSSf, bound(b)%FACE_OR(ifac), iga, mgp)

      ! Get local solution arrays
      do i = 1, NSHL
        j = IEN(iel, i)
        xl(i, :) = xg(j, :)
        dl(i, :) = dgAlpha(j, :)
        ul(i, :) = ugAlpha(j, :)
        acl(i, :) = acgAlpha(j, :)
        uml(i, :) = ugmAlpha(j, :)
        acml(i, :) = acgmAlpha(j, :)
        pl(i) = pgAlpha(j)
        phil(i) = phigAlpha(j)
        rphil(i) = rphigAlpha(j)
        wl(i) = wg(j)
        ! phi_bgl(i) = phi_bg(j)
      end do

      ! initialize local resistance matrix
      xKebe11 = 0.0d0
      xGebe = 0.0d0
      xDebe1 = 0.0d0
      xMebe = 0.0d0

      xLSUebe = 0d0
      xPLSebe = 0d0
      xULSebe = 0d0

      xLSebe = 0.0d0

      ! initialize local load vector
      Rhsu = 0.0d0
      Rhsp = 0.0d0
      Rhsphi = 0.0d0


      do igauss = 1, NGAUSSf

        call eval_faceshape(nshl, iel, gp(igauss, :), mgp(igauss, :), &
                            bound(b)%FACE_OR(ifac), &
                            xl, dl, wl, shlu, shgradgu, dxidx, &
                            Gij, Ginv, nor)

        ! Interpolate
        call e3int_qr(NSHL, NSD, 1, shlu, pl, pri)
        call e3int_qr(NSHL, NSD, 1, shlu, phil, phii)
        call e3int_qr(NSHL, NSD, 1, shlu, rphil, rphii)
        call e3int_qr(NSHL, NSD, NSD, shlu, ul, ui)
        call e3int_qr(NSHL, NSD, NSD, shlu, uml, umi)
        call e3int_qr(NSHL, NSD, NSD, shlu, xl, xi)

        call e3int_qr_grad(NSHL, NSD, 1, shgradgu, phil, dphidxi)
        call e3int_qr_grad(NSHL, NSD, NSD, shgradgu, ul, duidxi)

        call prop_interp(rhow, rhoa, phii, rhoi)
        call prop_interp(muw, mua, phii, mui)
        ! velocity on the rotor is mesh velocity, not zero...
        gi = umi
        gphi = 0d0

        unor = sum((ui - umi)*nor)  ! u \cdot n
        upos = 0.5d0*(unor + abs(unor))
        uneg = 0.5d0*(unor - abs(unor))

        ! Absolute normal for the rest
        gnor = sum(gi*nor)
        unor = sum(ui*nor)  ! u \cdot n

        lambda = -2d0/3d0 * mui

        drhodphii = rhoa - rhow
        dmudphii = mua - muw
        if(phii < 0d0 .or. phii > 1.0d0) then
          drhodphii = 0d0
          dmudphii = 0d0
        endif

        divu = duidxi(1,1) + duidxi(2,2) + duidxi(3,3)
        DetJgw = gw(igauss) * DetJb
        
        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
          call e3bSTAB_weak(tauB, tauNor, ui - umi, nor, dxidx, rhoi, mui)

          ! call e3bRHS_weak( &
          !   nshl, nor, tauB, tauNor, gw(igauss), &
          !   shlu, shgradgu, ui, umi, pri, duidxi, &
          !   gi, rhoi, mui, &
          !   Rhsu, Rhsp, ti, tmp1, tmp2)
          ti(:) = -pri * nor(:) + &
            mui * matmul(duidxi + transpose(duidxi), nor) + &
            lambda * divu * nor(:)
          do aa = 1, NSHL
            RHSu(:, aa) = RHSu(:, aa) - &
              shlu(aa) * (-ti(:) - rhoi * uneg * (ui(:) - gi(:)) + tauB * (ui(:)-gi(:))) * DetJgw
            RHSu(:, aa) = RHSu(:, aa) + &
                mui * dot_product(shgradgu(aa, :), nor(:)) * (ui(:) - gi(:)) * DetJgw
            RHSu(:, aa) = RHSu(:, aa) + &
                mui * dot_product(shgradgu(aa, :), ui(:) - gi(:)) * nor(:) * DetJgw
            RHSu(:, aa) = RHSu(:, aa) + &
                shgradgu(aa, :) * lambda * (unor - gnor) * DetJgw
            RHSp(aa) = RHSp(aa) + shlu(aa) * (unor - gnor) * DetJgw
          enddo

        end if

        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
          ! call e3bLHS_weak( &
          !   nshl, ui, umi, duidxi, rhoi, mui, &
          !   tauB, tauNor, &
          !   gw(igauss), shlu, shgradgu, xKebe11, &
          !   xGebe, xDebe1, nor)
          do bb = 1, NSHL
            do aa = 1, NSHL
              tmp(aa, bb) = -shlu(aa)*mui*sum(shgradgu(bb, :)*nor(:)) &
                            - sum(shgradgu(aa, :)*nor(:))*mui*shlu(bb) &
                            + shlu(aa)*tauB*shlu(bb) &
                            - shlu(aa)*uneg*rhoi*shlu(bb)

            enddo
          enddo

          do bb = 1, NSHL
            do aa = 1, NSHL
              do ii = 1, NSD
                do jj = 1, NSD
                  xKebe11((ii-1)*NSD+jj, aa, bb) = xKebe11((ii-1)*NSD+jj, aa, bb) + &
                    fact2 * (-shlu(aa) * mui * shgradgu(bb, ii) * nor(jj) &
                             -shlu(aa) * lambda * shgradgu(bb, jj) * nor(ii) &
                             -shgradgu(aa, jj) * mui * shlu(bb) * nor(ii) &
                             -shgradgu(aa, ii) * lambda * shlu(bb) * nor(jj)) * &
                             DetJgw
                enddo
                xKebe11((ii-1)*NSD+ii, aa, bb) = xKebe11((ii-1)*NSD+ii, aa, bb) +&
                    fact2 * tmp(aa, bb) * DetJgw
              enddo
              xDebe1(:, aa, bb) = xDebe1(:, aa, bb) - &
                fact2 * shlu(aa) * shlu(bb) * nor(:) * DetJgw

              xGebe(:, aa, bb) = xGebe(:, aa, bb) + &
                fact2 * shlu(aa) * shlu(bb) * nor(:) * DetJgw

              xULSebe(:, aa, bb) = xULSebe(:, aa, bb) - &
                fact2 * shlu(aa) * dmudphii * shlu(bb) * matmul(duidxi+transpose(duidxi), nor) * DetJgw - &
                fact2 * shlu(aa) * (-2d0/3) * dmudphii * shlu(bb) * divu * nor(:) * DetJgw - &
                fact2 * dot_product(shgradgu(aa, :), nor(:)) * (ui(:) - gi(:)) * dmudphii * shlu(bb) * DetJgw -&
                fact2 * dot_product(shgradgu(aa, :), ui(:) - gi(:)) * nor(:) * dmudphii * shlu(bb) * DetJgw -&
                fact2 * shgradgu(aa, :) * (-2d0/3) * dmudphii * shlu(bb) * (unor-gnor) * DetJgw - &
                fact2 * shlu(aa) * drhodphii * shlu(bb) * (ui(:)-gi(:)) * uneg * DetJgw + &
                fact2 * shlu(aa) * tauB/mui * dmudphii * shlu(bb) * (ui(:) - gi(:)) * DetJgw

            enddo
          enddo
        end if

      end do
      call BCLhs_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, Rhsu, Rhsp, &
                    xLSebe, xLSUebe, xULSebe, xPLSebe, Rhsphi)
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
        call FillSparseMat_3D(nshl, iel, xKebe11, xGebe, xDebe1, xMebe, &
                              xLSebe, xLSUebe, xULSebe, xPLSebe)
      endif
      if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
        call LocaltoGlobal_3D(nshl, iel, Rhsu, Rhsp)
      endif

    end do


  end do

  if (NSHL >= 0) then
    deallocate (shlu, shgradgu, shhessgu)
    deallocate (dl, ul, acl, uml)
    deallocate (acml, pl, xl, phil, rphil)
    deallocate (wl)
    deallocate (xKebe11, xMebe)
    deallocate (xGebe, xDebe1)
    deallocate (Rhsu, Rhsp, Rhsphi)
    deallocate (xLSebe, xLSUebe)
    deallocate (xULSebe, xPLSebe)
    deallocate (gp, gw, mgp)
    deallocate (tmp)

  end if

end subroutine FaceAssembly_NS_weak_CF
