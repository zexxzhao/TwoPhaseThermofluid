!======================================================================
!
!======================================================================
subroutine FaceAssembly_NS_weak(dgAlpha, ugAlpha, ugmAlpha, acgAlpha, &
                                acgmAlpha, pgAlpha, phigAlpha, rphigAlpha, &
                                assemble_tensor_flag)
  use aAdjKeep
  use mpi
  use commonvars
  use commonpars
  implicit none

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
  real(8) :: ui(NSD), pri, duidxi(NSD, NSD), phi, rphi, dphidxi(NSD), &
             Gij(NSD, NSD), Ginv(NSD, NSD), ti(NSD), nor(NSD), &
             bui(NSD), umi(NSD), dxidx(NSD, NSD), dphidx(NSD), xi(3)
  real(8) :: tmp1(NSD), tmp2(NSD, NSD)
  real(8) :: fact1, fact2, tauB, tauNor, gi(NSD), gphi, tauBLS
  real(8) :: gneg, wave_u, wave_phi
  real(8) :: e3(3)
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
        pri = sum(pl*shlu)
        phi = sum(phil*shlu)
        rphi = sum(rphil*shlu)
        do i = 1, NSD
          ui(i) = sum(ul(:, i)*shlu)
          umi(i) = sum(uml(:, i)*shlu)
          xi(i) = sum(xl(:, i)*shlu)
        end do
        rho = rhow * phi + rhoa * (1d0 - phi)
        mu = muw * phi + mua * (1d0 - phi)

        do j = 1, NSD
          dphidxi(j) = sum(phil*shgradgu(:, j))
          do i = 1, NSD
            duidxi(i, j) = sum(ul(:, i)*shgradgu(:, j))
          end do
        end do

        ! velocity on the rotor is mesh velocity, not zero...
        gi = umi
        gphi = 0d0


        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_VEC) > 0) then
          call e3bSTAB_weak(tauB, tauBLS, tauNor, ui - umi, nor, dxidx)

          call e3bRHS_weak(nshl, nor, tauB, tauNor, gw(igauss), &
                           shlu, shgradgu, ui, umi, pri, duidxi, &
                           gi, Rhsu, Rhsp, ti, tmp1, tmp2)

        end if

        if(iand(assemble_tensor_flag, ASSEMBLE_TENSOR_MAT) > 0) then
          call e3bLHS_weak(nshl, ui, umi, duidxi, tauB, tauBLS, tauNor, &
                           gw(igauss), shlu, shgradgu, xKebe11, &
                           xGebe, xDebe1, nor, &
                           xLSebe, xLSUebe, xULSebe, xPLSebe, &
                           Rhsphi, phi, dphidxi)
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
