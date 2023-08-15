!======================================================================
!
!======================================================================
subroutine e3LHS_3D_fluid_Old(nshl, ui, umi, aci, pri, duidxi, &
                              dpridxi, dphidxi, dphidxidxj, dphidti, &
                              rLi, rhoi, mui, tauM, tauP, tauLS, tauC, &
                              tauBar, tauBar1, kap_dc, kap_dc_phi, gwt, shgu, shgradgu, &
                              shhessgu, xKebe11, &
                              xGebe, xDebe1, xMebe, &
                              xLSebe, xLSUebe, &
                              xULSebe, xPLSebe)
  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  real(8), intent(in) :: ui(NSD), umi(NSD), aci(NSD), pri, kap_dc, kap_dc_phi, &
                         duidxi(NSD, NSD), dpridxi(NSD), rLi(NSD), &
                         rhoi, mui, &
                         dphidxi(NSD), dphidxidxj(NSD, NSD), dphidti, &
                         shgu(NSHL), shgradgu(NSHL, NSD), &
                         shhessgu(NSHL, NSD, NSD), &
                         gwt, tauM, tauP, tauLS, tauC, tauBar, tauBar1

  real(8), intent(inout) :: xKebe11(NSD*NSD, NSHL, NSHL), &
                            xGebe(NSD, NSHL, NSHL), &
                            xDebe1(NSD, NSHL, NSHL), &
                            xMebe(NSHL, NSHL), &
                            xLSebe(NSHL, NSHL), &
                            xLSUebe(NSD, NSHL, NSHL), &
                            xULSebe(NSD, NSHL, NSHL), &
                            xPLSebe(NSHL, NSHL)

  integer :: aa, bb, i, j

  real(8) :: fact1, fact2, fact3, &
             tmp1(NSHL), tmp1_ls(NSHL), tmp2(NSHL), tmp4(NSHL, NSHL), &
             advu1(NSD), advu2(NSD), mupkdc, kappadc, &
             res_phic, res_phid, res_phi, nu, tmp3(NSD), tmp5(NSHL), tmp6(NSD), shconvggu(NSHL), shconvggls(NSHL)

  real(8) :: advu_ls(3)
  real(8) :: shconv(NSHL), divu
  real(8) :: drhodphi, dmudphi
  ! loop over local shape functions in each direction
  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt

  drhodphi = rhow - rhoa
  dmudphi = muw - mua

  nu = mui/rhoi

  mupkdc = nu + kap_dc
  kappadc = kappa + kap_dc_phi

!  write(*,*) "kappadc:", kappadc
  tmp1 = 0.0d0
  tmp1_ls = 0.0d0
  tmp2 = 0.0d0
  tmp4 = 0.0d0

  advu1(:) = ui(:) - umi(:)
  advu2(:) = -tauM*rLi(:)
  advu_ls(:) = advu1(:) !+ gravvec(:)*usettle
!  advu_ls (3) = advu_ls(3) - usettle
  do aa = 1, NSHL
    shconv(aa) = sum(shgradgu(aa, :)*advu1)
  enddo
  divu = duidxi(1, 1) + duidxi(2, 2) + duidxi(3, 3)
  tmp1(:) = rhoi*(advu1(1)*shgradgu(:, 1) + &! Na,_j (u_j-um_j)
                  advu1(2)*shgradgu(:, 2) + &
                  advu1(3)*shgradgu(:, 3))

  tmp1_ls(:) = rhoi*(advu_ls(1)*shgradgu(:, 1) + &! Na,_j (u_j-um_j)
                     advu_ls(2)*shgradgu(:, 2) + &
                     advu_ls(3)*shgradgu(:, 3))

  tmp2(:) = advu2(1)*shgradgu(:, 1) + &! Na,_i (-tauM*Li)
            advu2(2)*shgradgu(:, 2) + &
            advu2(3)*shgradgu(:, 3)

  tmp5(:) = kappa*(shhessgu(:, 1, 1) &
                   + shhessgu(:, 2, 2) + shhessgu(:, 3, 3))

  tmp6(:) = gravvec(1)*duidxi(:, 1) + &! gravvec*duidxi
            gravvec(2)*duidxi(:, 2) + &
            gravvec(3)*duidxi(:, 3)

  do aa = 1, NSHL
    !shconvggu(aa) = sum(shgradgu(aa, :)*advu1)
    !shconvggls(aa) = sum(shgradgu(aa, :)*advu_ls)
  end do

  res_phic = dphidti + sum(advu_ls(:)*dphidxi(:)) !Convective part of the residual

  res_phi = res_phic - (kappa*dphidxidxj(1, 1) + &
                        kappa*dphidxidxj(2, 2) + &
                        kappa*dphidxidxj(3, 3))

  do bb = 1, NSHL      ! Diagonal blocks of K11
    do aa = 1, NSHL

      tmp4(aa, bb) = fact1*(shgu(aa)*rhoi*shgu(bb) + &
                            tmp1(aa)*tauM*rhoi*shgu(bb)) + &
                     fact2*(shgu(aa)*tmp1(bb) + &
                            tmp1(aa)*tauM*tmp1(bb) + &
                            mupkdc*(shgradgu(aa, 1)*shgradgu(bb, 1) + &
                                    shgradgu(aa, 2)*shgradgu(bb, 2) + &
                                    shgradgu(aa, 3)*shgradgu(bb, 3)) + &
                            tmp2(aa)*rhoi*tauBar*tmp2(bb))
    end do
  end do

  ! Physics-Physics Interaction
  do bb = 1, NSHL
    do aa = 1, NSHL
      xKebe11(1, aa, bb) = xKebe11(1, aa, bb) + &
                           (tmp4(aa, bb) + &
                            fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 1) + &
                                   shgradgu(aa, 1)*tauC*shgradgu(bb, 1)))*DetJ*gwt !rho*Vlsic*grad(Na)*Rc

      xKebe11(5, aa, bb) = xKebe11(5, aa, bb) + &
                           (tmp4(aa, bb) + &
                            fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 2) + &
                                   shgradgu(aa, 2)*tauC*shgradgu(bb, 2)))*DetJ*gwt

      xKebe11(9, aa, bb) = xKebe11(9, aa, bb) + &
                           (tmp4(aa, bb) + &
                            fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 3) + &
                                   shgradgu(aa, 3)*tauC*shgradgu(bb, 3)))*DetJ*gwt

      xKebe11(2, aa, bb) = xKebe11(2, aa, bb) + &
                           fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 1) + &
                                  shgradgu(aa, 1)*tauC*shgradgu(bb, 2))*DetJ*gwt

      xKebe11(4, aa, bb) = xKebe11(4, aa, bb) + &
                           fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 2) + &
                                  shgradgu(aa, 2)*tauC*shgradgu(bb, 1))*DetJ*gwt

      xKebe11(3, aa, bb) = xKebe11(3, aa, bb) + &
                           fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 1) + &
                                  shgradgu(aa, 1)*tauC*shgradgu(bb, 3))*DetJ*gwt

      xKebe11(7, aa, bb) = xKebe11(7, aa, bb) + &
                           fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 3) + &
                                  shgradgu(aa, 3)*tauC*shgradgu(bb, 1))*DetJ*gwt

      xKebe11(6, aa, bb) = xKebe11(6, aa, bb) + &
                           fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 2) + &
                                  shgradgu(aa, 2)*tauC*shgradgu(bb, 3))*DetJ*gwt

      xKebe11(8, aa, bb) = xKebe11(8, aa, bb) + &
                           fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 3) + &
                                  shgradgu(aa, 3)*tauC*shgradgu(bb, 2))*DetJ*gwt
    end do
  end do

  ! Physics-Mesh
  ! Divergence Matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
                         fact2*(-shgradgu(aa, 1)*shgu(bb) + &
                                tmp1(aa)*tauM*shgradgu(bb, 1))*DetJ*gwt
      xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
                         fact2*(-shgradgu(aa, 2)*shgu(bb) + &
                                tmp1(aa)*tauM*shgradgu(bb, 2))*DetJ*gwt
      xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
                         fact2*(-shgradgu(aa, 3)*shgu(bb) + &
                                tmp1(aa)*tauM*shgradgu(bb, 3))*DetJ*gwt
    end do
  end do

  ! Physics-LS
  ! Divergence Matrix

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(1) - &
  !                                 tauM*tmp1(aa)*gravvec(1)*shgu(bb))*DetJ*gwt
  !     xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(2) - &
  !                                 tauM*tmp1(aa)*gravvec(2)*shgu(bb))*DetJ*gwt
  !     xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(3) - &
  !                                 tauM*tmp1(aa)*gravvec(3)*shgu(bb))*DetJ*gwt
  !   end do
  ! end do

!  do bb = 1, NSHL
!    do aa = 1, NSHL
!      xULSebe(1,aa,bb) = xULSebe(1,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(1)) &
!        tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(1)+ &
!                fact2*shgu(aa)*gravvec(1)*tmp1_ls(bb)) + &
!        )*DetJ*gwt

!      xULSebe(2,aa,bb) = xULSebe(2,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(2)) &
  !tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(2)+ &
  !        fact2*shgu(aa)*gravvec(2)*tmp1_ls(bb))+&
!        )*DetJ*gwt

!      xULSebe(3,aa,bb) = xULSebe(3,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(3)) &
!        tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(3)+ &
!                fact2*shgu(aa)*gravvec(3)*tmp1_ls(bb))
!        )*DetJ*gwt
!    end do
!  end do

  ! Physics-LS
  ! Divergence Matrix, added by Jinhui
  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
  !                          tauBar1*(fact1*shconvggu(aa)*shgu(bb) + fact2*shconvggu(aa)*shconvggls(bb))*DetJ*gwt
  !   end do
  ! end do

  ! LS-Physics
  ! Divergence Matrix

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSUebe(1, aa, bb) = xLSUebe(1, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(1))*DetJ*gwt!+&
! !               shgradgu(aa,1)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !     xLSUebe(2, aa, bb) = xLSUebe(2, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(2))*DetJ*gwt!+&
! !               shgradgu(aa,2)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !     xLSUebe(3, aa, bb) = xLSUebe(3, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(3))*DetJ*gwt!+ &
! !               shgradgu(aa,3)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !   end do
  ! end do

  ! Physics-Physics
  ! Divergence Matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xDebe1(1, aa, bb) = xDebe1(1, aa, bb) + &
                          (fact2*(shgu(aa)*shgradgu(bb, 1) + &
                                  shgradgu(aa, 1)*tauP*tmp1(bb)) + &
                           fact1*(shgradgu(aa, 1)*tauP*rhoi*shgu(bb)))*DetJ*gwt

      xDebe1(2, aa, bb) = xDebe1(2, aa, bb) + &
                          (fact2*(shgu(aa)*shgradgu(bb, 2) + &
                                  shgradgu(aa, 2)*tauP*tmp1(bb)) + &
                           fact1*(shgradgu(aa, 2)*tauP*rhoi*shgu(bb)))*DetJ*gwt

      xDebe1(3, aa, bb) = xDebe1(3, aa, bb) + &
                          (fact2*(shgu(aa)*shgradgu(bb, 3) + &
                                  shgradgu(aa, 3)*tauP*tmp1(bb)) + &
                           fact1*(shgradgu(aa, 3)*tauP*rhoi*shgu(bb)))*DetJ*gwt
    end do
  end do

  ! Mass Matrix and LS matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xMebe(aa, bb) = xMebe(aa, bb) + &
                      fact2*tauP*(shgradgu(aa, 1)*shgradgu(bb, 1) + &
                                  shgradgu(aa, 2)*shgradgu(bb, 2) + &
                                  shgradgu(aa, 3)*shgradgu(bb, 3))*DetJ*gwt

    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      xLSebe(aa, bb) = xLSebe(aa, bb) + &
              (shgu(aa) + tauls * shconv(aa)) * (fact1 * shgu(bb) + fact2 * shconv(bb)) * DetJ * gwt

    end do
  end do
  do aa = 1, NSHL
    do bb = 1, NSHL
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * shgu(aa) * drhodphi * shgu(bb) * rLi(:) / rhoi * DetJ * gwt
      xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * duidxi(:, 1) + duidxi(1, :)) * dmudphi * shgu(bb) * DetJ * gwt
      xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * duidxi(:, 2) + duidxi(2, :)) * dmudphi * shgu(bb) * DetJ * gwt
      xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * duidxi(:, 3) + duidxi(3, :)) * dmudphi * shgu(bb) * DetJ * gwt
      
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * drhodphi * tauM * shconv(aa) * shgu(bb) * rLi(:) * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * rhoi * tauM * shconv(aa) * shgu(bb) * drhodphi * rLi(:) / rhoi * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * drhodphi * shgu(bb) * tauC * divu * shgradgu(aa, :) * DetJ * gwt
    enddo
  enddo
  do aa = 1, NSHL
    do bb = 1, NSHL
      xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
        fact2 * tauls * res_phic * shgu(bb) * shgradgu(aa, :) * DetJ * gwt
      xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
        fact2 * (shgu(aa) + tauls * shconv(aa)) * shgu(bb) * dphidxi(:) * DetJ * gwt
    enddo
  enddo
  do aa = 1, NSHL
    do bb = 1, NSHL
      xPLSebe(aa, bb) = xPLSebe(aa, bb) + &
        fact2 * tauM * drhodphi / rhoi * shgu(bb) * sum(rLi(:) * shgradgu(aa, :)) * DetJ * gwt
    enddo
  enddo


end subroutine e3LHS_3D_fluid_Old

!======================================================================
!
!======================================================================
!subroutine e3LHS_3D_fluid(nshl, ui, umi, aci, pri, duidxi, &
!                          dpridxi, rLi, tauM, tauC, kdc, &
!                          gwt, shlu, shgradlu, shgradgu, shhessgu, &
!                          xKebe11, xGebe, xDebe1, xMebe)
!
!  use aAdjKeep
!  use commonvars
!  implicit none
!
!  integer, intent(in) :: nshl
!
!  integer aa, bb, i, j, k
!
!  real(8) gwt, gwt0, tauM, tauP, tauC, tauBar, vval
!
!  real(8) shlu(NSHL), shgradlu(NSHL, NSD), &
!    shgradgu(NSHL, NSD), shhessgu(NSHL, NSD, NSD)
!
!  real(8) xKebe11(NSD*NSD, NSHL, NSHL), &
!    xKebe21(NSD*NSD, NSHL, NSHL), xKebe22(NSD*NSD, NSHL, NSHL), &
!    xGebe(NSD, NSHL, NSHL), &
!    xDebe1(NSD, NSHL, NSHL), &
!    xMebe(NSHL, NSHL)
!
!  real(8) ui(NSD), umi(NSD), aci(NSD), pri, duidxi(NSD, NSD), &
!    dpridxi(NSD), rLi(NSD)
!
!  real(8) fact1, fact2, fact3, kdc, mupkdc, &
!    tmp1(NSHL), tmp1_ls(NSHL), tmp2(NSHL), tmp3(NSHL), tmp4(NSHL, NSHL), &
!    tmp5(NSHL, NSD), &
!    advu1(NSD), advu2(NSD), &
!    mbulk, mshear, divu, numod, Emod, DetJgwt
!
!  DetJgwt = 1d0*DetJ*gwt !! q-weight
!!...  loop over local shape functions in each direction
!
!  fact1 = almi
!  fact2 = alfi*gami*Delt
!
!  mupkdc = mu + kdc
!
!  advu1(:) = ui(:) - umi(:)
!  advu2(:) = -tauM*rLi(:)
!
!  tmp1(:) = rho*(advu1(1)*shgradgu(:, 1) + & ! Na,_j (u_j-um_j)
!                 advu1(2)*shgradgu(:, 2) + &
!                 advu1(3)*shgradgu(:, 3))
!
!  tmp1(:) = rho*(advu1(1)*shgradgu(:, 1) + & ! Na,_j (u_j-um_j)
!                 advu1(2)*shgradgu(:, 2) + &
!                 advu1(3)*shgradgu(:, 3))
!
!  tmp2(:) = rho*(advu2(1)*shgradgu(:, 1) + & ! Na,_i (-tauM*Li)
!                 advu2(2)*shgradgu(:, 2) + &
!                 advu2(3)*shgradgu(:, 3))
!
!  tmp3(:) = -mu*(shhessgu(:, 1, 1) &
!                 + shhessgu(:, 2, 2) + shhessgu(:, 3, 3))
!
!  do bb = 1, NSHL      ! Diagonal blocks of K11
!    do aa = 1, NSHL
!
!      tmp4(aa, bb) = fact1*(shlu(aa)*rho*shlu(bb) + &
!                            tmp1(aa)*tauM*rho*shlu(bb)) + &
!                     fact2*(shlu(aa)*(tmp1(bb) + tmp2(bb)) + &
!                            (tmp2(aa) + tmp1(aa))*tauM*(tmp1(bb) + tmp3(bb)) + &
!                            mupkdc*(shgradgu(aa, 1)*shgradgu(bb, 1) + &
!                                    shgradgu(aa, 2)*shgradgu(bb, 2) + &
!                                    shgradgu(aa, 3)*shgradgu(bb, 3)))
!
!    end do
!  end do
!
!  do bb = 1, NSHL
!    do aa = 1, NSHL
!
!      xKebe11(1, aa, bb) = xKebe11(1, aa, bb) + &
!                           (tmp4(aa, bb) + &
!                            fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 1) + &
!                                   shgradgu(aa, 1)*tauC*shgradgu(bb, 1)))*DetJ*gwt
!
!      xKebe11(5, aa, bb) = xKebe11(5, aa, bb) + &
!                           (tmp4(aa, bb) + &
!                            fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 2) + &
!                                   shgradgu(aa, 2)*tauC*shgradgu(bb, 2)))*DetJ*gwt
!
!      xKebe11(9, aa, bb) = xKebe11(9, aa, bb) + &
!                           (tmp4(aa, bb) + &
!                            fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 3) + &
!                                   shgradgu(aa, 3)*tauC*shgradgu(bb, 3)))*DetJ*gwt
!
!      xKebe11(2, aa, bb) = xKebe11(2, aa, bb) + &
!                           fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 1) + &
!                                  shgradgu(aa, 1)*tauC*shgradgu(bb, 2))*DetJ*gwt
!
!      xKebe11(4, aa, bb) = xKebe11(4, aa, bb) + &
!                           fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 2) + &
!                                  shgradgu(aa, 2)*tauC*shgradgu(bb, 1))*DetJ*gwt
!
!      xKebe11(3, aa, bb) = xKebe11(3, aa, bb) + &
!                           fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 1) + &
!                                  shgradgu(aa, 1)*tauC*shgradgu(bb, 3))*DetJ*gwt
!
!      xKebe11(7, aa, bb) = xKebe11(7, aa, bb) + &
!                           fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 3) + &
!                                  shgradgu(aa, 3)*tauC*shgradgu(bb, 1))*DetJ*gwt
!
!      xKebe11(6, aa, bb) = xKebe11(6, aa, bb) + &
!                           fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 2) + &
!                                  shgradgu(aa, 2)*tauC*shgradgu(bb, 3))*DetJ*gwt
!
!      xKebe11(8, aa, bb) = xKebe11(8, aa, bb) + &
!                           fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 3) + &
!                                  shgradgu(aa, 3)*tauC*shgradgu(bb, 2))*DetJ*gwt
!
!    end do
!  end do
!
!  do bb = 1, NSHL
!    do aa = 1, NSHL
!
!      xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
!                         fact2*(-shgradgu(aa, 1)*shlu(bb) + &
!                                tmp1(aa)*tauM*shgradgu(bb, 1))*DetJ*gwt
!      xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
!                         fact2*(-shgradgu(aa, 2)*shlu(bb) + &
!                                tmp1(aa)*tauM*shgradgu(bb, 2))*DetJ*gwt
!      xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
!                         fact2*(-shgradgu(aa, 3)*shlu(bb) + &
!                                tmp1(aa)*tauM*shgradgu(bb, 3))*DetJ*gwt
!
!    end do
!  end do
!
!  do bb = 1, NSHL
!    do aa = 1, NSHL
!
!      xDebe1(1, aa, bb) = xDebe1(1, aa, bb) + &
!                          (fact2*(shlu(aa)*shgradgu(bb, 1) + &
!                                  shgradgu(aa, 1)*tauM*(tmp1(bb) + tmp3(bb))) + &
!                           fact1*(shgradgu(aa, 1)*tauM*rho*shlu(bb)))*DetJgwt
!
!      xDebe1(2, aa, bb) = xDebe1(2, aa, bb) + &
!                          (fact2*(shlu(aa)*shgradgu(bb, 2) + &
!                                  shgradgu(aa, 2)*tauM*(tmp1(bb) + tmp3(bb))) + &
!                           fact1*(shgradgu(aa, 2)*tauM*rho*shlu(bb)))*DetJgwt
!
!      xDebe1(3, aa, bb) = xDebe1(3, aa, bb) + &
!                          (fact2*(shlu(aa)*shgradgu(bb, 3) + &
!                                  shgradgu(aa, 3)*tauM*(tmp1(bb) + tmp3(bb))) + &
!                           fact1*(shgradgu(aa, 3)*tauM*rho*shlu(bb)))*DetJgwt
!
!    end do
!  end do
!
!  do bb = 1, NSHL
!    do aa = 1, NSHL
!      xMebe(aa, bb) = xMebe(aa, bb) + &
!                      fact2*tauM*(shgradgu(aa, 1)*shgradgu(bb, 1) + &
!                                  shgradgu(aa, 2)*shgradgu(bb, 2) + &
!                                  shgradgu(aa, 3)*shgradgu(bb, 3))*DetJgwt
!    end do
!  end do
!
!end subroutine e3LHS_3D_fluid

!======================================================================
!
!======================================================================
! subroutine e3LHS_3D_struct(nshl, gwt, shg, shgradg, Ftens, &
!                            Stens, Ctens, xKebe, dc)
! 
!   use aAdjKeep
!   use commonvars
! 
!   implicit none
! 
!   integer, intent(in) :: nshl
! 
!   integer aa, bb, i, j, k, l, m, o
! 
!   real(8) gwt, shg(NSHL), shgradg(NSHL, NSD), &
!     xKebe(NSD*NSD, NSHL, NSHL), fact1, fact2, pt5
! 
!   real(8) Ftens(NSD, NSD), Stens(NSD, NSD), Ctens(NSD, NSD, NSD, NSD)
! 
!   real(8) temp1(NSHL, NSHL), temp2(NSD, NSD, NSD, NSD), &
!     temp3(NSD, NSD, NSD, NSD), temp4(NSHL, NSHL, NSD, NSD), &
!     temp5(NSHL, NSHL, NSD, NSD), temp6(NSHL, NSHL)
! 
!   real(8) dC
! 
!   pt5 = 0.5d+0
! 
! !  dC = 6.0d+5
! !  dC = 5d+4
! !  dC = 0.0d0
! 
!   fact1 = alfi*beti*Delt*Delt
!   fact2 = alfi*gami*Delt
! 
!   temp1 = 0d+0
! 
!   do bb = 1, NSHL
!     do aa = 1, NSHL   ! Mass Contribution
!       temp1(aa, bb) = almi*shg(aa)*rho*shg(bb) &
!                       + fact2*shg(aa)*dC*shg(bb)
!     end do
!   end do
! 
! ! ijlo component
! ! Build material contribution
! 
!   temp4 = 0d+0
! 
!   do j = 1, NSD
!     do i = 1, NSD
!       do bb = 1, NSHL
!         do aa = 1, NSHL
!           temp4(aa, bb, i, j) &
!             = shgradg(aa, 1)*Ctens(i, j, 1, 1)*shgradg(bb, 1) &
!               + shgradg(aa, 2)*Ctens(i, j, 1, 2)*shgradg(bb, 1) &
!               + shgradg(aa, 3)*Ctens(i, j, 1, 3)*shgradg(bb, 1) &
!               + shgradg(aa, 1)*Ctens(i, j, 2, 1)*shgradg(bb, 2) &
!               + shgradg(aa, 2)*Ctens(i, j, 2, 2)*shgradg(bb, 2) &
!               + shgradg(aa, 3)*Ctens(i, j, 2, 3)*shgradg(bb, 2) &
!               + shgradg(aa, 1)*Ctens(i, j, 3, 1)*shgradg(bb, 3) &
!               + shgradg(aa, 2)*Ctens(i, j, 3, 2)*shgradg(bb, 3) &
!               + shgradg(aa, 3)*Ctens(i, j, 3, 3)*shgradg(bb, 3)
!         end do
!       end do
!     end do
!   end do
! 
!   temp4 = temp4*fact1
! 
!   ! Build geometric nonlinearity
! 
!   temp6 = 0d+0
! 
!   do bb = 1, NSHL
!     do aa = 1, NSHL
!       temp6(aa, bb) = shgradg(aa, 1)*Stens(1, 1)*shgradg(bb, 1) &
!                       + shgradg(aa, 2)*Stens(1, 2)*shgradg(bb, 1) &
!                       + shgradg(aa, 3)*Stens(1, 3)*shgradg(bb, 1) &
!                       + shgradg(aa, 1)*Stens(2, 1)*shgradg(bb, 2) &
!                       + shgradg(aa, 2)*Stens(2, 2)*shgradg(bb, 2) &
!                       + shgradg(aa, 3)*Stens(2, 3)*shgradg(bb, 2) &
!                       + shgradg(aa, 1)*Stens(3, 1)*shgradg(bb, 3) &
!                       + shgradg(aa, 2)*Stens(3, 2)*shgradg(bb, 3) &
!                       + shgradg(aa, 3)*Stens(3, 3)*shgradg(bb, 3)
!     end do
!   end do
! 
!   temp6 = temp6*fact1
! 
!   ! loop over elements in each direction
! 
!   do bb = 1, NSHL
!     do aa = 1, NSHL
! 
!       xKebe(1, aa, bb) = xKebe(1, aa, bb) + &
!                          (temp1(aa, bb) + & ! Mass
!                           temp4(aa, bb, 1, 1) + & ! Mat Stiff
!                           temp6(aa, bb))*DetJ*gwt ! Geom Stiff
! 
!       xKebe(5, aa, bb) = xKebe(5, aa, bb) + &
!                          (temp1(aa, bb) + &! Mass
!                           temp4(aa, bb, 2, 2) + &! Mat Stiff
!                           temp6(aa, bb))*DetJ*gwt ! Geom Stiff
! 
!       xKebe(9, aa, bb) = xKebe(9, aa, bb) + &
!                          (temp1(aa, bb) + &! Mass
!                           temp4(aa, bb, 3, 3) + &! Mat Stiff
!                           temp6(aa, bb))*DetJ*gwt ! Geom Stiff
! 
!       xKebe(2, aa, bb) = xKebe(2, aa, bb) + temp4(aa, bb, 1, 2)*DetJ*gwt ! Mat Stiff
!       xKebe(4, aa, bb) = xKebe(4, aa, bb) + temp4(aa, bb, 2, 1)*DetJ*gwt ! Mat Stiff
!       xKebe(3, aa, bb) = xKebe(3, aa, bb) + temp4(aa, bb, 1, 3)*DetJ*gwt ! Mat Stiff
!       xKebe(7, aa, bb) = xKebe(7, aa, bb) + temp4(aa, bb, 3, 1)*DetJ*gwt ! Mat Stiff
!       xKebe(6, aa, bb) = xKebe(6, aa, bb) + temp4(aa, bb, 2, 3)*DetJ*gwt ! Mat Stiff
!       xKebe(8, aa, bb) = xKebe(8, aa, bb) + temp4(aa, bb, 3, 2)*DetJ*gwt ! Mat Stiff
! 
!     end do
!   end do
! 
!   return
! end

!======================================================================
!
!======================================================================
subroutine e3LHS_3D_mesh(nshl, gwt, shgradgu, xKebe22)

  use aAdjKeep
  use commonvars

  implicit none

  integer, intent(in) :: nshl

  integer aa, bb, i, j

  real(8) gwt

  real(8) shgradgu(NSHL, NSD)
  real(8) xKebe22(NSD*NSD, NSHL, NSHL)

  real(8) fact1, fact2, fact3, &
    tmp1(NSHL), tmp2(NSD), &
    tmp3(NSHL), tmp4(NSHL, NSHL), &
    tmp5(NSHL, NSD), &
    advu1(NSD), advu2(NSD), &
    mbulk, mshear, divu, numod, Emod

  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt

  ! Mesh "elastic" parameters
  numod = 3d-1
  Emod = 1d0
  mbulk = numod*Emod/((1d+0 + numod)*(1d+0 - 2d+0*numod))
  mshear = Emod/(2d+0*(1d+0 + numod)) ! PDE

  ! K22 - Mesh-Mesh Interaction - Use Linear Elasticity
  do bb = 1, NSHL
    do aa = 1, NSHL

      xKebe22(1, aa, bb) = xKebe22(1, aa, bb) + &
                           fact3*((mbulk + 2d+0*mshear)*shgradgu(aa, 1)*shgradgu(bb, 1) + &  !St
                                  mshear*shgradgu(aa, 2)*shgradgu(bb, 2) + &
                                  mshear*shgradgu(aa, 3)*shgradgu(bb, 3))*gwt

      xKebe22(5, aa, bb) = xKebe22(5, aa, bb) + &
                           fact3*(mshear*shgradgu(aa, 1)*shgradgu(bb, 1) + &! Stiff
                                  (mbulk + 2d+0*mshear)*shgradgu(aa, 2)*shgradgu(bb, 2) + &
                                  mshear*shgradgu(aa, 3)*shgradgu(bb, 3))*gwt

      xKebe22(9, aa, bb) = xKebe22(9, aa, bb) + &
                           fact3*(mshear*shgradgu(aa, 1)*shgradgu(bb, 1) + &! Stiff
                                  mshear*shgradgu(aa, 2)*shgradgu(bb, 2) + &
                                  (mbulk + 2d+0*mshear)*shgradgu(aa, 3)*shgradgu(bb, 3))*gwt

      xKebe22(2, aa, bb) = xKebe22(2, aa, bb) + &
                           fact3*(mbulk*shgradgu(aa, 1)*shgradgu(bb, 2) + &
                                  mshear*shgradgu(aa, 2)*shgradgu(bb, 1))*gwt

      xKebe22(4, aa, bb) = xKebe22(4, aa, bb) + &
                           fact3*(mshear*shgradgu(aa, 1)*shgradgu(bb, 2) + &
                                  mbulk*shgradgu(aa, 2)*shgradgu(bb, 1))*gwt

      xKebe22(3, aa, bb) = xKebe22(3, aa, bb) + &
                           fact3*(mbulk*shgradgu(aa, 1)*shgradgu(bb, 3) + &
                                  mshear*shgradgu(aa, 3)*shgradgu(bb, 1))*gwt

      xKebe22(7, aa, bb) = xKebe22(7, aa, bb) + &
                           fact3*(mshear*shgradgu(aa, 1)*shgradgu(bb, 3) + &
                                  mbulk*shgradgu(aa, 3)*shgradgu(bb, 1))*gwt

      xKebe22(6, aa, bb) = xKebe22(6, aa, bb) + &
                           fact3*(mbulk*shgradgu(aa, 2)*shgradgu(bb, 3) + &
                                  mshear*shgradgu(aa, 3)*shgradgu(bb, 2))*gwt

      xKebe22(8, aa, bb) = xKebe22(8, aa, bb) + &
                           fact3*(mshear*shgradgu(aa, 2)*shgradgu(bb, 3) + &
                                  mbulk*shgradgu(aa, 3)*shgradgu(bb, 2))*gwt

    end do

  end do

end subroutine e3LHS_3D_mesh

!======================================================================
! LHS for weak BC
!======================================================================
subroutine e3bLHS_weak(nshl, ui, umi, duidxi, rhoi, mui, &
                      tauB, tauNor, gwt, &
                      shlu, shgradgu, xKebe, xGebe, xDebe, nor)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(inout) :: nshl

  real(8), intent(in) :: ui(NSD), umi(NSD), duidxi(NSD, NSD), &
                         rhoi, mui, &
                         tauB, tauNor, gwt, nor(NSD), &
                         shlu(NSHL), shgradgu(NSHL, NSD)
  real(8), intent(inout) :: xKebe(NSD*NSD, NSHL, NSHL), &
                            xGebe(NSD, NSHL, NSHL), &
                            xDebe(NSD, NSHL, NSHL)

  integer :: aa, bb
  real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL, NSHL), &
             unor, umul, munor, gmul, uneg, nu

  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt

  tmp1 = 0.0d0
  tmp2 = 0.0d0

  tmp1(:) = shgradgu(:, 1)*nor(1) + shgradgu(:, 2)*nor(2) &
            + shgradgu(:, 3)*nor(3)

  unor = sum((ui - umi)*nor(:))  ! u \cdot n
  uneg = 0.5d0*(unor - abs(unor))
  munor = tauNor - tauB

  nu = mui/rhoi

  ! gmul =  1d0 => sym
  ! gmul = -1d0 => skew
  gmul = 1.0d0
  do bb = 1, NSHL      ! Diagonal blocks of K
    do aa = 1, NSHL

      tmp2(aa, bb) = -shlu(aa)*mui*tmp1(bb) &
                     - gmul*tmp1(aa)*mui*shlu(bb) &
                     + shlu(aa)*tauB*shlu(bb) &
                     - shlu(aa)*uneg*rhoi*shlu(bb)
    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL

      xKebe(1, aa, bb) = xKebe(1, aa, bb) + &
                         fact2*(tmp2(aa, bb) &
                                - shlu(aa)*mui*shgradgu(bb, 1)*nor(1) &
                                - gmul*shgradgu(aa, 1)*nor(1)*mui*shlu(bb) &
                                + shlu(aa)*nor(1)*munor*nor(1)*shlu(bb))*DetJb*gwt

      xKebe(5, aa, bb) = xKebe(5, aa, bb) + &
                         fact2*(tmp2(aa, bb) &
                                - shlu(aa)*mui*shgradgu(bb, 2)*nor(2) &
                                - gmul*shgradgu(aa, 2)*nor(2)*mui*shlu(bb) &
                                + shlu(aa)*nor(2)*munor*nor(2)*shlu(bb))*DetJb*gwt

      xKebe(9, aa, bb) = xKebe(9, aa, bb) + &
                         fact2*(tmp2(aa, bb) &
                                - shlu(aa)*mui*shgradgu(bb, 3)*nor(3) &
                                - gmul*shgradgu(aa, 3)*nor(3)*mui*shlu(bb) &
                                + shlu(aa)*nor(3)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(2, aa, bb) = xKebe(2, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 1)*nor(2) &
                                - gmul*shgradgu(aa, 2)*nor(1)*mui*shlu(bb) &
                                + shlu(aa)*nor(1)*munor*nor(2)*shlu(bb))*DetJb*gwt

      xKebe(4, aa, bb) = xKebe(4, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 2)*nor(1) &
                                - gmul*shgradgu(aa, 1)*nor(2)*mui*shlu(bb) &
                                + shlu(aa)*nor(2)*munor*nor(1)*shlu(bb))*DetJb*gwt

      xKebe(3, aa, bb) = xKebe(3, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 1)*nor(3) &
                                - gmul*shgradgu(aa, 3)*nor(1)*mui*shlu(bb) &
                                + shlu(aa)*nor(1)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(7, aa, bb) = xKebe(7, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 3)*nor(1) &
                                - gmul*shgradgu(aa, 1)*nor(3)*mui*shlu(bb) &
                                + shlu(aa)*nor(3)*munor*nor(1)*shlu(bb))*DetJb*gwt

      xKebe(6, aa, bb) = xKebe(6, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 2)*nor(3) &
                                - gmul*shgradgu(aa, 3)*nor(2)*mui*shlu(bb) &
                                + shlu(aa)*nor(2)*munor*nor(3)*shlu(bb))*DetJb*gwt

      xKebe(8, aa, bb) = xKebe(8, aa, bb) + &
                         fact2*(-shlu(aa)*mui*shgradgu(bb, 3)*nor(2) &
                                - gmul*shgradgu(aa, 2)*nor(3)*mui*shlu(bb) &
                                + shlu(aa)*nor(3)*munor*nor(2)*shlu(bb))*DetJb*gwt

    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      xDebe(1, aa, bb) = xDebe(1, aa, bb) - &
                         fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xDebe(2, aa, bb) = xDebe(2, aa, bb) - &
                         fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xDebe(3, aa, bb) = xDebe(3, aa, bb) - &
                         fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
                         fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
                         fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
                         fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
    end do
  end do

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSebe(aa, bb) = xLSebe(aa, bb) + &
  !                      fact2*(-shlu(aa)*kappa*tmp1(bb) + &
  !                             shlu(aa)*shlu(bb)*tauBLS - &
  !                             tmp1(aa)*kappa*shlu(bb) - &
  !                             shlu(aa)*uneg*shlu(bb))*DetJb*gwt
  !   end do
  ! end do

  ! if (unor < 0.0d0) then
  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSUebe(1, aa, bb) = xLSUebe(1, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(1)*(phi - gphi))*DetJb*gwt
  !     xLSUebe(2, aa, bb) = xLSUebe(2, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(2)*(phi - gphi))*DetJb*gwt
  !     xLSUebe(3, aa, bb) = xLSUebe(3, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(3)*(phi - gphi))*DetJb*gwt
  !   end do
  ! end do
  ! end if

  ! do aa = 1, NSHL
  !   Rhsphi(aa) = Rhsphi(aa) - &
  !                (shlu(aa)*kappa*sum(dphidxi(:)*nor(:)) + &
  !                 shlu(aa)*tauB*(phi - gphi) + &
  !                 tmp1(aa)*kappa*(phi - gphi) - &
  !                 shlu(aa)*uneg*(phi - gphi))*DetJb*gwt
  ! end do

end subroutine e3bLHS_weak
!======================================================================
! LHS of weak BC for compressible flow
!======================================================================
subroutine e3bLHS_weak_CF(nshl, ui, umi, duidxi, rhoi, mui, &
                          tauB, tauNor, gwt, &
                          shlu, shgradgu, xKebe, xGebe, xDebe, nor)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(inout) :: nshl

  real(8), intent(in) :: ui(NSD), umi(NSD), duidxi(NSD, NSD), &
                         rhoi, mui, &
                         tauB, tauNor, gwt, nor(NSD), &
                         shlu(NSHL), shgradgu(NSHL, NSD)
  real(8), intent(inout) :: xKebe(NSD*NSD, NSHL, NSHL), &
                            xGebe(NSD, NSHL, NSHL), &
                            xDebe(NSD, NSHL, NSHL)

  integer :: aa, bb, ii, jj
  real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL, NSHL), &
             unor, umul, munor, gmul, uneg, nu
  real(8) :: lambda

  ! loop over local shape functions in each direction

  fact1 = almi
  fact2 = alfi*gami*Delt

  tmp1 = 0.0d0
  tmp2 = 0.0d0

  tmp1(:) = shgradgu(:, 1)*nor(1) + shgradgu(:, 2)*nor(2) &
            + shgradgu(:, 3)*nor(3)

  unor = sum((ui - umi)*nor(:))  ! u \cdot n
  uneg = 0.5d0*(unor - abs(unor))
  munor = tauNor - tauB

  nu = mui/rhoi
  lambda = -2d0/3d0 * mui

  ! gmul =  1d0 => sym
  ! gmul = -1d0 => skew
  gmul = 1.0d0
  do bb = 1, NSHL      ! Diagonal blocks of K
    do aa = 1, NSHL

      tmp2(aa, bb) = -shlu(aa)*mui*tmp1(bb) &
                     - gmul*tmp1(aa)*mui*shlu(bb) &
                     + shlu(aa)*tauB*shlu(bb) &
                     - shlu(aa)*uneg*rhoi*shlu(bb)
    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL

      do ii = 1,NSD
        do jj = 1,NSD
          xKebe((ii - 1) * NSD + jj, aa, bb) = xKebe((ii - 1) * NSD + jj, aa, bb) &
              + fact2 * (- shlu(aa)*mui*shgradgu(bb, ii)*nor(jj) &
                         - shlu(aa)*lambda*shgradgu(bb, jj)*nor(ii) &
                         - gmul*shgradgu(aa, jj)*nor(ii)*mui*shlu(bb) &
                         - gmul*shgradgu(aa, ii)*nor(jj)*lambda*shlu(bb) &
                          + shlu(aa)*nor(ii)*munor*nor(jj)*shlu(bb)) * DetJb * gwt
        enddo
        xKebe((ii - 1) * NSD + ii, aa, bb) = xKebe((ii - 1) * NSD + ii, aa, bb) &
              + fact2 * tmp2(aa, bb) * DetJb * gwt
      enddo
      ! xKebe(1, aa, bb) = xKebe(1, aa, bb) + &
      !                    fact2*(tmp2(aa, bb) &
      !                           - shlu(aa)*mui*shgradgu(bb, 1)*nor(1) &
      !                           - gmul*shgradgu(aa, 1)*nor(1)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(1)*munor*nor(1)*shlu(bb))*DetJb*gwt

      ! xKebe(5, aa, bb) = xKebe(5, aa, bb) + &
      !                    fact2*(tmp2(aa, bb) &
      !                           - shlu(aa)*mui*shgradgu(bb, 2)*nor(2) &
      !                           - gmul*shgradgu(aa, 2)*nor(2)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(2)*munor*nor(2)*shlu(bb))*DetJb*gwt

      ! xKebe(9, aa, bb) = xKebe(9, aa, bb) + &
      !                    fact2*(tmp2(aa, bb) &
      !                           - shlu(aa)*mui*shgradgu(bb, 3)*nor(3) &
      !                           - gmul*shgradgu(aa, 3)*nor(3)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(3)*munor*nor(3)*shlu(bb))*DetJb*gwt

      ! xKebe(2, aa, bb) = xKebe(2, aa, bb) + &
      !                    fact2*(- shlu(aa)*mui*shgradgu(bb, 1)*nor(2) &
      !                           - gmul*shgradgu(aa, 2)*nor(1)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(1)*munor*nor(2)*shlu(bb))*DetJb*gwt

      ! xKebe(4, aa, bb) = xKebe(4, aa, bb) + &
      !                    fact2*(-shlu(aa)*mui*shgradgu(bb, 2)*nor(1) &
      !                           - gmul*shgradgu(aa, 1)*nor(2)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(2)*munor*nor(1)*shlu(bb))*DetJb*gwt

      ! xKebe(3, aa, bb) = xKebe(3, aa, bb) + &
      !                    fact2*(-shlu(aa)*mui*shgradgu(bb, 1)*nor(3) &
      !                           - gmul*shgradgu(aa, 3)*nor(1)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(1)*munor*nor(3)*shlu(bb))*DetJb*gwt

      ! xKebe(7, aa, bb) = xKebe(7, aa, bb) + &
      !                    fact2*(-shlu(aa)*mui*shgradgu(bb, 3)*nor(1) &
      !                           - gmul*shgradgu(aa, 1)*nor(3)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(3)*munor*nor(1)*shlu(bb))*DetJb*gwt

      ! xKebe(6, aa, bb) = xKebe(6, aa, bb) + &
      !                    fact2*(-shlu(aa)*mui*shgradgu(bb, 2)*nor(3) &
      !                           - gmul*shgradgu(aa, 3)*nor(2)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(2)*munor*nor(3)*shlu(bb))*DetJb*gwt

      ! xKebe(8, aa, bb) = xKebe(8, aa, bb) + &
      !                    fact2*(-shlu(aa)*mui*shgradgu(bb, 3)*nor(2) &
      !                           - gmul*shgradgu(aa, 2)*nor(3)*mui*shlu(bb) &
      !                           + shlu(aa)*nor(3)*munor*nor(2)*shlu(bb))*DetJb*gwt

    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      ! xDebe(1, aa, bb) = xDebe(1, aa, bb) - &
      !                    fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      ! xDebe(2, aa, bb) = xDebe(2, aa, bb) - &
      !                    fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      ! xDebe(3, aa, bb) = xDebe(3, aa, bb) - &
      !                    fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
      xDebe(:, aa, bb) = xDebe(:, aa, bb) - &
                         fact2*shlu(aa)*shlu(bb)*nor(:)*DetJb*gwt
    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      ! xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
      !                    fact2*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
      ! xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
      !                    fact2*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
      ! xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
      !                    fact2*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
      xGebe(:, aa, bb) = xGebe(:, aa, bb) + &
                         fact2*shlu(aa)*shlu(bb)*nor(:)*DetJb*gwt
    end do
  end do

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSebe(aa, bb) = xLSebe(aa, bb) + &
  !                      fact2*(-shlu(aa)*kappa*tmp1(bb) + &
  !                             shlu(aa)*shlu(bb)*tauBLS - &
  !                             tmp1(aa)*kappa*shlu(bb) - &
  !                             shlu(aa)*uneg*shlu(bb))*DetJb*gwt
  !   end do
  ! end do

  ! if (unor < 0.0d0) then
  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSUebe(1, aa, bb) = xLSUebe(1, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(1)*(phi - gphi))*DetJb*gwt
  !     xLSUebe(2, aa, bb) = xLSUebe(2, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(2)*(phi - gphi))*DetJb*gwt
  !     xLSUebe(3, aa, bb) = xLSUebe(3, aa, bb) - &
  !                          fact2*(shlu(aa)*shlu(bb)*nor(3)*(phi - gphi))*DetJb*gwt
  !   end do
  ! end do
  ! end if

  ! do aa = 1, NSHL
  !   Rhsphi(aa) = Rhsphi(aa) - &
  !                (shlu(aa)*kappa*sum(dphidxi(:)*nor(:)) + &
  !                 shlu(aa)*tauB*(phi - gphi) + &
  !                 tmp1(aa)*kappa*(phi - gphi) - &
  !                 shlu(aa)*uneg*(phi - gphi))*DetJb*gwt
  ! end do

end subroutine e3bLHS_weak_CF

!======================================================================
! LHS for weak BC
!======================================================================
subroutine e3bLHS_outflow(nshl, ui, umi, rhoi, mui, &
                          gwt, shlu, xKebe, nor, phi, gphi, xLSebe, xLSUebe)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(inout) :: nshl

  real(8), intent(in) :: ui(NSD), umi(NSD), gwt, nor(NSD), shlu(NSHL), phi, gphi
  real(8), intent(in) :: rhoi, mui
  real(8), intent(inout) :: xKebe(NSD*NSD, NSHL, NSHL), xLSebe(NSHL, NSHL), xLSUebe(NSD, NSHl, NSHL)

  integer :: aa, bb
  real(8) :: fact2, tmp2(NSHL, NSHL), tmp3(NSHL, NSHL), unor, uneg

  ! loop over local shape functions in each direction

  fact2 = alfi*gami*Delt

  tmp2 = 0.0d0
  tmp3 = 0.0d0

  unor = sum((ui - umi)*nor(:))  ! u \cdot n
  uneg = 0.5d0*(unor - abs(unor))

  do bb = 1, NSHL      ! Diagonal blocks of K
    do aa = 1, NSHL
      tmp2(aa, bb) = -shlu(aa)*uneg*rhoi*shlu(bb)
      tmp3(aa, bb) = -shlu(aa)*uneg*shlu(bb)
    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL

      xKebe(1, aa, bb) = xKebe(1, aa, bb) + &
                         fact2*(tmp2(aa, bb))*DetJb*gwt

      xKebe(5, aa, bb) = xKebe(5, aa, bb) + &
                         fact2*(tmp2(aa, bb))*DetJb*gwt

      xKebe(9, aa, bb) = xKebe(9, aa, bb) + &
                         fact2*(tmp2(aa, bb))*DetJb*gwt

      xLSebe(aa, bb) = xLSebe(aa, bb) + &
                       fact2*(tmp3(aa, bb))*DetJb*gwt

      xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb)
    end do
  end do

end subroutine e3bLHS_outflow

!======================================================================
! LHS for DG
!======================================================================
! subroutine e3bLHS_DG(nshl, shlu, shgradgu, ui, umi, duidxi, tauB, gwt, &
!                      nor, xKebe, xGebe, xDebe)
!   use aAdjKeep
!   use commonvars
!   implicit none
! 
!   integer, intent(in) :: nshl
!   real(8), intent(in) :: shlu(NSHL), shgradgu(NSHL, NSD), &
!                          ui(NSD), umi(NSD), duidxi(NSD, NSD), &
!                          tauB, gwt, nor(NSD)
! 
!   real(8), intent(inout) :: xKebe(NSD*NSD, NSHL, NSHL), &
!                             xGebe(NSD, NSHL, NSHL), &
!                             xDebe(NSD, NSHL, NSHL)
! 
!   integer :: aa, bb
!   real(8) :: fact1, fact2, tmp1(NSHL), tmp2(NSHL, NSHL), &
!              unor, uneg
! 
!   fact1 = almi
!   fact2 = alfi*gami*Delt
! 
!   tmp1 = 0.0d0
!   tmp2 = 0.0d0
! 
!   tmp1(:) = shgradgu(:, 1)*nor(1) + shgradgu(:, 2)*nor(2) &
!             + shgradgu(:, 3)*nor(3)
! 
!   ! Diagonal blocks of K
! 
! !!$  ! No viscous contributions...
! !!$  do bb = 1, NSHL
! !!$    do aa = 1, NSHL
! !!$      tmp2(aa,bb) = shlu(aa)*tauB*shlu(bb)
! !!$    end do
! !!$  end do
! !!$
! !!$  do bb = 1, NSHL
! !!$    do aa = 1, NSHL
! !!$      xKebe(1,aa,bb) = xKebe(1,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
! !!$      xKebe(5,aa,bb) = xKebe(5,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
! !!$      xKebe(9,aa,bb) = xKebe(9,aa,bb) + fact2*(tmp2(aa,bb))*DetJb*gwt
! !!$    end do
! !!$  end do
! 
!   ! with viscous contributions
!   do bb = 1, NSHL
!     do aa = 1, NSHL
!       tmp2(aa, bb) = -0.5d0*(shlu(aa)*mu*tmp1(bb) &
!                              + tmp1(aa)*mu*shlu(bb)) &
!                      + shlu(aa)*tauB*shlu(bb)
!     end do
!   end do
! 
!   do bb = 1, NSHL
!     do aa = 1, NSHL
!       xKebe(1, aa, bb) = xKebe(1, aa, bb) + &
!                          fact2*(tmp2(aa, bb) &
!                                 - 0.5d0*(mu*shgradgu(bb, 1)*nor(1)*shlu(aa) &
!                                          + mu*shgradgu(aa, 1)*nor(1)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(5, aa, bb) = xKebe(5, aa, bb) + &
!                          fact2*(tmp2(aa, bb) &
!                                 - 0.5d0*(mu*shgradgu(bb, 2)*nor(2)*shlu(aa) &
!                                          + mu*shgradgu(aa, 2)*nor(2)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(9, aa, bb) = xKebe(9, aa, bb) + &
!                          fact2*(tmp2(aa, bb) &
!                                 - 0.5d0*(mu*shgradgu(bb, 3)*nor(3)*shlu(aa) &
!                                          + mu*shgradgu(aa, 3)*nor(3)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(2, aa, bb) = xKebe(2, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 1)*nor(2)*shlu(aa) &
!                                         + mu*shgradgu(aa, 2)*nor(1)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(4, aa, bb) = xKebe(4, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 2)*nor(1)*shlu(aa) &
!                                         + mu*shgradgu(aa, 1)*nor(2)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(3, aa, bb) = xKebe(3, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 1)*nor(3)*shlu(aa) &
!                                         + mu*shgradgu(aa, 3)*nor(1)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(7, aa, bb) = xKebe(7, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 3)*nor(1)*shlu(aa) &
!                                         + mu*shgradgu(aa, 1)*nor(3)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(6, aa, bb) = xKebe(6, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 2)*nor(3)*shlu(aa) &
!                                         + mu*shgradgu(aa, 3)*nor(2)*shlu(bb)))*DetJb*gwt
! 
!       xKebe(8, aa, bb) = xKebe(8, aa, bb) + &
!                          fact2*(-0.5d0*(mu*shgradgu(bb, 3)*nor(2)*shlu(aa) &
!                                         + mu*shgradgu(aa, 2)*nor(3)*shlu(bb)))*DetJb*gwt
!     end do
!   end do
! 
!   ! weighting function and pressure
!   do bb = 1, NSHL
!     do aa = 1, NSHL
!       xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
!       xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
!       xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
!     end do
!   end do
! 
!   ! continuity weighting (q) and velocity
!   do bb = 1, NSHL
!     do aa = 1, NSHL
!       xDebe(1, aa, bb) = xDebe(1, aa, bb) - &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(1)*DetJb*gwt
!       xDebe(2, aa, bb) = xDebe(2, aa, bb) - &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(2)*DetJb*gwt
!       xDebe(3, aa, bb) = xDebe(3, aa, bb) - &
!                          fact2*0.5d0*shlu(aa)*shlu(bb)*nor(3)*DetJb*gwt
!     end do
!   end do
! 
! end subroutine e3bLHS_DG

!======================================================================
!
!======================================================================
subroutine e3LHS_3D_fluid_quenching(&
    nshl, ui, umi, aci, pri, duidxi, &
    dpridxi, dphidxi, dphidxidxj, dphidti, phii, &
    rLi, fi, rhoi, mui, Ti, &
    tauM, tauP, tauLS, tauC, &
    tauBar, tauBar1, kap_dc, kap_dc_phi, gwt, shgu, shgradgu, &
    shhessgu, xKebe11, &
    xGebe, xDebe1, xMebe, &
    xLSebe, xLSUebe, &
    xULSebe, xPLSebe)

  use aAdjKeep
  use commonvars
  implicit none

  integer, intent(in) :: nshl

  real(8), intent(in) :: ui(NSD), umi(NSD), aci(NSD), pri, kap_dc, kap_dc_phi, &
                         duidxi(NSD, NSD), dpridxi(NSD), rLi(NSD), fi(NSD), &
                         dphidxi(NSD), dphidxidxj(NSD, NSD), dphidti, phii, &
                         rhoi, mui, Ti,&
                         shgu(NSHL), shgradgu(NSHL, NSD), &
                         shhessgu(NSHL, NSD, NSD), &
                         gwt, tauM, tauP, tauLS, tauC, tauBar, tauBar1

  real(8), intent(inout) :: xKebe11(NSD*NSD, NSHL, NSHL), &
                            xGebe(NSD, NSHL, NSHL), &
                            xDebe1(NSD, NSHL, NSHL), &
                            xMebe(NSHL, NSHL), &
                            xLSebe(NSHL, NSHL), &
                            xLSUebe(NSD, NSHL, NSHL), &
                            xULSebe(NSD, NSHL, NSHL), &
                            xPLSebe(NSHL, NSHL)

  integer :: aa, bb, i, j

  real(8) :: fact1, fact2, fact3, &
             advci(NSD), advfi(NSD), uprime(NSD), mupkdc, kappadc, &
             res_phic, res_phid, res_phi, nu, tmp3(NSD), tmp5(NSHL), tmp6(NSD), shconvggu(NSHL), shconvggls(NSHL)

  ! real(8) :: advu_ls(3)
  real(8) :: shconv(NSHL), shconv_full(NSHL), divu, drmdu(NSHL)
  real(8) :: drhodphi, dmudphi, tmp, tmp1, mdot, vdot, rm(NSD), dmdphii, dmdTi
  real(8) :: diag(NSHL, NSHL), phic
  ! loop over local shape functions in each direction
  fact1 = almi
  fact2 = alfi*gami*Delt
  fact3 = alfi*beti*Delt*Delt

  drhodphi = rhow - rhoa
  dmudphi = muw - mua

  nu = mui/rhoi

  mupkdc = mui + kap_dc
  kappadc = kappa + kap_dc_phi

  phic = max(min(phii, 1.0d0), 0.0d0)
  if(Ti > Ts) then
    mdot = c_evap * (1-phic) * rhow * (Ti - Ts) / Ts
    dmdphii = -c_evap * rhow * (Ti - Ts) / Ts
    dmdTi = c_evap * (1-phic) * rhow / Ts
  else
    mdot = c_cond * (phic) * rhoa * (Ti - Ts) / Ts
    dmdphii = c_cond * rhoa * (Ti - Ts) / Ts
    dmdTi = c_cond * (phic) * rhoa / Ts
  endif
  vdot = mdot / rhoa - mdot / rhow
!  write(*,*) "kappadc:", kappadc
  ! tmp1_ls = 0.0d0
  ! tmp2 = 0.0d0
  !tmp4 = 0.0d0

  advci(:) = ui(:) - umi(:)
  uprime(:) = -tauM*rLi(:)
  advfi(:) = advci(:) !+ uprime(:)

  do aa = 1, NSHL
    !shconv_full(aa) = sum(shgradgu(aa, :)*advfi(:))
    shconv(aa) = sum(shgradgu(aa, :)*advci(:))
  enddo
  shconv_full(:) = shconv(:)
  divu = duidxi(1, 1) + duidxi(2, 2) + duidxi(3, 3)

  ! tmp1_ls(:) = rhoi*(advu_ls(1)*shgradgu(:, 1) + &! Na,_j (u_j-um_j)
  !                   advu_ls(2)*shgradgu(:, 2) + &
  !                   advu_ls(3)*shgradgu(:, 3))

  tmp5(:) = kappa*(shhessgu(:, 1, 1) &
                   + shhessgu(:, 2, 2) + shhessgu(:, 3, 3))

  tmp6(:) = gravvec(1)*duidxi(:, 1) + &! gravvec*duidxi
            gravvec(2)*duidxi(:, 2) + &
            gravvec(3)*duidxi(:, 3)
  drmdu(:) = rhoi * (fact1*shgu(:)+fact2*shconv(:))
  do aa = 1, NSHL
    !shconvggu(aa) = sum(shgradgu(aa, :)*advu1)
    !shconvggls(aa) = sum(shgradgu(aa, :)*advu_ls)
  end do

  res_phic = dphidti + sum(advci(:)*dphidxi(:)) + phii * vdot - mdot / rhoa 
  ! res_phic = dphidti + sum(advci(:)*dphidxi(:)) - mdot / rhoa 

  ! res_phic = res_phic - (kappa*dphidxidxj(1, 1) + &
  !                       kappa*dphidxidxj(2, 2) + &
  !                       kappa*dphidxidxj(3, 3))

  tmp1 = rhoi * vdot + sum(advci(:) * dphidxi(:)) * drhodphi
  do bb = 1, NSHL      ! Diagonal blocks of K11
    do aa = 1, NSHL

      ! tmp4(aa, bb) = fact1*(shgu(aa)*rhoi*shgu(bb) + &
      !                       shconv(aa)*tauM*rhoi*rhoi*shgu(bb)) + &
      !                fact2*(shgu(aa)*rhoi*shconv(bb) + &
      !                       rhoi*shconv(aa)*tauM*rhoi*shconv(bb) + &
      !                       mupkdc * sum(shgradgu(aa, :)*shgradgu(bb, :)) + &
      !                       uprimegradN(aa)*rhoi*tauBar*uprimegradN(bb)) + &
      !                fact1*tauM*shgu(aa)*rhoi*shgu(bb)*tmp + &
      !                fact2*tauM*shgu(aa)*rhoi*shconv(bb) * tmp
      tmp = 0d0
      tmp = tmp + fact1 * shgu(aa) * rhoi * shgu(bb) + fact2 * shgu(aa) * rhoi * shconv_full(bb)
      tmp = tmp + fact2 * mupkdc * sum(shgradgu(aa, :)*shgradgu(bb, :))
      tmp = tmp + rhoi * tauM * shconv_full(aa) * drmdu(bb)
      tmp = tmp + shgu(aa) * tauM * drmdu(bb) * tmp1
      diag(aa, bb) = tmp
    end do
  end do

  ! Physics-Physics Interaction
  do bb = 1, NSHL
    do aa = 1, NSHL
      do i = 1, NSD
        do j = 1, NSD
          xKebe11((i - 1) * NSD + j, aa, bb) = xKebe11((i - 1) * NSD + j, aa, bb) + &
                (fact2*(shgradgu(aa, j)*mupkdc*shgradgu(bb, i) + &
                       shgradgu(aa, i)*(tauC-mui*2d0/3d0)*shgradgu(bb, j)))*DetJ*gwt !rho*Vlsic*grad(Na)*Rc
        end do
        xKebe11((i-1) * NSD + i, aa, bb) = xKebe11((i-1) * NSD + i, aa, bb) + &
                                          diag(aa, bb)*DetJ*gwt !rho*Vlsic*grad(Na)*Rc
      end do 
      !xKebe11(1, aa, bb) = xKebe11(1, aa, bb) + &
      !                     (tmp4(aa, bb) + &
      !                      fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 1) + &
      !                             shgradgu(aa, 1)*(tauC+mui/3d0)*shgradgu(bb, 1)))*DetJ*gwt !rho*Vlsic*grad(Na)*Rc

      !xKebe11(5, aa, bb) = xKebe11(5, aa, bb) + &
      !                     (tmp4(aa, bb) + &
      !                      fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 2) + &
      !                             shgradgu(aa, 2)*(tauC+mui/3d0)*shgradgu(bb, 2)))*DetJ*gwt

      !xKebe11(9, aa, bb) = xKebe11(9, aa, bb) + &
      !                     (tmp4(aa, bb) + &
      !                      fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 3) + &
      !                             shgradgu(aa, 3)*(tauC+mui/3d0)*shgradgu(bb, 3)))*DetJ*gwt

      !xKebe11(2, aa, bb) = xKebe11(2, aa, bb) + &
      !                     fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 1) + &
      !                            shgradgu(aa, 1)*(tauC+mui/3d0)*shgradgu(bb, 2))*DetJ*gwt

      !xKebe11(4, aa, bb) = xKebe11(4, aa, bb) + &
      !                     fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 2) + &
      !                            shgradgu(aa, 2)*(tauC+mui/3d0)*shgradgu(bb, 1))*DetJ*gwt

      !xKebe11(3, aa, bb) = xKebe11(3, aa, bb) + &
      !                     fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 1) + &
      !                            shgradgu(aa, 1)*(tauC+mui/3d0)*shgradgu(bb, 3))*DetJ*gwt

      !xKebe11(7, aa, bb) = xKebe11(7, aa, bb) + &
      !                     fact2*(shgradgu(aa, 1)*mupkdc*shgradgu(bb, 3) + &
      !                            shgradgu(aa, 3)*(tauC+mui/3d0)*shgradgu(bb, 1))*DetJ*gwt

      !xKebe11(6, aa, bb) = xKebe11(6, aa, bb) + &
      !                     fact2*(shgradgu(aa, 3)*mupkdc*shgradgu(bb, 2) + &
      !                            shgradgu(aa, 2)*(tauC+mui/3d0)*shgradgu(bb, 3))*DetJ*gwt

      !xKebe11(8, aa, bb) = xKebe11(8, aa, bb) + &
      !                     fact2*(shgradgu(aa, 2)*mupkdc*shgradgu(bb, 3) + &
      !                            shgradgu(aa, 3)*(tauC+mui/3d0)*shgradgu(bb, 2))*DetJ*gwt
    end do
  end do

  ! Physics-Mesh
  ! Divergence Matrix
  ! dRm/dP
  do bb = 1, NSHL
    do aa = 1, NSHL
      xGebe(:, aa, bb) = xGebe(:, aa, bb) &
                       - fact2 * shgradgu(aa, :)*shgu(bb)*DetJ*gwt &
                       + fact2 * rhoi*shconv_full(aa)*tauM*shgradgu(bb, :)*DetJ*gwt &
                       + fact2 * shgu(aa) * tauM * shgradgu(bb, :)*tmp1*DetJ*gwt
      
      ! xGebe(1, aa, bb) = xGebe(1, aa, bb) + &
      !                    fact2*(-shgradgu(aa, 1)*shgu(bb) + &
      !                           rhoi*ugradN(aa)*tauM*shgradgu(bb, 1))*DetJ*gwt
      ! xGebe(2, aa, bb) = xGebe(2, aa, bb) + &
      !                    fact2*(-shgradgu(aa, 2)*shgu(bb) + &
      !                           ugradN(aa)*tauM*shgradgu(bb, 2))*DetJ*gwt
      ! xGebe(3, aa, bb) = xGebe(3, aa, bb) + &
      !                    fact2*(-shgradgu(aa, 3)*shgu(bb) + &
      !                           ugradN(aa)*tauM*shgradgu(bb, 3))*DetJ*gwt
    end do
  end do

  ! Physics-LS
  ! Divergence Matrix

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(1) - &
  !                                 tauM*tmp1(aa)*gravvec(1)*shgu(bb))*DetJ*gwt
  !     xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(2) - &
  !                                 tauM*tmp1(aa)*gravvec(2)*shgu(bb))*DetJ*gwt
  !     xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
  !                          fact2*(-shgu(aa)*shgu(bb)*gravvec(3) - &
  !                                 tauM*tmp1(aa)*gravvec(3)*shgu(bb))*DetJ*gwt
  !   end do
  ! end do

!  do bb = 1, NSHL
!    do aa = 1, NSHL
!      xULSebe(1,aa,bb) = xULSebe(1,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(1)) &
!        tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(1)+ &
!                fact2*shgu(aa)*gravvec(1)*tmp1_ls(bb)) + &
!        )*DetJ*gwt

!      xULSebe(2,aa,bb) = xULSebe(2,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(2)) &
  !tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(2)+ &
  !        fact2*shgu(aa)*gravvec(2)*tmp1_ls(bb))+&
!        )*DetJ*gwt

!      xULSebe(3,aa,bb) = xULSebe(3,aa,bb) + &
!       (fact2*(-shgu(aa)*shgu(bb)*gravvec(3)) &
!        tauLS*( fact1*shgu(aa)*shgu(bb)*gravvec(3)+ &
!                fact2*shgu(aa)*gravvec(3)*tmp1_ls(bb))
!        )*DetJ*gwt
!    end do
!  end do

  ! Physics-LS
  ! Divergence Matrix, added by Jinhui
  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
  !                          tauBar1*(fact1*shconvggu(aa)*shgu(bb) + fact2*shconvggu(aa)*shconvggls(bb))*DetJ*gwt
  !   end do
  ! end do

  ! LS-Physics
  ! Divergence Matrix

  ! do bb = 1, NSHL
  !   do aa = 1, NSHL
  !     xLSUebe(1, aa, bb) = xLSUebe(1, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(1))*DetJ*gwt!+&
! !               shgradgu(aa,1)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !     xLSUebe(2, aa, bb) = xLSUebe(2, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(2))*DetJ*gwt!+&
! !               shgradgu(aa,2)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !     xLSUebe(3, aa, bb) = xLSUebe(3, aa, bb) + &
  !                          fact2*(shgu(aa)*shgu(bb)*dphidxi(3))*DetJ*gwt!+ &
! !               shgradgu(aa,3)*tauLS*shgu(bb)*res_phi)*DetJ*gwt
  !   end do
  ! end do

  ! Physics-Physics
  ! Divergence Matrix
  ! dRc/dU
  do bb = 1, NSHL
    do aa = 1, NSHL
      xDebe1(:, aa, bb) = xDebe1(:, aa, bb) + &
                          (fact2*(shgu(aa)*shgradgu(bb, :) + &
                                  rhoi*shgradgu(aa, :)*tauM*shconv(bb)) + &
                           fact1*(shgradgu(aa, :)*tauM*rhoi*shgu(bb)))*DetJ*gwt

      ! xDebe1(1, aa, bb) = xDebe1(1, aa, bb) + &
      !                     (fact2*(shgu(aa)*shgradgu(bb, 1) + &
      !                             shgradgu(aa, 1)*tauP*ugradN(bb)) + &
      !                      fact1*(shgradgu(aa, 1)*tauP*rhoi*shgu(bb)))*DetJ*gwt

      ! xDebe1(2, aa, bb) = xDebe1(2, aa, bb) + &
      !                     (fact2*(shgu(aa)*shgradgu(bb, 2) + &
      !                             shgradgu(aa, 2)*tauP*ugradN(bb)) + &
      !                      fact1*(shgradgu(aa, 2)*tauP*rhoi*shgu(bb)))*DetJ*gwt

      ! xDebe1(3, aa, bb) = xDebe1(3, aa, bb) + &
      !                     (fact2*(shgu(aa)*shgradgu(bb, 3) + &
      !                             shgradgu(aa, 3)*tauP*ugradN(bb)) + &
      !                      fact1*(shgradgu(aa, 3)*tauP*rhoi*shgu(bb)))*DetJ*gwt
    end do
  end do

  ! Mass Matrix and LS matrix
  do bb = 1, NSHL
    do aa = 1, NSHL
      xMebe(aa, bb) = xMebe(aa, bb) + &
                      fact2 * tauP * sum(shgradgu(aa, :)*shgradgu(bb, :)) * DetJ * gwt

    end do
  end do

  do bb = 1, NSHL
    do aa = 1, NSHL
      tmp = 0d0
      tmp = tmp + fact1 * shgu(bb) + fact2 * shconv(bb)
      tmp = tmp + fact2 * shgu(bb) * vdot 
      tmp = tmp + fact2 * shgu(bb) * dmdphii * phii * (1.0d0/rhoa - 1.0d0/rhow)
      tmp = tmp - fact2 * shgu(bb) * dmdphii / rhoa
      xLSebe(aa, bb) = (shgu(aa) + tauls * shconv(aa)) * tmp * DetJ * gwt
      xLSebe(aa, bb) = xLSebe(aa, bb) + &
              fact2 * kappadc * sum(shgradgu(aa, :)*shgradgu(bb, :)) * DetJ * gwt

    end do
  end do
  rm(1) = rhoi * (aci(1) + sum(advci(:) * duidxi(:, 1))) - fi(1)
  rm(2) = rhoi * (aci(2) + sum(advci(:) * duidxi(:, 2))) - fi(2)
  rm(3) = rhoi * (aci(3) + sum(advci(:) * duidxi(:, 3))) - fi(3)
  do aa = 1, NSHL
    do bb = 1, NSHL
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * shgu(aa) * drhodphi * shgu(bb) * rm(:) / rhoi * DetJ * gwt
      xULSebe(1, aa, bb) = xULSebe(1, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * (duidxi(:, 1) + duidxi(1, :))) * dmudphi * shgu(bb) * DetJ * gwt
      xULSebe(2, aa, bb) = xULSebe(2, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * (duidxi(:, 2) + duidxi(2, :))) * dmudphi * shgu(bb) * DetJ * gwt
      xULSebe(3, aa, bb) = xULSebe(3, aa, bb) + &
        fact2 * sum(shgradgu(aa, :) * (duidxi(:, 3) + duidxi(3, :))) * dmudphi * shgu(bb) * DetJ * gwt
      
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * drhodphi * tauM * shconv_full(aa) * shgu(bb) * rLi(:) * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * rhoi * tauM * shconv_full(aa) * shgu(bb) * drhodphi * rm(:) / rhoi * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * shgu(bb) * (tauC*drhodphi/rhoi - dmudphi*2d0/3d0) * divu * shgradgu(aa, :) * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * drhodphi * rm(:)/rhoi * tauM * shgu(aa) * shgu(bb) * tmp1 * DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * drhodphi * tauM * shgu(aa) * rLi(:) * (shgu(bb)*vdot + shconv_full(bb))* DetJ * gwt
      xULSebe(:, aa, bb) = xULSebe(:, aa, bb) + &
        fact2 * shgu(aa) * tauM * rLi(:) * rhoi * (1/rhoa - 1/rhow) * dmdphii * shgu(bb) * DetJ * gwt
    enddo
  enddo
  do aa = 1, NSHL
    do bb = 1, NSHL
      xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
        fact2 * tauls * res_phic * shgu(bb) * shgradgu(aa, :) * DetJ * gwt
      xLSUebe(:, aa, bb) = xLSUebe(:, aa, bb) + &
        fact2 * (shgu(aa) + tauls * shconv(aa)) * shgu(bb) * dphidxi(:) * DetJ * gwt
    enddo
  enddo
  do aa = 1, NSHL
    do bb = 1, NSHL
      xPLSebe(aa, bb) = xPLSebe(aa, bb) + &
        fact2 * tauM * drhodphi / rhoi * shgu(bb) * sum(rm(:) * shgradgu(aa, :)) * DetJ * gwt + &
        shgu(aa) * fact2 * dmdphii * (1/rhoa - 1/rhow) * shgu(bb) * DetJ * gwt
    enddo
  enddo
  !xULSebe = 0
  !xLSUebe = 0
  !xPLSebe = 0

end subroutine e3LHS_3D_fluid_quenching
