!======================================================================
!
!======================================================================
subroutine genGPandGW(gp, gw, ngp)

  use commonvars
  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 3), gw(ngp)

  if (iga) then
    call genGPandGWcube(gp, gw, ngp)
  else
    call genGPandGWtetpri(gp, gw, ngp)
  end if

end subroutine genGPandGW

!======================================================================
!
!======================================================================
subroutine genGPandGWb(gp, gw, ngp)

  use commonvars
  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 2), gw(ngp)

  if (iga) then
    call genGPandGWquad(gp, gw, ngp)
  else
    call genGPandGWtri(gp, gw, ngp)
  end if

end subroutine genGPandGWb

!======================================================================
! Generate surface gauss points mapped to volume
!======================================================================
subroutine genGPMap(nshl, ngp, for, iga, gp)
  implicit none

  integer, intent(in)  :: nshl, ngp, for
  real(8), intent(out) :: gp(ngp, 3)
  logical, intent(in)  :: iga

  if (iga) then
    gp = 0.0d0
  else
    if (nshl == 4) then
      call genGPtritet(ngp, for, gp)
    else if (nshl == 6) then
      call genGPtripri(ngp, for, gp)
    else
      write (*, *) "ERROR: undefined nshl in genGPMap"
      stop
    end if
  end if

end subroutine genGPMap

!======================================================================
!
!======================================================================
subroutine genGPandGWcube(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 3), gw(ngp)

  integer :: ngp1, i, j, k, ip
  real(8), allocatable :: gp1(:), gw1(:)

  if (ngp == 8) then
    ngp1 = 2
  else if (ngp == 27) then
    ngp1 = 3
  else if (ngp == 64) then
    ngp1 = 4
  end if

  allocate (gp1(ngp1), gw1(ngp1))
  call genGPandGWline(gp1, gw1, ngp1)

  ip = 0
  do k = 1, ngp1
    do j = 1, ngp1
      do i = 1, ngp1
        ip = ip + 1
        gp(ip, 1) = gp1(i)
        gp(ip, 2) = gp1(j)
        gp(ip, 3) = gp1(k)
        gw(ip) = gw1(i)*gw1(j)*gw1(k)
      end do
    end do
  end do

  deallocate (gp1, gw1)

end subroutine genGPandGWcube

!======================================================================
!
!======================================================================
subroutine genGPandGWquad(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 2), gw(ngp)

  integer :: ngp1, i, j, ip
  real*8, allocatable :: gp1(:), gw1(:)

  if (ngp == 4) then
    ngp1 = 2
  else if (ngp == 9) then
    ngp1 = 3
  else if (ngp == 16) then
    ngp1 = 4
  end if

  allocate (gp1(ngp1), gw1(ngp1))
  call genGPandGWline(gp1, gw1, ngp1)

  ip = 0
  do j = 1, ngp1
    do i = 1, ngp1
      ip = ip + 1
      gp(ip, 1) = gp1(i)
      gp(ip, 2) = gp1(j)
      gw(ip) = gw1(i)*gw1(j)
    end do
  end do

  deallocate (gp1, gw1)

end subroutine genGPandGWquad

!======================================================================
! Subroutine generates Gaussian Quadrature points and weights
!======================================================================
subroutine genGPandGWline(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp), gw(ngp)

  if (ngp == 1) then
    gp(1) = 0.0d0
    gw(1) = 2.0d0
  else if (ngp == 2) then
    gp(1) = -0.5773502691896257645091488d0
    gp(2) = 0.5773502691896257645091488d0
    gw(1) = 1.0d0
    gw(2) = 1.0d0
  else if (ngp == 3) then
    gp(1) = -0.77459666924148337703585308d0
    gp(2) = 0.0d0
    gp(3) = 0.77459666924148337703585308d0
    gw(1) = 0.55555555555555555555555556d0
    gw(2) = 0.88888888888888888888888889d0
    gw(3) = 0.55555555555555555555555556d0
  else if (ngp == 4) then
    gp(1) = -0.86113631159405257524d0
    gp(2) = -0.33998104358485626481d0
    gp(3) = 0.33998104358485626481d0
    gp(4) = 0.86113631159405257524d0
    gw(1) = 0.34785484513745385736d0
    gw(2) = 0.65214515486254614264d0
    gw(3) = 0.65214515486254614264d0
    gw(4) = 0.34785484513745385736d0
  else if (ngp == 5) then
    gp(1) = -0.90617984593866399282d0
    gp(2) = -0.53846931010568309105d0
    gp(3) = 0.0d0
    gp(4) = 0.53846931010568309105d0
    gp(5) = 0.90617984593866399282d0
    gw(1) = 0.23692688505618908749d0
    gw(2) = 0.47862867049936646808d0
    gw(3) = 0.56888888888888888888d0
    gw(4) = 0.47862867049936646808d0
    gw(5) = 0.23692688505618908749d0
  else if (ngp == 6) then
    gp(1) = -0.9324695142031520d0
    gp(2) = -0.6612093864662645d0
    gp(3) = -0.2386191860831969d0
    gp(4) = 0.2386191860831969d0
    gp(5) = 0.6612093864662645d0
    gp(6) = 0.9324695142031520d0
    gw(1) = 0.1713244923791703d0
    gw(2) = 0.3607615730481386d0
    gw(3) = 0.4679139345726911d0
    gw(4) = 0.4679139345726911d0
    gw(5) = 0.3607615730481386d0
    gw(6) = 0.1713244923791703d0
  else if (ngp == 7) then
    gp(1) = -0.9491079123427585d0
    gp(2) = -0.7415311855993944d0
    gp(3) = -0.4058451513773972d0
    gp(4) = 0.0d0
    gp(5) = 0.4058451513773972d0
    gp(6) = 0.7415311855993944d0
    gp(7) = 0.9491079123427585d0
    gw(1) = 0.1294849661688697d0
    gw(2) = 0.2797053914892767d0
    gw(3) = 0.3818300505051189d0
    gw(4) = 0.4179591836734694d0
    gw(5) = 0.3818300505051189d0
    gw(6) = 0.2797053914892767d0
    gw(7) = 0.1294849661688697d0
  else if (ngp == 8) then
    gp(1) = -.9602898564975362d0
    gp(2) = -.7966664774136267d0
    gp(3) = -.5255324099163290d0
    gp(4) = -.1834346424956498d0
    gp(5) = .1834346424956498d0
    gp(6) = .5255324099163290d0
    gp(7) = .7966664774136267d0
    gp(8) = .9602898564975362d0
    gw(1) = .1012285362903763d0
    gw(2) = .2223810344533745d0
    gw(3) = .3137066458778873d0
    gw(4) = .3626837833783620d0
    gw(5) = .3626837833783620d0
    gw(6) = .3137066458778873d0
    gw(7) = .2223810344533745d0
    gw(8) = .1012285362903763d0
  else if (ngp == 9) then
    gp(1) = -.9681602395076261d0
    gp(2) = -.8360311073266358d0
    gp(3) = -.6133714327005904d0
    gp(4) = -.3242534234038089d0
    gp(5) = 0.0d0
    gp(6) = .3242534234038089d0
    gp(7) = .6133714327005904d0
    gp(8) = .8360311073266358d0
    gp(9) = .9681602395076261d0
    gw(1) = .0812743883615744d0
    gw(2) = .1806481606948574d0
    gw(3) = .2606106964029354d0
    gw(4) = .3123470770400029d0
    gw(5) = .3302393550012598d0
    gw(6) = .3123470770400028d0
    gw(7) = .2606106964029355d0
    gw(8) = .1806481606948574d0
    gw(9) = .0812743883615744d0
  else if (ngp == 10) then
    gp(1) = -0.97390653d0
    gp(2) = -0.86506337d0
    gp(3) = -0.67940957d0
    gp(4) = -0.43339539d0
    gp(5) = -0.14887434d0
    gp(6) = 0.14887434d0
    gp(7) = 0.43339539d0
    gp(8) = 0.67940957d0
    gp(9) = 0.86506337d0
    gp(10) = 0.97390653d0
    gw(1) = 0.06667134d0
    gw(2) = 0.14945135d0
    gw(3) = 0.21908636d0
    gw(4) = 0.26926672d0
    gw(5) = 0.29552422d0
    gw(6) = 0.29552422d0
    gw(7) = 0.26926672d0
    gw(8) = 0.21908636d0
    gw(9) = 0.14945135d0
    gw(10) = 0.06667134d0
  else if (ngp == 15) then
    gp(1) = -0.9879925180204854d0
    gp(2) = -0.9372733924007059d0
    gp(3) = -0.8482065834104272d0
    gp(4) = -0.7244177313601700d0
    gp(5) = -0.5709721726085388d0
    gp(6) = -0.3941513470775634d0
    gp(7) = -0.2011940939974345d0
    gp(8) = 0.0d0
    gp(9) = 0.2011940939974345d0
    gp(10) = 0.3941513470775634d0
    gp(11) = 0.5709721726085388d0
    gp(12) = 0.7244177313601700d0
    gp(13) = 0.8482065834104272d0
    gp(14) = 0.9372733924007059d0
    gp(15) = 0.9879925180204854d0
    gw(1) = 0.03075324199611807d0
    gw(2) = 0.07036604748811134d0
    gw(3) = 0.1071592204671351d0
    gw(4) = 0.1395706779261761d0
    gw(5) = 0.1662692058169852d0
    gw(6) = 0.1861610000155741d0
    gw(7) = 0.1984314853271374d0
    gw(8) = 0.2025782419255562d0
    gw(9) = 0.1984314853271374d0
    gw(10) = 0.1861610000155741d0
    gw(11) = 0.1662692058169852d0
    gw(12) = 0.1395706779261761d0
    gw(13) = 0.1071592204671351d0
    gw(14) = 0.07036604748811134d0
    gw(15) = 0.03075324199611807d0
  else
    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPandGWline", ngp
    stop
  end if

end subroutine genGPandGWline

!======================================================================
! Subroutine generates Gaussian Quadrature points and weights
! tet - ngp = 1/4/8
! pri - ngp = 6
!======================================================================
subroutine genGPandGWtetpri(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 3), gw(ngp)

  if (ngp == 1) then

    gp(1, 1) = 0.25d0
    gp(1, 2) = 0.25d0
    gp(1, 3) = 0.25d0

    gw(1) = 1d0/6d0

  else if (ngp == 4) then   ! tet

    gp(1, 1) = 0.5854101966249685d0
    gp(1, 2) = 0.1381966011250105d0
    gp(1, 3) = 0.1381966011250105d0

    gp(2, 1) = 0.1381966011250105d0
    gp(2, 2) = 0.5854101966249685d0
    gp(2, 3) = 0.1381966011250105d0

    gp(3, 1) = 0.1381966011250105d0
    gp(3, 2) = 0.1381966011250105d0
    gp(3, 3) = 0.5854101966249685d0

    gp(4, 1) = 0.1381966011250105d0
    gp(4, 2) = 0.1381966011250105d0
    gp(4, 3) = 0.1381966011250105d0

    gw(1) = 0.2500000000000000d0/6d0
    gw(2) = 0.2500000000000000d0/6d0
    gw(3) = 0.2500000000000000d0/6d0
    gw(4) = 0.2500000000000000d0/6d0

  else if (ngp == 8) then

    gp(1, 1) = -0.57735026918963d0
    gp(1, 2) = -0.57735026918963d0
    gp(1, 3) = -0.57735026918963d0

    gp(2, 1) = 0.57735026918963d0
    gp(2, 2) = -0.57735026918963d0
    gp(2, 3) = -0.57735026918963d0

    gp(3, 1) = -0.57735026918963d0
    gp(3, 2) = 0.57735026918963d0
    gp(3, 3) = -0.57735026918963d0

    gp(4, 1) = 0.57735026918963d0
    gp(4, 2) = 0.57735026918963d0
    gp(4, 3) = -0.57735026918963d0

    gp(5, 1) = -0.57735026918963d0
    gp(5, 2) = -0.57735026918963d0
    gp(5, 3) = 0.57735026918963d0

    gp(6, 1) = 0.57735026918963d0
    gp(6, 2) = -0.57735026918963d0
    gp(6, 3) = 0.57735026918963d0

    gp(7, 1) = -0.57735026918963d0
    gp(7, 2) = 0.57735026918963d0
    gp(7, 3) = 0.57735026918963d0

    gp(8, 1) = 0.57735026918963d0
    gp(8, 2) = 0.57735026918963d0
    gp(8, 3) = 0.57735026918963d0

    gw(1:8) = 1d0

  else if (ngp == 6) then         ! prism

    gp(1, 1) = 0.6666666666666667d0
    gp(1, 2) = 0.1666666666666667d0
    gp(1, 3) = -0.5773502691896257645091488d0

    gp(2, 1) = 0.1666666666666667d0
    gp(2, 2) = 0.6666666666666667d0
    gp(2, 3) = -0.5773502691896257645091488d0

    gp(3, 1) = 0.1666666666666667d0
    gp(3, 2) = 0.1666666666666667d0
    gp(3, 3) = -0.5773502691896257645091488d0

    gp(4, 1) = 0.6666666666666667d0
    gp(4, 2) = 0.1666666666666667d0
    gp(4, 3) = 0.5773502691896257645091488d0

    gp(5, 1) = 0.1666666666666667d0
    gp(5, 2) = 0.6666666666666667d0
    gp(5, 3) = 0.5773502691896257645091488d0

    gp(6, 1) = 0.1666666666666667d0
    gp(6, 2) = 0.1666666666666667d0
    gp(6, 3) = 0.5773502691896257645091488d0

    gw(1) = 0.1666666666666667d0
    gw(2) = 0.1666666666666667d0
    gw(3) = 0.1666666666666667d0
    gw(4) = 0.1666666666666667d0
    gw(5) = 0.1666666666666667d0
    gw(6) = 0.1666666666666667d0

  else
    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPandGWtetpri", ngp
    stop

  end if

end subroutine genGPandGWtetpri

!======================================================================
! Subroutine generates Gaussian Quadrature points and weights
!======================================================================
subroutine genGPandGWtri(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 2), gw(ngp)

  if (ngp == 1) then

    gp(1, 1) = 0.3333333333333333d0
    gp(1, 2) = 0.3333333333333333d0

    gw(1) = 0.5d0

  else if (ngp == 3) then

    gp(1, 1) = 0.6666666666666667d0
    gp(1, 2) = 0.1666666666666667d0

    gp(2, 1) = 0.1666666666666667d0
    gp(2, 2) = 0.6666666666666667d0

    gp(3, 1) = 0.1666666666666667d0
    gp(3, 2) = 0.1666666666666667d0

    gw(1) = 0.166666666666666666667d0
    gw(2) = 0.166666666666666666667d0
    gw(3) = 0.166666666666666666667d0

  else

    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPandGWtri", ngp
    stop

  end if

end subroutine genGPandGWtri

!======================================================================
!
!======================================================================
subroutine genGPtritet(ngp, for, gp)
  implicit none

  integer, intent(in)  :: ngp, for
  real(8), intent(out) :: gp(ngp, 3)

  if (ngp == 1) then

    if (for == 1) then
      gp(1, 1) = 0.3333333333333333d0
      gp(1, 2) = 0.3333333333333333d0
      gp(1, 3) = 0.3333333333333333d0

    else if (for == 2) then
      gp(1, 1) = 0.0d0
      gp(1, 2) = 0.3333333333333333d0
      gp(1, 3) = 0.3333333333333333d0

    else if (for == 3) then
      gp(1, 1) = 0.3333333333333333d0
      gp(1, 2) = 0.0d0
      gp(1, 3) = 0.3333333333333333d0

    else if (for == 4) then
      gp(1, 1) = 0.3333333333333333d0
      gp(1, 2) = 0.3333333333333333d0
      gp(1, 3) = 0.0d0

    end if

  else if (ngp == 3) then

    if (for == 1) then
      gp(1, 1) = 0.166666666666666666667d0
      gp(1, 2) = 0.666666666666666666667d0
      gp(1, 3) = 0.166666666666666666667d0

      gp(2, 1) = 0.166666666666666666667d0
      gp(2, 2) = 0.166666666666666666667d0
      gp(2, 3) = 0.666666666666666666667d0

      gp(3, 1) = 0.666666666666666666667d0
      gp(3, 2) = 0.166666666666666666667d0
      gp(3, 3) = 0.166666666666666666667d0

    else if (for == 2) then
      gp(1, 1) = 0.0d0
      gp(1, 2) = 0.166666666666666666667d0
      gp(1, 3) = 0.166666666666666666667d0

      gp(2, 1) = 0.0d0
      gp(2, 2) = 0.166666666666666666667d0
      gp(2, 3) = 0.666666666666666666667d0

      gp(3, 1) = 0.0d0
      gp(3, 2) = 0.666666666666666666667d0
      gp(3, 3) = 0.166666666666666666667d0

    else if (for == 3) then
      gp(1, 1) = 0.666666666666666666667d0
      gp(1, 2) = 0.0d0
      gp(1, 3) = 0.166666666666666666667d0

      gp(2, 1) = 0.166666666666666666667d0
      gp(2, 2) = 0.0d0
      gp(2, 3) = 0.666666666666666666667d0

      gp(3, 1) = 0.166666666666666666667d0
      gp(3, 2) = 0.0d0
      gp(3, 3) = 0.166666666666666666667d0

    else if (for == 4) then
      gp(1, 1) = 0.166666666666666666667d0
      gp(1, 2) = 0.666666666666666666667d0
      gp(1, 3) = 0.0d0

      gp(2, 1) = 0.666666666666666666667d0
      gp(2, 2) = 0.166666666666666666667d0
      gp(2, 3) = 0.0d0

      gp(3, 1) = 0.166666666666666666667d0
      gp(3, 2) = 0.166666666666666666667d0
      gp(3, 3) = 0.0d0

    end if

  else

    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPtritet"
    stop

  end if

end subroutine genGPtritet

!======================================================================
!
!======================================================================
subroutine genGPtripri(ngp, for, gp)
  implicit none

  integer, intent(in)  :: ngp, for
  real(8), intent(out) :: gp(ngp, 3)

  if (ngp == 1) then

    if (for == 1) then
      gp(1, 1) = 0.3333333333333333d0
      gp(1, 2) = 0.3333333333333333d0
      gp(1, 3) = -1.0d0
    else if (for == 2) then
      gp(1, 1) = 0.3333333333333333d0
      gp(1, 2) = 0.3333333333333333d0
      gp(1, 3) = 1.0d0
    end if

  else if (ngp == 3) then

    if (for == 1) then
      gp(1, 1) = 0.1666666666666667d0
      gp(1, 2) = 0.1666666666666667d0
      gp(1, 3) = -1.0d0

      gp(2, 1) = 0.1666666666666667d0
      gp(2, 2) = 0.6666666666666667d0
      gp(2, 3) = -1.0d0

      gp(3, 1) = 0.6666666666666667d0
      gp(3, 2) = 0.1666666666666667d0
      gp(3, 3) = -1.0d0

    else if (for == 2) then
      gp(1, 1) = 0.6666666666666667d0
      gp(1, 2) = 0.1666666666666667d0
      gp(1, 3) = 1.0d0

      gp(2, 1) = 0.1666666666666667d0
      gp(2, 2) = 0.6666666666666667d0
      gp(2, 3) = 1.0d0

      gp(3, 1) = 0.1666666666666667d0
      gp(3, 2) = 0.1666666666666667d0
      gp(3, 3) = 1.0d0
    end if

  else

    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPtriprm"
    stop

  end if

end subroutine genGPtripri

!======================================================================
! Subroutine generates Gaussian Quadrature points and weights : Prism
!======================================================================
subroutine genGPandGWpri(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 3), gw(ngp)

  if (ngp == 6) then

    gp(1, 1) = 0.6666666666666667d0
    gp(1, 2) = 0.1666666666666667d0
    gp(1, 3) = -0.5773502691896257645091488d0

    gp(2, 1) = 0.1666666666666667d0
    gp(2, 2) = 0.6666666666666667d0
    gp(2, 3) = -0.5773502691896257645091488d0

    gp(3, 1) = 0.1666666666666667d0
    gp(3, 2) = 0.1666666666666667d0
    gp(3, 3) = -0.5773502691896257645091488d0

    gp(4, 1) = 0.6666666666666667d0
    gp(4, 2) = 0.1666666666666667d0
    gp(4, 3) = 0.5773502691896257645091488d0

    gp(5, 1) = 0.1666666666666667d0
    gp(5, 2) = 0.6666666666666667d0
    gp(5, 3) = 0.5773502691896257645091488d0

    gp(6, 1) = 0.1666666666666667d0
    gp(6, 2) = 0.1666666666666667d0
    gp(6, 3) = 0.5773502691896257645091488d0

    gw(1) = 0.1666666666666667d0
    gw(2) = 0.1666666666666667d0
    gw(3) = 0.1666666666666667d0
    gw(4) = 0.1666666666666667d0
    gw(5) = 0.1666666666666667d0
    gw(6) = 0.1666666666666667d0

  else

    write (*, *) "ERROR: Undefined NGAUSS", ngp
    write (*, *) "genGPandGWpri"
    stop

  end if

end subroutine genGPandGWpri
