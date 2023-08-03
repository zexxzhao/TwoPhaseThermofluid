!======================================================================
! Subroutine generates Gaussian Quadrature points and weights for
! NURBS or T-spline
!======================================================================
subroutine genGPandGW_shell(gp, gw, NGAUSS)

  implicit none

  integer, intent(in) :: NGAUSS
  real(8) :: gp(NGAUSS), gw(NGAUSS)

  if (NGAUSS == 2) then
    gp(1) = -0.5773502691896257645091488d0
    gp(2) = 0.5773502691896257645091488d0
    gw(1) = 1.0d0
    gw(2) = 1.0d0
  else if (NGAUSS == 3) then
    gp(1) = -0.77459666924148337703585308d0
    gp(2) = 0.0d0
    gp(3) = 0.77459666924148337703585308d0
    gw(1) = 0.55555555555555555555555556d0
    gw(2) = 0.88888888888888888888888889d0
    gw(3) = 0.55555555555555555555555556d0
  else if (NGAUSS == 4) then
    gp(1) = -0.86113631159405257524d0
    gp(2) = -0.33998104358485626481d0
    gp(3) = 0.33998104358485626481d0
    gp(4) = 0.86113631159405257524d0
    gw(1) = 0.34785484513745385736d0
    gw(2) = 0.65214515486254614264d0
    gw(3) = 0.65214515486254614264d0
    gw(4) = 0.34785484513745385736d0
  else if (NGAUSS == 5) then
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
  else if (NGAUSS == 6) then
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
  else if (NGAUSS == 7) then
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
  else if (NGAUSS == 8) then
    gp(1) = -0.9602898564975362d0
    gp(2) = -0.7966664774136267d0
    gp(3) = -0.5255324099163290d0
    gp(4) = -0.1834346424956498d0
    gp(5) = 0.1834346424956498d0
    gp(6) = 0.5255324099163290d0
    gp(7) = 0.7966664774136267d0
    gp(8) = 0.9602898564975362d0
    gw(1) = 0.1012285362903763d0
    gw(2) = 0.2223810344533745d0
    gw(3) = 0.3137066458778873d0
    gw(4) = 0.3626837833783620d0
    gw(5) = 0.3626837833783620d0
    gw(6) = 0.3137066458778873d0
    gw(7) = 0.2223810344533745d0
    gw(8) = 0.1012285362903763d0
  else if (NGAUSS == 9) then
    gp(1) = -0.9681602395076261d0
    gp(2) = -0.8360311073266358d0
    gp(3) = -0.6133714327005904d0
    gp(4) = -0.3242534234038089d0
    gp(5) = 0d+0
    gp(6) = 0.3242534234038089d0
    gp(7) = 0.6133714327005904d0
    gp(8) = 0.8360311073266358d0
    gp(9) = 0.9681602395076261d0
    gw(1) = 0.0812743883615744d0
    gw(2) = 0.1806481606948574d0
    gw(3) = 0.2606106964029354d0
    gw(4) = 0.3123470770400029d0
    gw(5) = 0.3302393550012598d0
    gw(6) = 0.3123470770400028d0
    gw(7) = 0.2606106964029355d0
    gw(8) = 0.1806481606948574d0
    gw(9) = 0.0812743883615744d0
  else if (NGAUSS == 10) then
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
  else if (NGAUSS == 15) then
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
  end if

end subroutine genGPandGW_shell

!======================================================================
! Subroutine generates Gaussian Quadrature points and weights for
! linear triangle FEM
!======================================================================
subroutine genGPandGW_tri(gp, gw, ngp)

  implicit none

  integer, intent(in)  :: ngp
  real(8), intent(out) :: gp(ngp, 2), gw(ngp)

  if (ngp == 1) then

    gp(1, 1) = 0.333333333333333333333d0
    gp(1, 2) = 0.333333333333333333333d0

    gw(1) = 0.5d0

  else if (ngp == 3) then

    gp(1, 1) = 0.666666666666666666667d0
    gp(1, 2) = 0.166666666666666666667d0

    gp(2, 1) = 0.166666666666666666667d0
    gp(2, 2) = 0.166666666666666666667d0

    gp(3, 1) = 0.166666666666666666667d0
    gp(3, 2) = 0.166666666666666666667d0

    gw(1) = 0.166666666666666666667d0
    gw(2) = 0.166666666666666666667d0
    gw(3) = 0.166666666666666666667d0

  else

    write (*, *) "ERROR: Undefined NGAUSS", ngp
    stop

  end if
end subroutine genGPandGW_tri