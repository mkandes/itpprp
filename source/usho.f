! =====================================================================
! subroutine USHO 
!
! This subroutine constructs a simple harmonic oscillator ground state
! wavefunction.
! --------------------------------------------------------------------- 
      subroutine USHO(nRP,nR,nP,x0,y0,kX,kY,R,P,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: x0
      double precision :: y0
      double precision :: kX
      double precision :: kY
      double precision :: R(nR)
      double precision :: P(nP)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double precision :: pi 
      double complex :: a 
      double complex :: b 
      double complex :: bX 
      double complex :: bY

      pi = 3.14159265358979323846d0
      a = dcmplx((kX*kY/pi**2)**(0.25d0),0.0d0)

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            bX = dcmplx(-0.5d0*kX*(R(i)*dcos(P(j))-x0)**2,0.0d0) 
            bY = dcmplx(-0.5d0*kY*(R(i)*dsin(P(j))-y0)**2,0.0d0) 
            b = bX + bY
            U(k) = a*cdexp(b)
         enddo
      enddo

      return
      end
! =====================================================================
