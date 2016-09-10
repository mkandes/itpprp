! =====================================================================
! subroutine USHOR 
!
! This subroutine constructs an approximate ground state
! wavefunction for the simple harmonic oscillator ring.
! ---------------------------------------------------------------------
      subroutine USHOR(nRP,nR,nP,r0,kR,R,P,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: r0
      double precision :: kR
      double precision :: R(nR)
      double precision :: P(nP)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double precision :: pi
      double complex :: a
      double complex :: b

      pi = 3.14159265358979323846d0
      a = dcmplx((kR/pi)**(0.25d0),0.0d0)

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            b = dcmplx(-0.5d0*kR*(R(i)-r0)**2,0.0d0)
            U(k) = a*cdexp(b)
         enddo
      enddo

      return
      end
! =====================================================================
