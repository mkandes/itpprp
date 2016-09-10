! =====================================================================
! subroutine USHO2
!
! This subroutine computes the analytic solution of the time-independent
! Schrodinger equation for a two-dimensional, symmetric simple harmonic
! oscillator potential in polar coordinates.
! --------------------------------------------------------------------- 
      subroutine USHO2(nRP,nR,nP,qnR,qmL,kR,R,P,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      integer :: qnR
      integer :: qmL
      double precision :: kR
      double precision :: R(nR)
      double precision :: P(nP)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double precision :: pi 

      integer :: factorial
      double precision :: alaguerre

      pi = 3.14159265358979323846d0

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            U(k) = CMPLX( &
               & SQRT((kR**(ABS(qmL)+1) * REAL(factorial(qnR))) / & 
               & (pi * REAL(factorial(qnR + ABS(qmL))))) * & 
               & R(i)**ABS(qmL) * &
               & alaguerre(qnR, ABS(qmL), kR * R(i)**2) * & 
               & EXP(-0.5 * kR *R(i)**2) , 0.0 ) * & 
               & EXP(CMPLX(0.0, REAL(qmL)*P(j))) 
         enddo
      enddo

      return
      end
! =====================================================================
