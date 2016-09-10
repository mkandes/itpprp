! =====================================================================
! subroutine VSHO 
!
! This subroutine constructs a simple harmonic oscillator potential.
! --------------------------------------------------------------------- 
      subroutine VSHO(nRP,nR,nP,x0,y0,kX,kY,R,P,V)
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
      double precision :: V(nRP)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            V(k) = 0.5d0*(kX*(R(i)*dcos(P(j))-x0))**2
            V(k) = V(k) + 0.5d0*(kY*(R(i)*dsin(P(j))-y0))**2
         enddo
      enddo

      return
      end
! =====================================================================
