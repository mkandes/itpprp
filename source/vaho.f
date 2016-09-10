! =====================================================================
! subroutine VAHO 
!
! This subroutine constructs a polar simple harmonic oscillator
! potential. It may be used to generate an inital wavefuction that is
! correctly aligned with the potential minimum of the simple harmonic
! oscillator ring.
! ---------------------------------------------------------------------
      subroutine VAHO(nRP,nR,nP,r0,p0,kR,kP,R,P,V)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: pi
      double precision :: r0
      double precision :: p0
      double precision :: kR
      double precision :: kP
      double precision :: R(nR)
      double precision :: P(nP)
      double precision :: V(nRP)

      integer :: i
      integer :: j
      integer :: k

      pi = 3.14159265358979323846d0

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            if (P(j).le.(1.5*pi)) then
               V(k) = 0.5d0*(kR*(R(i)-r0))**2+0.5d0*(kP*r0*(P(j)-p0))**2
            else
               V(k) = 0.5d0*(kR*(R(i)-r0))**2+0.5d0*(kP*r0*(P(j)-p0- 2.0d0*pi))**2
            endif
         enddo
      enddo

      return
      end
! =====================================================================
