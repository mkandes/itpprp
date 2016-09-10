! =====================================================================
! function P2XV
!
! This functions computes the expectation value of the angular position
! operator squared.
! --------------------------------------------------------------------- 
      double precision function P2XV(nRP,nR,nP,P,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: P(nP)
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      P2XV = 0.0d0

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            P2XV = P2XV + (P(j)*cdabs(U(k)))**2*RdRP(i)
         enddo
      enddo

      return
      end
! =====================================================================
