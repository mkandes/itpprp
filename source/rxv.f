! =====================================================================
! function RXV 
!
! This function computes the expectation value of the radial position
! operator.
! ---------------------------------------------------------------------
      double precision function RXV(nRP,nR,nP,R,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: R(nR)
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      RXV = 0.0d0

      k = 0
      do j = 1, nP
         do i = 1, nR 
            k = k + 1
            RXV = RXV + R(i)*cdabs(U(k))**2*RdRP(i)
         enddo
      enddo

      return
      end
! =====================================================================
