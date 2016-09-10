! =====================================================================
! function VXV 
!
! This function computes the expectation value of any potential operator,
! except for the nonlinear, mean-field interaction potential of the GPE.
! ---------------------------------------------------------------------
      double precision function VXV(nRP,nR,nP,RdRP,V,U) 
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: RdRP(nR)
      double precision :: V(nRP)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double precision :: XV

      VXV = 0.0d0

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            XV = V(k)*cdabs(U(k))**2
            VXV = VXV + XV*RdRP(i)
         enddo
      enddo

      return
      end
! =====================================================================
