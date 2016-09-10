! =====================================================================
! function NRM 
!
! This function computes the normalization constant of a wavefunction.
! ---------------------------------------------------------------------
      double precision function NRM(nRP,nR,nP,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      NRM = 0.0d0 

      k = 0
      do j = 1, nP
         do i = 1, nR
            k = k + 1
            NRM = NRM + cdabs(U(k))**2*RdRP(i)
         enddo
      enddo

      return
      end
! =====================================================================
