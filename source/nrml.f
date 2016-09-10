! =====================================================================
! subroutine NRML 
!
! This subroutine normalizes a wavefunction.
! --------------------------------------------------------------------- 
      subroutine NRML(nRP,nR,nP,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      double precision :: A

      double precision :: NRM

      A = NRM(nRP,nR,nP,RdRP,U)
      A = 1.0d0/dsqrt(A)

      do i = 1, nRP
         U(i) = dcmplx(A,0.0d0)*U(i)
      enddo

      return
      end
! =====================================================================
