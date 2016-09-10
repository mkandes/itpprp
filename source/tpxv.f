! =====================================================================
! function TPXV 
!
! This function computes the expectation value of the angular kinetic
! energy operator.
! ---------------------------------------------------------------------
      double precision function TPXV(nRP,nR,nP,RdRP,TPD,TPO,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: RdRP(nR)
      double precision :: TPD(nR)
      double precision :: TPO(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      double complex :: a
      double complex :: b

      TPXV = 0.0d0

      k = 0
      do j = 1, nR
         k = k + 1
         a = dcmplx(TPD(j),0.0d0)
         b = dcmplx(TPO(j),0.0d0)
         TPXV = TPXV + dble(dconjg(U(k))*(a*U(k)-b*(U(k+nP-1)+ U(k+1))))*RdRP(j)
         do i = 2, nP-1
            k = k + 1
            TPXV = TPXV + dble(dconjg(U(k))*(a*U(k)-b*(U(k-1)+ U(k+1))))*RdRP(j)
         enddo
         k = k + 1
         TPXV = TPXV + dble(dconjg(U(k))*(a*U(k)-b*(U(k-1)+ U(k-nP+1))))*RdRP(j)
      enddo

      return
      end
! =====================================================================
