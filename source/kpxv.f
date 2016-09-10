! =====================================================================
! function KPXV
!
! This function computes the expectation value of the angular momentum
! operator.
! --------------------------------------------------------------------- 
      double precision function KPXV(nRP,nR,nP,dP,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: dP
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double complex :: a

      KPXV = 0.0d0
      a = dcmplx(0.0d0,-0.5d0/dP)

      k = 0
      do j = 1, nR
         k = k + 1
         KPXV = KPXV + dble(dconjg(U(k))*a*(U(k+1)-U(k+nP-1)))*RdRP(j)
         do i = 2, nP-1
            k = k + 1
            KPXV = KPXV + dble(dconjg(U(k))*a*(U(k+1)-U(k-1)))*RdRP(j)
         enddo
         k = k + 1
         KPXV = KPXV + dble(dconjg(U(k))*a*(U(k-nP+1)-U(k-1)))*RdRP(j)
      enddo

      return
      end
! =====================================================================
