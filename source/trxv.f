! =====================================================================
! function TRXV 
!
! This function computes the expectation value of the radial kinetic
! energy operator.
! ---------------------------------------------------------------------
      double precision function TRXV(nRP,nR,nP,RdRP,TRL,TRD,TRU,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: RdRP(nR)
      double precision :: TRL(nR-1)
      double precision :: TRD(nR)
      double precision :: TRU(nR-1)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      double complex :: a
      double complex :: b
      double complex :: c

      TRXV = 0.0d0

      k = 0
      do j = 1, nP
         k = k + 1
         b = dcmplx(TRD(1),0.0d0)
         c = dcmplx(TRU(1),0.0d0)
         TRXV = TRXV + dble(dconjg(U(k))*(b*U(k)-c*U(k+1)))*RdRP(1)
         do i = 2, nR-1
            k = k + 1
            a = dcmplx(TRL(i-1),0.0d0)
            b = dcmplx(TRD(i),0.0d0)
            c = dcmplx(TRU(i),0.0d0)
            TRXV = TRXV + dble(dconjg(U(k))*(b*U(k)-a*U(k-1)- c*U(k+1)))*RdRP(i)
         enddo
         k = k + 1
         a = dcmplx(TRL(nR-1),0.0d0)
         b = dcmplx(TRD(nR),0.0d0)
         TRXV = TRXV + dble(dconjg(U(k))*(b*U(k)-a*U(k-1)))*RdRP(nR)
      enddo

      return
      end
! =====================================================================
