! =====================================================================
!
! function KRXV 
!
! This function computes the expectation value of the radial momentum
! operator.
!
! ---------------------------------------------------------------------
      double precision function KRXV(nRP,nR,nP,dR,R,RdRP,U)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: dR
      double precision :: R(nR)
      double precision :: RdRP(nR)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k
      double complex :: a
      double complex :: b

      KRXV = 0.0d0
      a = dcmplx(0.0d0,-0.5d0/dR)

      k = 0
      do j = 1, nP
         k = k + 1
         b = dcmplx(0.0d0,-0.5d0/R(1)**2)
         KRXV = KRXV + dble(dconjg(U(k))*(a*U(k+1)+b*U(k)))*RdRP(1)
         do i = 2, nR-1
            k = k + 1
            b = dcmplx(0.0d0,-0.5d0/R(i)**2)
            KRXV = KRXV + dble(dconjg(U(k))*(a*(U(k+1)-U(k-1))+ b*U(k)))*RdRP(i)
         enddo
         k = k + 1
         b = dcmplx(0.0d0,-0.5d0/R(nR)**2)
         KRXV = KRXV + dble(dconjg(U(k))*(a*-U(k-1)+b*U(k)))*RdRP(nR)
      enddo

      return
      end
! =====================================================================

