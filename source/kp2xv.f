! =====================================================================
! function KP2XV
!
! This function computes the expectation value of the angular momentum
! operator squared.
! --------------------------------------------------------------------- 
      double precision function KP2XV(nRP,nR,nP,dP,RdRP,U)
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
      double precision :: XV
      double complex :: a 
      double complex :: b

      KP2XV = 0.0d0
      a = dcmplx(2.0d0/dP**2,0.0d0)
      b = dcmplx(1.0d0/dP**2,0.0d0)

      k = 0
      do j = 1, nR
      k = k + 1
         XV = dble(dconjg(U(k))*(a*U(k)-b*(U(k+nP-1)+U(k+1))))
         KP2XV = KP2XV + XV*RdRP(j)
         do i = 2, nP-1
            k = k + 1
            XV = dble(dconjg(U(k))*(a*U(k)-b*(U(k-1)+U(k+1))))
            KP2XV = KP2XV + XV*RdRP(j)
         enddo
         k = k + 1
         XV = dble(dconjg(U(k))*(a*U(k)-b*(U(k-1)+U(k-nP+1)))) 
         KP2XV = KP2XV + XV*RdRP(j)
      enddo

      return
      end
! =====================================================================

