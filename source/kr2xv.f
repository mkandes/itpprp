! =====================================================================
! function KR2XV 
!
! This function computes the expectation value of the radial momentum
! operator squared.
! ---------------------------------------------------------------------
      double precision function KR2XV(nRP,nR,nP,dR,R,RdRP,U)
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
      double precision :: XV
      double precision :: XV1
      double precision :: XV2
      double precision :: XV3
      double precision :: XV4
      double complex :: a
      double complex :: b
      double complex :: c
      double complex :: d
      double complex :: e

      KR2XV = 0.0d0
      a = dcmplx(2.0d0/dR**2,0.0d0)
      b = dcmplx(1.0d0/dR**2,0.0d0)

      k = 0
      do j = 1, nP
         k = k + 1
         c = dcmplx(-0.25d0/(dR*R(1)**2),0.0d0)
         d = dcmplx(-0.25d0/R(1)**4,0.0d0)
         e = dcmplx(1.0d0/R(1)**3,0.0d0)
         XV1 = dble(dconjg(U(k))*(a*U(k)-b*U(k+1)))
         XV2 = dble(dconjg(U(k))*c*U(k+1))
         XV3 = dble(dconjg(U(k))*d*U(k))
         XV4 = dble(dconjg(U(k))*e*U(k))
         XV = XV1 + XV2 + XV3 + XV4
         KR2XV = KR2XV + XV*RdRP(1)
         do i = 2, nR-1
            k = k + 1
            c = dcmplx(-0.25d0/(dR*R(i)**2),0.0d0)
            d = dcmplx(-0.25d0/R(i)**4,0.0d0)
            e = dcmplx(1.0d0/R(i)**3,0.0d0)
            XV1 = dble(dconjg(U(k))*(a*U(k)-b*(U(k-1)+U(k+1)))) 
            XV2 = dble(dconjg(U(k))*c*(U(k+1)-U(k-1)))
            XV3 = dble(dconjg(U(k))*d*U(k))
            XV4 = dble(dconjg(U(k))*e*U(k))
            XV = XV1 + XV2 + XV3 + XV4
            KR2XV = KR2XV + XV*RdRP(i)
         enddo
         k = k + 1
         c = dcmplx(-0.25d0/(dR*R(nR)**2),0.0d0)
         d = dcmplx(-0.25d0/R(nR)**4,0.0d0)
         e = dcmplx(1.0d0/R(nR)**3,0.0d0)
         XV1 = dble(dconjg(U(k))*(a*U(k)-b*U(k-1)))
         XV2 = dble(dconjg(U(k))*c*-U(k-1))
         XV3 = dble(dconjg(U(k))*d*U(k))
         XV4 = dble(dconjg(U(k))*e*U(k))
         XV = XV1 + XV2 + XV3 + XV4
         KR2XV = KR2XV + XV*RdRP(nR)
      enddo

      return
      end
! =====================================================================
