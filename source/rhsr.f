! =====================================================================
! subroutine RHSR 
!
! This subroutine computes the matrix-vector multiplication required to
! construct the right-hand side vector of the Peaceman-Rachford step
! involving the radial kinetic energy terms. It has been modified for
! imaginary time propagation.
! ---------------------------------------------------------------------
      subroutine RHSR(nRP,nR,nP,dT,TRL,TRD,TRU,V,U,UT)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: dT
      double precision :: TRL(nR-1)
      double precision :: TRD(nR)
      double precision :: TRU(nR-1)
      double precision :: V(nRP)
      double complex :: U(nRP) 
      double complex :: UT(nRP)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      do j = 1, nP
         k = k + 1
         UT(k) = dcmplx(1.0d0-0.5d0*dT*TRD(1)-0.25d0*dT*V(k),0.0d0)*U(k)
         UT(k) = UT(k)+dcmplx(0.5d0*dT*TRU(1),0.0d0)*U(k+1)
         do i = 2, nR-1
            k = k + 1
            UT(k) = dcmplx(0.5d0*dT*TRL(i-1),0.0d0)*U(k-1)
            UT(k) = UT(k)+dcmplx(1.0d0-0.5d0*dT*TRD(i)- 0.25d0*dT*V(k),0.0d0)*U(k)
            UT(k) = UT(k)+dcmplx(0.5d0*dT*TRU(i),0.0d0)*U(k+1)
         enddo
         k = k + 1
         UT(k) = dcmplx(0.5d0*dT*TRL(nR-1),0.0d0)*U(k-1)
         UT(k) = UT(k)+dcmplx(1.0d0-0.5d0*dT*TRD(nR)- 0.25d0*dT*V(k),0.0d0)*U(k)
      enddo

      do k = 1, nRP
         U(k) = UT(k)
      enddo

      return
      end
! =====================================================================
