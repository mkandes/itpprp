! =====================================================================
! subroutine RHSP 
!
! This subroutine computes the matrix-vector multiplication required to
! construct the right-hand side vector of the Peaceman-Rachford step
! involving the angular kinetic energy terms. It has been modified for
! imaginary time propagation.
! ---------------------------------------------------------------------
      subroutine RHSP(nRP,nR,nP,dT,TPD,TPO,V,U,UT)
      implicit none

      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: dT
      double precision :: TPD(nR)
      double precision :: TPO(nR)
      double precision :: V(nRP)
      double complex :: U(nRP)
      double complex :: UT(nRP)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      do j = 1, nR
         k = k + 1
         UT(k) = dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k+nP-1)
         UT(k) = UT(k)+dcmplx(1.0d0-0.5d0*dT*TPD(j)- 0.25d0*dT*V(k),0.0d0)*U(k)
         UT(k) = UT(k)+dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k+1)
         do i = 2, nP-1
            k = k + 1
            UT(k) = dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k-1)
            UT(k) = UT(k)+dcmplx(1.0d0-0.5d0*dT*TPD(j)- 0.25d0*dT*V(k),0.0d0)*U(k)
            UT(k) = UT(k)+dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k+1)
         enddo
         k = k + 1
         UT(k) = dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k-1)
         UT(k) = UT(k)+dcmplx(1.0d0-0.5d0*dT*TPD(j)- 0.25d0*dT*V(k),0.0d0)*U(k)
         UT(k) = UT(k)+dcmplx(0.5d0*dT*TPO(j),0.0d0)*U(k-nP+1)
      enddo

      do k = 1, nRP
         U(k) = UT(k)
      enddo

      return
      end
! =====================================================================
