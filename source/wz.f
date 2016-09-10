! =====================================================================
! subroutine WZ 
!
! This subroutine constructs the W and Z vectors associated with
! applying the Sherman-Morrison algorithm (as outline within the
! numerical methods chapter of this thesis). It has been modifed for
! imaginary time propagation.
! ---------------------------------------------------------------------
      subroutine WZ(nP,dT,TPO,W,Z)
      implicit none

      integer :: nP
      double precision :: dT 
      double precision :: TPO
      double complex :: W(nP) 
      double complex :: Z(nP)

      integer :: i

      W(1) = dcmplx(1.0d0,0.0d0)
      Z(1) = dcmplx(0.5d0*dT*TPO,0.0d0) 
      do i = 2, nP-1
         W(i) = dcmplx(0.0d0,0.0d0)
         Z(i) = dcmplx(0.0d0,0.0d0)
      enddo
      W(nP) = dcmplx(1.0d0,0.0d0)
      Z(nP) = dcmplx(0.5d0*dT*TPO,0.0d0)

      return
      end
! =====================================================================
