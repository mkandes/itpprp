! =====================================================================
! subroutine LHSP 
!
! This subroutine constructs the left-hand side matrix for the
! Peaceman-Rachford step involving the angular kinetic energy
! terms. It has been modified for imaginary time propagation.
! ---------------------------------------------------------------------
      subroutine LHSP(nP,dT,TPD,TPO,VS,DL,D,DU)
      implicit none

      integer :: nP
      double precision :: dT 
      double precision :: TPD 
      double precision :: TPO 
      double precision :: VS(nP)
      double complex :: DL(nP-1)
      double complex :: D(nP) 
      double complex :: DU(nP-1)

      integer :: i

      DL(1) = dcmplx(-0.5d0*dT*TPO,0.0d0)
      D(1) = dcmplx(1.0d0+0.5d0*dT*TPD+0.25d0*dT*VS(1)+0.5d0*dT*TPO,0.0d0)
      DU(1) = dcmplx(-0.5d0*dT*TPO,0.0d0)
      do i = 2, nP-1
         DL(i) = dcmplx(-0.5d0*dT*TPO,0.0d0)
         D(i) = dcmplx(1.0d0+0.5d0*dT*TPD+0.25d0*dT*VS(i),0.0d0)
         DU(i) = dcmplx(-0.5d0*dT*TPO,0.0d0)
      enddo
      D(nP) = dcmplx(1.0d0+0.5d0*dT*TPD+0.25d0*dT*VS(nP)+0.5d0*dT*TPO,0.0d0)

      return
      end
! =====================================================================
