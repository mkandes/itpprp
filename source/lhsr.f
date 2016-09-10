! =====================================================================
! subroutine LHSR 
!
! This subroutine constructs the left-hand side matrix of the
! Peaceman-Rachford step involving the radial kinetic energy
! terms. It has been modified for imaginary time propagation.
! ---------------------------------------------------------------------
      subroutine LHSR(nR,dT,TRL,TRD,TRU,VS,DL,D,DU)
      implicit none

      integer :: nR
      double precision :: dT
      double precision :: TRL(nR-1)
      double precision :: TRD(nR)
      double precision :: TRU(nR-1)
      double precision :: VS(nR)
      double complex :: DL(nR-1)
      double complex :: D(nR)
      double complex :: DU(nR-1)

      integer :: i

      do i = 1, nR-1
         DL(i) = dcmplx(-0.5d0*dT*TRL(i),0.0d0)
         D(i) = dcmplx(1.0d0+0.5d0*dT*TRD(i)+0.25d0*dT*VS(i),0.0d0)
         DU(i) = dcmplx(-0.5d0*dT*TRU(i),0.0d0)
      enddo
      D(nR) = dcmplx(1.0d0+0.5d0*dT*TRD(nR)+0.25d0*dT*VS(nR),0.0d0)

      return
      end
! =====================================================================
