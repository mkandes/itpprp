! =====================================================================
! subroutine SM 
!
! This subroutine corrects the wavefunction at each time step using the
! Sherman-Morrison algorithm to account for the periodic boundary
! conditions of the polar coordinate system. Note that this subroutine
! does not directly solve the set of linear equations required by the
! Sherman-Morrison algorithm, it merely uses the already computed
! solutions to complete the wavefunction correction.
! --------------------------------------------------------------------- 
      subroutine SM(nP,W,Z,US)
      implicit none

      integer :: nP
      double complex :: W(nP)
      double complex :: Z(nP)
      double complex :: US(nP)

      integer :: i
      double complex :: d
      double complex :: dDEN
      double complex :: dNUM

      dDEN = dcmplx(0.0d0,0.0d0)
      dNUM = dcmplx(0.0d0,0.0d0)

      do i = 1, nP
         dDEN = dDEN + Z(i)*W(i)
         dNUM = dNUM + Z(i)*US(i)
      enddo
      dDEN = dcmplx(1.0d0,0.0d0) - dDEN

      d = dNUM/dDEN

      do i = 1, nP
         US(i) = US(i) + d*W(i)
      enddo

      return
      end
! =====================================================================
