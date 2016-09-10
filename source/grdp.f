! =====================================================================
! subroutine GRDP
!
! This subroutine generates the angular coordinate grid points.
! ---------------------------------------------------------------------
      subroutine GRDP(nP,dP,P)
      implicit none

      integer :: nP
      double precision :: dP
      double precision :: P(nP)

      integer :: i

      do i = 1, nP
         P(i) = dfloat(i-1)*dP
      enddo

      return
      end
! =====================================================================

