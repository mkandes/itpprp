! ===================================================================== 
! subroutine MKRdRP
!
! This subroutine generates a one-dimensional array that stores the
! values of the differential area unit defined for the polar grid.
! ---------------------------------------------------------------------
      subroutine MKRdRP(nR,dRP,R,RdRP)
      implicit none

      integer :: nR
      double precision :: dRP
      double precision :: R(nR)
      double precision :: RdRP(nR)

      integer :: i

      do i = 1, nR
         RdRP(i) = R(i)*dRP
      enddo

      return
      end
! =====================================================================
