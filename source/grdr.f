! =====================================================================
! subroutine GRDR
!
! This subroutine generates the radial coordinate grid points.
! --------------------------------------------------------------------- 
      subroutine GRDR(nR,rA,rB,dR,R,RH)
      implicit none

      integer :: nR
      double precision :: rA
      double precision :: rB 
      double precision :: dR
      double precision :: R(nR)
      double precision :: RH(nR+1)

      integer :: i

      do i = 1, nR
         R(i) = rA + dfloat(i)*dR
         RH(i) = R(i) - 0.5d0*dR
      enddo
      RH(nR+1) = rB - 0.5d0*dR

      return
      end
! =====================================================================

