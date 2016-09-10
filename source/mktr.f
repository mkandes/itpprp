! =====================================================================
! subroutine MKTR 
!
! This subroutine constructs the radial kinetic energy terms that
! appear in the finite difference operators.
! ---------------------------------------------------------------------
      subroutine MKTR(nR,dR,R,RH,TRL,TRD,TRU)
      implicit none

      integer :: nR
      double precision :: dR
      double precision :: R(nR)
      double precision :: RH(nR+1)
      double precision :: TRL(nR-1)
      double precision :: TRD(nR)
      double precision :: TRU(nR-1)

      integer :: i

      do i = 1, nR-1
         TRL(i) = 0.5d0*RH(i+1)/(R(i+1)*dR**2)
         TRD(i) = 0.5d0*(RH(i)+RH(i+1))/(R(i)*dR**2)
         TRU(i) = 0.5d0*RH(i+1)/(R(i)*dR**2)
      enddo
      TRD(nR) = 0.5d0*(RH(nR)+RH(nR+1))/(R(nR)*dR**2)

      return
      end
! =====================================================================

