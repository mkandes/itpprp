! =====================================================================
! subroutine MKTP 
!
! This subroutine constructs the angular kinetic energy terms that
! appear in the finite difference operators.
! ---------------------------------------------------------------------
      subroutine MKTP(nR,dP,R,TPD,TPO)
      implicit none

      integer :: nR
      double precision :: dP 
      double precision :: R(nR)
      double precision :: TPD(nR)
      double precision :: TPO(nR)

      integer :: i

      do i = 1, nR
         TPD(i) = 1.0d0/(R(i)*dP)**2
         TPO(i) = 0.5d0/(R(i)*dP)**2
      enddo

      return
      end
! =====================================================================
