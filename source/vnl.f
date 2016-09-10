! =====================================================================
! subroutine VNL 
!
! This subroutine constructs the nonlinear, mean-field interaction
! potential in the Gross-Pitaevskii equation.
! ---------------------------------------------------------------------
      subroutine VNL(nQQ,g,V,U) 
      implicit none

      integer :: nQQ
      double precision :: g 
      double precision :: V(nQQ) 
      double complex :: U(nQQ)

      integer :: i

      do i = 1, nQQ
         V(i) = g*cdabs(U(i))**2
      enddo

      return
      end
! =====================================================================
