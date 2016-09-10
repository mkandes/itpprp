! ===================================================================== 
! subroutine ZSWP
!
! This subroutine swaps the indexed storage a of double precision, 
! complex-valued function defined on a two-dimensional grid that has 
! been stored in a one-dimensional array. 
! ---------------------------------------------------------------------
      subroutine ZSWP(nQQ,nQ1,nQ2,Z,ZT)
      implicit none

      integer :: nQQ 
      integer :: nQ1 
      integer :: nQ2 
      double complex :: Z(nQQ)
      double complex :: ZT(nQQ)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      do j = 1, nQ2
         do i = 1, nQ1
            k = k + 1
            ZT((i-1)*nQ2+j) = Z(k)
         enddo
      enddo

      do k = 1, nQQ
         Z(k) = ZT(k)
      enddo

      return
      end
! =====================================================================
