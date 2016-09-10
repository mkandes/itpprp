! =====================================================================
! subroutine DSWP
! 
! This subroutine swaps the indexed storage of a double precision, 
! real-valued function defined on a two-dimensional grid that has been 
! stored in a one-dimensional array.
! ---------------------------------------------------------------------
      subroutine DSWP(nQQ,nQ1,nQ2,D,DT)
      implicit none

      integer :: nQQ
      integer :: nQ1
      integer :: nQ2
      double precision :: D(nQQ)
      double precision :: DT(nQQ)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      do j = 1, nQ2
         do i = 1, nQ1
            k = k + 1
            DT((i-1)*nQ2+j) = D(k)
         enddo
      enddo

      do k = 1, nQQ
         D(k) = DT(k)
      enddo

      return
      end
! =====================================================================
