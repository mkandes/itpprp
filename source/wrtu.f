! =====================================================================
! subroutine WRTU 
!
! This subroutine writes to file the wavefuction.
! --------------------------------------------------------------------- 
      subroutine WRTU(nFile,nRP,nR,nP,R,P,U)
      implicit none

      integer :: nFile
      integer :: nRP
      integer :: nR
      integer :: nP
      double precision :: R(nR)
      double precision :: P(nP)
      double complex :: U(nRP)

      integer :: i
      integer :: j
      integer :: k

      k = 0
      open(unit=nFile)
      do j = 1, nP
         do i = 1, nR 
            k = k + 1
            write(nFile,*) U(k)
         enddo
      enddo
      close(unit=nFile)

      return
      end
! =====================================================================
