! =====================================================================
! subroutine WRTUU
!
! This subroutine writes to file the density of the wavefunction defined
! on the polar grid. Note that the output is formatted to accommodate the
! splot function of gnuplot (using cylindrical mapping).
! --------------------------------------------------------------------- 
      subroutine WRTUU(nFile,nRP,nR,nP,R,P,U)
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
            write(nFile,*) P(j), cdabs(U(k))**2, R(i)
         enddo
         write(nFile,*) ' '
      enddo
      close(unit=nFile)

      return
      end
! =====================================================================
