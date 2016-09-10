! =====================================================================
! program ITPPRP 
!
! This program finds the ground state wavefunction of the
! time-independent Gross-Pitaevskii equation in a polar coordinate
! system using imaginary time propagation. ITPPRP is a modified version
! of NLSPRP written to propagate in imaginary time. All modified
! functions and/or subroutines required for this task are given below.
! All other functions and/or subroutines used by ITPPRP are listed with
! NLSPRP.
! ---------------------------------------------------------------------
      program ITPPRP
      implicit none

      double precision, parameter :: pi = 3.14159265358979323846d0

      integer :: nT
      integer :: nTskp
      integer :: nR
      integer :: nP
      integer :: qnR
      integer :: qmL
      double precision :: Ta
      double precision :: Tb 
      double precision :: Ra
      double precision :: Rb
      double precision :: Ro
      double precision :: kR
      double precision :: Po
      double precision :: kP
      double precision :: Xo
      double precision :: kX
      double precision :: Yo
      double precision :: kY
      double precision :: g2

      integer :: i, j, k, m
      integer :: nRP
      integer :: nS
      integer :: nFile
      integer :: nInfo
      double precision :: T
      double precision :: dT
      double precision :: dR
      double precision :: dP 
      double precision :: dRP

      double precision, allocatable :: R(:)
      double precision, allocatable :: RH(:)
      double precision, allocatable :: P(:)
      double precision, allocatable :: RdRP(:)
      double precision, allocatable :: TRL(:)
      double precision, allocatable :: TRD(:)
      double precision, allocatable :: TRU(:)
      double precision, allocatable :: TPD(:)
      double precision, allocatable :: TPO(:)
      double precision, allocatable :: V(:)
      double precision, allocatable :: VT(:)
      double precision, allocatable :: VS(:)
      double precision, allocatable :: X(:)
      double complex, allocatable :: U(:) 
      double complex, allocatable :: UT(:)
      double complex, allocatable :: US(:)
      double complex, allocatable :: UR(:)
      double complex, allocatable :: DL(:)
      double complex, allocatable :: D(:) 
      double complex, allocatable :: DU(:) 
      double complex, allocatable :: W(:)
      double complex, allocatable :: Z(:)

      double precision :: NRM
      double precision :: RXV
      double precision :: R2XV
      double precision :: PXV
      double precision :: P2XV
      double precision :: KRXV
      double precision :: KR2XV
      double precision :: KPXV
      double precision :: KP2XV
      double precision :: TRXV 
      double precision :: TPXV
      double precision :: VXV

      read(5,*) nT, nTskp
      read(5,*) nR, nP
      read(5,*) Ta, Tb 
      read(5,*) Ra, Rb
      read(5,*) Ro, kR 
      read(5,*) Po, kP
      read(5,*) Xo, kX
      read(5,*) Yo, kY
      read(5,*) g2
      read(5,*) qnR, qmL

      nRP = nR*nP
      if (nR.ge.nP) then
         nS = nR
      else
         nS = nP
      endif
      nFile = 999
      nInfo = 0
      Po = 2.0d0*pi*Po
      T = 0.0d0
      dT = (Tb-Ta)/dfloat(nT)
      dR = (Rb-Ra)/dfloat(nR)
      dP = 2.0d0*pi/dfloat(nP)
      dRP = dR*dP

      allocate(R(nR),RH(nR+1),P(nP),RdRP(nR),TRL(nR-1),TRD(nR),TRU(nR-1))
      allocate(TPD(nR),TPO(nR),V(nRP),VT(nRP),VS(nS),X(15))
      allocate(U(nRP),UT(nRP),US(nS),UR(nP),DL(nS-1),D(nS),DU(nS))
      allocate(W(nP),Z(nP))

      call GRDR(nR,rA,rB,dR,R,RH)
      call GRDP(nP,dP,P)
      call MKRdRP(nR,dRP,R,RdRP)
      call MKTR(nR,dR,R,RH,TRL,TRD,TRU) 
      call MKTP(nR,dP,R,TPD,TPO)

      call USHO(nRP,nR,nP,Xo,Yo,kX,kY,R,P,U)
      !call USHOR(nRP,nR,nP,Ro,kR,R,P,U) 
      !call USHO2(nRP,nR,nP,qnR,qmL,kR,R,P,U)
      call NRML(nRP,nR,nP,RdRP,U)

      do i = 1, nT
         T = dfloat(i-1)*dT
         if (mod(i-1,nTskp).eq.0) then
            X(1) = NRM(nRP,nR,nP,RdRP,U)
            X(2) = RXV(nRP,nR,nP,R,RdRP,U)
            X(3) = R2XV(nRP,nR,nP,R,RdRP,U)
            X(4) = PXV(nRP,nR,nP,P,RdRP,U)
            X(5) = P2XV(nRP,nR,nP,P,RdRP,U)
            X(6) = KRXV(nRP,nR,nP,dR,R,RdRP,U)
            X(7) = KR2XV(nRP,nR,nP,dR,R,RdRP,U)
            X(10) = TRXV(nRP,nR,nP,RdRP,TRL,TRD,TRU,U)
            call ZSWP(nRP,nR,nP,U,UT)
            X(8) = KPXV(nRP,nR,nP,dP,RdRP,U)
            X(9) = KP2XV(nRP,nR,nP,dP,RdRP,U)
            X(11) = TPXV(nRP,nR,nP,RdRP,TPD,TPO,U)
            call ZSWP(nRP,nP,nR,U,UT)
            !call VAHO(nRP,nR,nP,Ro,Po,kR,kP,R,P,VT)
            call VSHO(nRP,nR,nP,Xo,Yo,kX,kY,R,P,VT)
            !call VSHOR(nRP,nR,nP,Ro,kR,R,VT)
            X(12) = VXV(nRP,nR,nP,RdRP,VT,U)
            call VNL(nRP,g2,VT,U)
            X(14) = VXV(nRP,nR,nP,RdRP,VT,U)
            X(15) = X(10) + X(11) + X(12) + X(14)
            write(900,*) T, X(1), X(2), X(3), X(4), X(5), X(6), X(7),&
               & X(8), X(9), X(10), X(11), X(12), X(13), X(14), X(15)
            call WRTUU(nFile,nRP,nR,nP,R,P,U) 
            nFile = nFile + 1
         endif
         !call VAHO(nRP,nR,nP,Ro,Po,kR,kP,R,P,VT)
         call VSHO(nRP,nR,nP,Xo,Yo,kX,kY,R,P,VT)
         !call VSHOR(nRP,nR,nP,Ro,kR,R,VT)
         V = VT
         call VNL(nRP,g2,VT,U)
         V = V + VT
         call RHSR(nRP,nR,nP,dT,TRL,TRD,TRU,V,U,UT)
         call DSWP(nRP,nR,nP,V,VT)
         call ZSWP(nRP,nR,nP,U,UT)
         do j = 1, nR
            m = (j-1)*nP
            do k = 1, nP
               VS(k) = V(m+k)
               US(k) = U(m+k)
            enddo
            call LHSP(nP,dT,TPD(j),TPO(j),VS,DL,D,DU)
            call ZGTSV(nP,1,DL,D,DU,US,nP,nInfo)
            call LHSP(nP,dT,TPD(j),TPO(j),VS,DL,D,DU)
            call WZ(nP,dT,TPO(j),W,Z)
            call ZGTSV(nP,1,DL,D,DU,W,nP,nInfo) 
            call SM(nP,W,Z,US)
            do k = 1, nP
               U(m+k) = US(k)
            enddo
         enddo
         call RHSP(nRP,nR,nP,dT,TPD,TPO,V,U,UT)
         call DSWP(nRP,nP,nR,V,VT)
         call ZSWP(nRP,nP,nR,U,UT)
         do j = 1, nP
            m = (j-1)*nR
            do k = 1, nR
               VS(k) = V(m+k)
               US(k) = U(m+k)
            enddo
            call LHSR(nR,dT,TRL,TRD,TRU,VS,DL,D,DU)
            call ZGTSV(nR,1,DL,D,DU,US,nR,nInfo)
            do k = 1, nR
               U(m+k) = US(k)
            enddo
         enddo
         call NRML(nRP,nR,nP,RdRP,U)
      enddo

      call WRTU(999,nRP,nR,nP,R,P,U)
      call WRTUU(1999,nRP,nR,nP,R,P,U)

900   format(1X,16(F23.15))

      deallocate(R,RH,P,RdRP,TRL,TRD,TRU,TPD,TPO,V,VT,VS,X)
      deallocate(U,UT,US,UR,DL,D,DU,W,Z)

      stop
      end
! =====================================================================
