!=======================================================================
      SUBROUTINE SGSCALC
!=======================================================================
!
!     Calculate sub-grid scale eddy viscosity, TNU
!       Park, Lee, Lee & Choi, 2006, A dynamic subgrid-scale eddy viscosity model with a global
!         model coefficient, Phys. Fluids, 18, 125109.
!       Lee, Choi & Park, 2010, Dynamic global model for large eddy simulation of transient
!         flow, Phys. Fluids, 22, 075106.
!
!     Option #1
!       ISGSMOD = 1, Smagorinsky model (Smagorinsky, 1963)
!       ISGSMOD = 3, Vreman model (Vreman, 2004)
!
!     Option #2
!       IDYNMON = 0, Constant model
!       IDYNMON = 1, Dynamic procedure (Germano et al., 199; Lilly et al. 1992)
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : TNU

      IMPLICIT NONE

      IF(MSUB .EQ. 1) THEN
        IF (IREAD.EQ.0 .AND. NTIME.EQ.1) THEN
          TNU = 0D0
        ELSE
          IF (ISGSMOD .EQ. 1) THEN
           !IF (IDYNMON.EQ.1) THEN
           !    IF (IFILTER .EQ. 1) CALL SGS_DSM_FZ
           !    IF (IFILTER .EQ. 3) CALL SGS_DSM_FXZ
           !    IF (IFILTER .EQ. 4) CALL SGS_DSM_FXZ_HZ
           !ELSE
           !   CALL SGS_CSM
           !ENDIF
          ELSE IF (ISGSMOD .EQ. 3) THEN
            IF (IDYNMON.EQ.1) THEN                  !* Dynamic global model
              IF (IFILTER .EQ. 1) CALL SGS_DVMG_FZ  !  Filtering in z
              IF (IFILTER .EQ. 2) CALL SGS_DVMG_FX  !            in x
              IF (IFILTER .EQ. 3) CALL SGS_DVMG_FXZ !            in x and z
            ELSE
              CALL SGS_CVM                          !* Constant Vreman model
            ENDIF
          ENDIF
        ENDIF

        CALL NUTZERO                                !* Set TNU to zero in the solid body
        CALL NUTINTERPOL                            !* Interpolation of TNU
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE NUTZEROREAD
!=======================================================================
!
!     Read TNU zero points
!     : In the solid body, sub-grid scale eddy viscosity, TNU, is set to
!       zero and these points are calculated during IBM pre-processing.
!
!     NNZ: Number of TNU zero point
!     INZ, JNZ, KNZ: Index of TNU zero point
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : INZ,JNZ,KNZ,NWALL_DVM

      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8 :: N

      OPEN(14,FILE='ibmpre_nutzero.bin',FORM='UNFORMATTED')
      READ(14) NNZ
      READ(14) (INZ(N),N=1,NNZ)
      READ(14) (JNZ(N),N=1,NNZ)
      READ(14) (KNZ(N),N=1,NNZ)
      CLOSE(14)

      OPEN(15,FILE='ibmpre_wallfdvm.bin',FORM='UNFORMATTED')
      READ(15) ((( NWALL_DVM(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      CLOSE(15)

      WRITE(*,30) NNZ
   30 FORMAT(' # of Nut=0  Pts = ',I10)

      RETURN
      END

!=======================================================================
      SUBROUTINE SGSFILTERINIT
!=======================================================================
!
!     Filter coefficients for the Simpson's rule
!     : For dynamic procedures to determine the model coefficient,
!       a test filter is applied based on the Simpson's rule.
!
!     CFX1, CFX2  : Filter coefficients in x-direction
!     CFZP1, CFZP2:                        z-
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFX1,CFX2,CFZP1,CFZP2

      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: H1,H2,H3,H4

      !* Initialization
      CFX1  = 0.
      CFX2  = 0.
      CFZP1 = 0.
      CFZP2 = 0.

      !* x-direction
      DO I=1,N1M
        H1=VDX(I)
        H2=VDX(I+1)
        CFX1(I,-1)=(2.*H1-H2)/6./H1
        CFX1(I, 0)=((H1+H2)**2)/6./H1/H2
        CFX1(I, 1)=(2.*H2-H1)/6./H2
      ENDDO
      DO I=3,N1M-2
        H1=VDX(I-1)
        H2=VDX(I  )
        H3=VDX(I+1)
        H4=VDX(I+2)
        CFX2(I,-2)=(2.*H1-H2)/6./H1/2.
        CFX2(I,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFX2(I, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFX2(I, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFX2(I, 2)=(2.*H4-H3)/6./H4/2.
      ENDDO
      H1=VDX(2)
      H2=VDX(3)
      CFX2(2,-1)=(2.*H1-H2)/6./H1
      CFX2(2, 0)=((H1+H2)**2)/6./H1/H2
      CFX2(2, 1)=(2.*H2-H1)/6./H2
      H1=VDX(N1M-1)
      H2=VDX(N1M)
      CFX2(N1M-1,-1)=(2.*H1-H2)/6./H1
      CFX2(N1M-1, 0)=((H1+H2)**2)/6./H1/H2
      CFX2(N1M-1, 1)=(2.*H2-H1)/6./H2

      !* z-direction
      IF (N3M.EQ.1) GOTO 101

      DO K=1,N3M
        H1=VDZ(K)
        H2=VDZ(K+1)
        CFZP1(K,-1)=(2.*H1-H2)/6./H1
        CFZP1(K, 0)=((H1+H2)**2)/6./H1/H2
        CFZP1(K, 1)=(2.*H2-H1)/6./H2
      ENDDO
      DO K=2,N3M-1
        H1=VDZ(K-1)
        H2=VDZ(K  )
        H3=VDZ(K+1)
        H4=VDZ(K+2)
        CFZP2(K,-2)=(2.*H1-H2)/6./H1/2.
        CFZP2(K,-1)=((H1+H2)**2)/6./H1/H2/2.
        CFZP2(K, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
        CFZP2(K, 1)=((H3+H4)**2)/6./H3/H4/2.
        CFZP2(K, 2)=(2.*H4-H3)/6./H4/2.
      ENDDO
      H1=VDZ(N3M)
      H2=VDZ(1)
      H3=VDZ(2)
      H4=VDZ(3)
      CFZP2(1,-2)=(2.*H1-H2)/6./H1/2.
      CFZP2(1,-1)=((H1+H2)**2)/6./H1/H2/2.
      CFZP2(1, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
      CFZP2(1, 1)=((H3+H4)**2)/6./H3/H4/2.
      CFZP2(1, 2)=(2.*H4-H3)/6./H4/2.
      H1=VDZ(N3M-1)
      H2=VDZ(N3M)
      H3=VDZ(1)
      H4=VDZ(2)
      CFZP2(N3M,-2)=(2.*H1-H2)/6./H1/2.
      CFZP2(N3M,-1)=((H1+H2)**2)/6./H1/H2/2.
      CFZP2(N3M, 0)=((2.*H2-H1)/6./H2+(2.*H3-H4)/6./H3)/2.
      CFZP2(N3M, 1)=((H3+H4)**2)/6./H3/H4/2.
      CFZP2(N3M, 2)=(2.*H4-H3)/6./H4/2.

  101 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE SGSINIT
!=======================================================================
!
!     Initialization of variables related to sub-grid scale stress tensors
!       TNU : sub-grid scale eddy viscosity at the cell center
!       TNUX,TNUY,TNUZ: Interpolated eddy viscosity for x-,y- and z-mtm
!       NWALL_DVM: 0 - inside or around the IB body, 1 - fluid (no body)
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : TNU,TNUX,TNUY,TNUZ,NWALL_DVM

      IMPLICIT NONE

      TNU  = 0.
      TNUX = 0.
      TNUY = 0.
      TNUZ = 0.
      IF ((ILES.EQ.1).AND.(IBMON.EQ.0)) NWALL_DVM = 1

      RETURN
      END

!=======================================================================
      SUBROUTINE SGS_CVM !* Constant Vreman model
!=======================================================================
!
!     CALCULATE SUB-GRID SCALE VISCOSITY TNU (NUT->INTEGER VARIABLE)
!       Vreman, 2004, An eddy-viscosity subgrid-scale model for turbulent
!         shear flow Algebraic theory and applications, Phys. Fluids, 16, pp. 3670-3681.
!
!     TNU=CSGS*SQRT(BBB/AAA)
!     CSGS   : VREMAN CONSTANT
!     BBB    : B_{11}*B_{22}-B_{12}*B_{12}+B_{11}*B_{33}
!             -B_{13}*B_{13}+B_{22}*B_{33}-B_{23}*B_{23}
!     AAA    : A_{ij}*A_{ij}
!     B_{ij} : SIGMA(DEL_{m}^2*A_{mi}*A_{mj}) where m=1, 2, 3
!     A_{ij} : DUj/DXi
!     DEL_{m}: FILTER SIZE (IN OTHER VIEW, GRID SIZE)
!     VG     : Velocity Gradiant
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,TNU
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,KP,KM,JP,JM,IP,IM
      REAL*8    :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8    :: UP,UM,VP,VM,WP,WM
      REAL*8    :: A(3,3),B(6)
      REAL*8    :: AAA,BBB

!$OMP PARALLEL DO private(I,J,KP,KM,JP,JM,IP,IM) &
!$OMP private(UP,UM,VP,VM,WP,WM) &
!$OMP private(VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33) &
!$OMP private(A,B,AAA,BBB)
      DO 10 K=1,N3M
      KP=KPV(K)
      KM=KMV(K)
      DO 10 J=1,N2M
      JP=JPV(J)
      JM=JMV(J)
      DO 10 I=1,N1M
      IP=IPV(I)
      IM=IMV(I)

      VG11=SSDX(I)*(U(IP,J,K)-U(I,J,K))

      UP=VVDY(JP)*0.25                                                &
        *(SDY(JP)*(U(I,J,K)+U(IP,J,K))+SDY(J)*(U(I,JP,K)+U(IP,JP,K))) &
        *(1.-FIXJU(J))                                                &
        +0.5*(U(I,N2,K)+U(IP,N2,K))*FIXJU(J)
      UM=VVDY(J)*0.25                                                 &
        *(SDY(J)*(U(I,JM,K)+U(IP,JM,K))+SDY(JM)*(U(I,J,K)+U(IP,J,K))) &
        *(1.-FIXJL(J))                                                &
        +0.5*(U(I,0,K)+U(IP,0,K))*FIXJL(J)
      VG12=SSDY(J)*(UP-UM)

      UP=VVDZ(KP)*0.25                                                 &
         *(SDZ(KP)*(U(I,J,K)+U(IP,J,K))+SDZ(K)*(U(I,J,KP)+U(IP,J,KP))) &
         *(1.-FIXKU(K))                                                &
         +0.5*(U(I,J,N3)+U(IP,J,N3))*FIXKU(K)
      UM=VVDZ(K)*0.25                                                  &
         *(SDZ(K)*(U(I,J,KM)+U(IP,J,KM))+SDZ(KM)*(U(I,J,K)+U(IP,J,K))) &
         *(1.-FIXKL(K))                                                &
         +0.5*(U(I,J,0)+U(IP,J,0))*FIXKL(K)
      VG13=SSDZ(K)*(UP-UM)

      VP=VVDX(IP)*0.25                                                 &
         *(SDX(IP)*(V(I,J,K)+V(I,JP,K))+SDX(I)*(V(IP,J,K)+V(IP,JP,K))) &
         *(1.-FIXIU(I))                                                &
         +0.5*(V(N1,J,K)+V(N1,JP,K))*FIXIU(I)
      VM=VVDX(I)*0.25                                                  &
         *(SDX(I)*(V(IM,J,K)+V(IM,JP,K))+SDX(IM)*(V(I,J,K)+V(I,JP,K))) &
         *(1.-FIXIL(I))                                                &
         +0.5*(V(0,J,K)+V(0,JP,K))*FIXIL(I)
      VG21=SSDX(I)*(VP-VM)

      VG22=SSDY(J)*(V(I,JP,K)-V(I,J,K))

      VP=VVDZ(KP)*0.25                                                 &
         *(SDZ(KP)*(V(I,J,K)+V(I,JP,K))+SDZ(K)*(V(I,J,KP)+V(I,JP,KP))) &
         *(1.-FIXKU(K))                                                &
         +0.5*(V(I,J,N3)+V(I,JP,N3))*FIXKU(K)
      VM=VVDZ(K)*0.25                                                  &
         *(SDZ(K)*(V(I,J,KM)+V(I,JP,KM))+SDZ(KM)*(V(I,J,K)+V(I,JP,K))) &
         *(1.-FIXKL(K))                                                &
         +0.5*(V(I,J,0)+V(I,JP,0))*FIXKL(K)
      VG23=SSDZ(K)*(VP-VM)

      WP=VVDX(IP)*0.25                                                 &
         *(SDX(IP)*(W(I,J,K)+W(I,J,KP))+SDX(I)*(W(IP,J,K)+W(IP,J,KP))) &
         *(1.-FIXIU(I))                                                &
         +0.5*(W(N1,J,K)+W(N1,J,KP))*FIXIU(I)
      WM=VVDX(I)*0.25                                                  &
         *(SDX(I)*(W(IM,J,K)+W(IM,J,KP))+SDX(IM)*(W(I,J,K)+W(I,J,KP))) &
         *(1.-FIXIL(I))                                                &
         +0.5*(W(0,J,K)+W(0,J,KP))*FIXIL(I)
      VG31=SSDX(I)*(WP-WM)

      WP=VVDY(JP)*0.25                                                 &
         *(SDY(JP)*(W(I,J,K)+W(I,J,KP))+SDY(J)*(W(I,JP,K)+W(I,JP,KP))) &
         *(1.-FIXJU(J))                                                &
         +0.5*(W(I,N2,K)+W(I,N2,KP))*FIXJU(J)
      WM=VVDY(J)*0.25                                                  &
         *(SDY(J)*(W(I,JM,K)+W(I,JM,KP))+SDY(JM)*(W(I,J,K)+W(I,J,KP))) &
         *(1.-FIXJL(J))                                                &
         +0.5*(W(I,0,K)+W(I,0,KP))*FIXJL(J)
      VG32=SSDY(J)*(WP-WM)

      VG33=SSDZ(K)*(W(I,J,KP)-W(I,J,K))

      A(1,1)=VG11
      A(1,2)=VG21
      A(1,3)=VG31
      A(2,1)=VG12
      A(2,2)=VG22
      A(2,3)=VG32
      A(3,1)=VG13
      A(3,2)=VG23
      A(3,3)=VG33

      B(1)=SDX(I)**2.*A(1,1)*A(1,1)  &
          +SDY(J)**2.*A(2,1)*A(2,1)  &
          +SDZ(K)**2.*A(3,1)*A(3,1)
      B(2)=SDX(I)**2.*A(1,1)*A(1,2)  &
          +SDY(J)**2.*A(2,1)*A(2,2)  &
          +SDZ(K)**2.*A(3,1)*A(3,2)
      B(3)=SDX(I)**2.*A(1,1)*A(1,3)  &
          +SDY(J)**2.*A(2,1)*A(2,3)  &
          +SDZ(K)**2.*A(3,1)*A(3,3)
      B(4)=SDX(I)**2.*A(1,2)*A(1,2)  &
          +SDY(J)**2.*A(2,2)*A(2,2)  &
          +SDZ(K)**2.*A(3,2)*A(3,2)
      B(5)=SDX(I)**2.*A(1,2)*A(1,3)  &
          +SDY(J)**2.*A(2,2)*A(2,3)  &
          +SDZ(K)**2.*A(3,2)*A(3,3)
      B(6)=SDX(I)**2.*A(1,3)*A(1,3)  &
          +SDY(J)**2.*A(2,3)*A(2,3)  &
          +SDZ(K)**2.*A(3,3)*A(3,3)

      AAA=A(1,1)**2+A(1,2)**2+A(1,3)**2+A(2,1)**2+A(2,2)**2+A(2,3)**2+A(3,1)**2+A(3,2)**2+A(3,3)**2
      BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)-B(2)**2  -B(3)**2  -B(5)**2

      BBB = MAX(BBB,0D0)
      IF (AAA.GT.TINY(1D0)) THEN
        TNU(I,J,K)= CSGS*SQRT(BBB/AAA)
      ELSE
        TNU(I,J,K)= 0D0 ! Uniform flow has zero eddy viscosity
      ENDIF

   10 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE SGS_DVMG_FXZ !* Dymanic global model with filtering in x- and z-
!=======================================================================
!
!     Calculate sub-grid scale viscosity, TNU
!       Park, Lee, Lee & Choi, 2006, A dynamic subgrid-scale eddy viscosity model with a global
!         model coefficient, Phys. Fluids, 18, 125109.
!       Lee, Choi & Park, 2010, Dynamic global model for large eddy simulation of transient
!         flow, Phys. Fluids, 22, 075106.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,TNU,NWALL_DVM
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8  :: LIJMIJ_V,MIJMIJ_V
      REAL*8  :: SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME
      REAL*8, DIMENSION (:),     ALLOCATABLE :: B,STR
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SR,ALP
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: LIJ,MIJ,UI

      LIJMIJ_V = 0.
      MIJMIJ_V = 0.
      CSGS     = 0.
      CSGSPS   = 0.

!$OMP PARALLEL &
!$OMP private(B,STR,SR,ALP,LIJ,MIJ,UI)
      allocate (B(6),STR(6))
      allocate (SR(M1M,M3M,6),ALP(M1M,M3M,9))
      allocate (LIJ(M1M,M3M,6),MIJ(M1M,M3M,6),UI(M1M,M3M,3))
!$OMP DO &
!$OMP private(SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME) &
!$OMP reduction(+:LIJMIJ_V) &
!$OMP reduction(+:MIJMIJ_V)
      DO J=1,N2M

        CALL STRAIN(SR,ALP,J)

        DO K=1,N3M
        DO I=1,N1M
          SDXF2=SDX(I)**2.
          SDYF2=SDY(J)**2.
          SDZF2=SDZ(K)**2.
          B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)  &
              +SDYF2*ALP(I,K,4)*ALP(I,K,4)  &
              +SDZF2*ALP(I,K,7)*ALP(I,K,7)
          B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)  &
              +SDYF2*ALP(I,K,4)*ALP(I,K,5)  &
              +SDZF2*ALP(I,K,7)*ALP(I,K,8)
          B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)  &
              +SDYF2*ALP(I,K,4)*ALP(I,K,6)  &
              +SDZF2*ALP(I,K,7)*ALP(I,K,9)
          B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)  &
              +SDYF2*ALP(I,K,5)*ALP(I,K,5)  &
              +SDZF2*ALP(I,K,8)*ALP(I,K,8)
          B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)  &
              +SDYF2*ALP(I,K,5)*ALP(I,K,6)  &
              +SDZF2*ALP(I,K,8)*ALP(I,K,9)
          B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)  &
              +SDYF2*ALP(I,K,6)*ALP(I,K,6)  &
              +SDZF2*ALP(I,K,9)*ALP(I,K,9)
          AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2.  &
             +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2.  &
             +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
          BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)  &
             -B(2)**2. -B(3)**2. -B(5)**2.
          BBB = AMAX1(BBB,1.0E-16)
          TNU(I,J,K) = SQRT(BBB/AAA)
          MIJ(I,K,1) = TNU(I,J,K)*SR(I,K,1)
          MIJ(I,K,2) = TNU(I,J,K)*SR(I,K,2)
          MIJ(I,K,3) = TNU(I,J,K)*SR(I,K,3)
          MIJ(I,K,4) = TNU(I,J,K)*SR(I,K,4)
          MIJ(I,K,5) = TNU(I,J,K)*SR(I,K,5)
          MIJ(I,K,6) = TNU(I,J,K)*SR(I,K,6)
        ENDDO
        ENDDO

        CALL TEST_FILTERXZ(MIJ(1,1,1))
        CALL TEST_FILTERXZ(MIJ(1,1,2))
        CALL TEST_FILTERXZ(MIJ(1,1,3))
        CALL TEST_FILTERXZ(MIJ(1,1,4))
        CALL TEST_FILTERXZ(MIJ(1,1,5))
        CALL TEST_FILTERXZ(MIJ(1,1,6))
        CALL TEST_FILTERXZ(ALP(1,1,1))
        CALL TEST_FILTERXZ(ALP(1,1,2))
        CALL TEST_FILTERXZ(ALP(1,1,3))
        CALL TEST_FILTERXZ(ALP(1,1,4))
        CALL TEST_FILTERXZ(ALP(1,1,5))
        CALL TEST_FILTERXZ(ALP(1,1,6))
        CALL TEST_FILTERXZ(ALP(1,1,7))
        CALL TEST_FILTERXZ(ALP(1,1,8))
        CALL TEST_FILTERXZ(ALP(1,1,9))

        DO K=1,N3M
        DO I=1,N1M
          SDXF2=(0.5*SDX(I-1)*(1.-FIXIU(I))+SDX(I)  &
                +0.5*SDX(I+1)*(1.-FIXIL(I)))**2.
          SDYF2=SDY(J)**2.
          SDZF2=(0.5*SDZ(K-1)*(1.-FIXKU(K))+SDZ(K)  &
                +0.5*SDZ(K+1)*(1.-FIXKL(K)))**2.
          B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)&
              +SDYF2*ALP(I,K,4)*ALP(I,K,4)&
              +SDZF2*ALP(I,K,7)*ALP(I,K,7)
          B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)&
              +SDYF2*ALP(I,K,4)*ALP(I,K,5)&
              +SDZF2*ALP(I,K,7)*ALP(I,K,8)
          B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)&
              +SDYF2*ALP(I,K,4)*ALP(I,K,6)&
              +SDZF2*ALP(I,K,7)*ALP(I,K,9)
          B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)&
              +SDYF2*ALP(I,K,5)*ALP(I,K,5)&
              +SDZF2*ALP(I,K,8)*ALP(I,K,8)
          B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)&
              +SDYF2*ALP(I,K,5)*ALP(I,K,6)&
              +SDZF2*ALP(I,K,8)*ALP(I,K,9)
          B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)&
              +SDYF2*ALP(I,K,6)*ALP(I,K,6)&
              +SDZF2*ALP(I,K,9)*ALP(I,K,9)
          STR(1)=ALP(I,K,1)
          STR(2)=0.5*(ALP(I,K,2)+ALP(I,K,4))
          STR(3)=0.5*(ALP(I,K,3)+ALP(I,K,7))
          STR(4)=ALP(I,K,5)
          STR(5)=0.5*(ALP(I,K,6)+ALP(I,K,8))
          STR(6)=ALP(I,K,9)
          AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2. &
             +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2. &
             +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
          BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)  &
             -B(2)**2. -B(3)**2. -B(5)**2.
          BBB = AMAX1(BBB,1.0E-16)
          AMI = SQRT(BBB/AAA)

          MIJ(I,K,1)=AMI*STR(1)-MIJ(I,K,1)
          MIJ(I,K,2)=AMI*STR(2)-MIJ(I,K,2)
          MIJ(I,K,3)=AMI*STR(3)-MIJ(I,K,3)
          MIJ(I,K,4)=AMI*STR(4)-MIJ(I,K,4)
          MIJ(I,K,5)=AMI*STR(5)-MIJ(I,K,5)
          MIJ(I,K,6)=AMI*STR(6)-MIJ(I,K,6)
        ENDDO
        ENDDO

        DO K=1,N3M
        DO I=1,N1M
          UI (I,K,1)=0.5*(U(I,J,K)+U(I+1,J,K))
          UI (I,K,2)=0.5*(V(I,J,K)+V(I,J+1,K))
          UI (I,K,3)=0.5*(W(I,J,K)+W(I,J,K+1))
          LIJ(I,K,1)=UI(I,K,1)*UI(I,K,1)
          LIJ(I,K,2)=UI(I,K,1)*UI(I,K,2)
          LIJ(I,K,3)=UI(I,K,1)*UI(I,K,3)
          LIJ(I,K,4)=UI(I,K,2)*UI(I,K,2)
          LIJ(I,K,5)=UI(I,K,2)*UI(I,K,3)
          LIJ(I,K,6)=UI(I,K,3)*UI(I,K,3)
        ENDDO
        ENDDO

        CALL TEST_FILTERXZ(UI(1,1,1))
        CALL TEST_FILTERXZ(UI(1,1,2))
        CALL TEST_FILTERXZ(UI(1,1,3))
        CALL TEST_FILTERXZ(LIJ(1,1,1))
        CALL TEST_FILTERXZ(LIJ(1,1,2))
        CALL TEST_FILTERXZ(LIJ(1,1,3))
        CALL TEST_FILTERXZ(LIJ(1,1,4))
        CALL TEST_FILTERXZ(LIJ(1,1,5))
        CALL TEST_FILTERXZ(LIJ(1,1,6))

        DO K=1,N3M
        DO I=1,N1M
          VOLUME=SDX(I)*SDY(J)*SDZ(K)
          LIJ(I,K,1)=LIJ(I,K,1)-UI(I,K,1)*UI(I,K,1)
          LIJ(I,K,2)=LIJ(I,K,2)-UI(I,K,1)*UI(I,K,2)
          LIJ(I,K,3)=LIJ(I,K,3)-UI(I,K,1)*UI(I,K,3)
          LIJ(I,K,4)=LIJ(I,K,4)-UI(I,K,2)*UI(I,K,2)
          LIJ(I,K,5)=LIJ(I,K,5)-UI(I,K,2)*UI(I,K,3)
          LIJ(I,K,6)=LIJ(I,K,6)-UI(I,K,3)*UI(I,K,3)
          LIJMIJ_V=LIJMIJ_V                                     &
               +(2.*LIJ(I,K,2)*MIJ(I,K,2)+LIJ(I,K,1)*MIJ(I,K,1) &
                +2.*LIJ(I,K,3)*MIJ(I,K,3)+LIJ(I,K,4)*MIJ(I,K,4) &
                +2.*LIJ(I,K,5)*MIJ(I,K,5)+LIJ(I,K,6)*MIJ(I,K,6))&
               *FLOAT(NWALL_DVM(I,J,K))*VOLUME
          MIJMIJ_V=MIJMIJ_V                                     &
                 +(2.*MIJ(I,K,2)**2.+MIJ(I,K,1)**2.             &
                  +2.*MIJ(I,K,3)**2.+MIJ(I,K,4)**2.             &
                  +2.*MIJ(I,K,5)**2.+MIJ(I,K,6)**2.)            &
                  *FLOAT(NWALL_DVM(I,J,K))*VOLUME
        ENDDO
        ENDDO

      ENDDO
!$OMP END DO
      DEALLOCATE (B,STR)
      DEALLOCATE (SR,ALP)
      DEALLOCATE (LIJ,MIJ,UI)
!$OMP END PARALLEL

      CSGS = AMAX1(0.,-0.5*LIJMIJ_V/MIJMIJ_V)

!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        TNU(I,J,K)=CSGS*TNU(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE SGS_DVMG_FX !* Dymanic global model with filtering in x-
!=======================================================================
!
!     Calculate sub-grid scale viscosity, TNU
!       Park, Lee, Lee & Choi, 2006, A dynamic subgrid-scale eddy viscosity model with a global
!         model coefficient, Phys. Fluids, 18, 125109.
!       Lee, Choi & Park, 2010, Dynamic global model for large eddy simulation of transient
!         flow, Phys. Fluids, 22, 075106.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,TNU,NWALL_DVM
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: LIJMIJ_V,MIJMIJ_V
      REAL*8    :: SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME
      REAL*8, DIMENSION (:),     ALLOCATABLE :: B,STR
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SR,ALP
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: LIJ,MIJ,UI

      LIJMIJ_V = 0.
      MIJMIJ_V = 0.
      CSGS     = 0.
      CSGSPS   = 0.

!$OMP PARALLEL &
!$OMP private(B,STR,SR,ALP,LIJ,MIJ,UI)
      allocate (B(6),STR(6))
      allocate (SR(M1M,M3M,6),ALP(M1M,M3M,9))
      allocate (LIJ(M1M,M3M,6),MIJ(M1M,M3M,6),UI(M1M,M3M,3))
!$OMP DO &
!$OMP private(SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME) &
!$OMP reduction(+:LIJMIJ_V) &
!$OMP reduction(+:MIJMIJ_V)
      DO J=1,N2M

      CALL STRAIN(SR,ALP,J)

      DO K=1,N3M
      DO I=1,N1M
      SDXF2=SDX(I)**2.
      SDYF2=SDY(J)**2.
      SDZF2=SDZ(K)**2.
      B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,4)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,7)
      B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,5)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,8)
      B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,9)
      B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)  &
          +SDYF2*ALP(I,K,5)*ALP(I,K,5)  &
          +SDZF2*ALP(I,K,8)*ALP(I,K,8)
      B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,5)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,8)*ALP(I,K,9)
      B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,6)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,9)*ALP(I,K,9)
      AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2.  &
         +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2.  &
         +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
      BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)-B(2)**2. -B(3)**2. -B(5)**2.
      BBB=AMAX1(BBB,1.0E-16)
      TNU(I,J,K)= SQRT(BBB/AAA)
      MIJ(I,K,1)= TNU(I,J,K)*SR(I,K,1)
      MIJ(I,K,2)= TNU(I,J,K)*SR(I,K,2)
      MIJ(I,K,3)= TNU(I,J,K)*SR(I,K,3)
      MIJ(I,K,4)= TNU(I,J,K)*SR(I,K,4)
      MIJ(I,K,5)= TNU(I,J,K)*SR(I,K,5)
      MIJ(I,K,6)= TNU(I,J,K)*SR(I,K,6)
      ENDDO
      ENDDO

      CALL TEST_FILTERX(MIJ(1,1,1))
      CALL TEST_FILTERX(MIJ(1,1,2))
      CALL TEST_FILTERX(MIJ(1,1,3))
      CALL TEST_FILTERX(MIJ(1,1,4))
      CALL TEST_FILTERX(MIJ(1,1,5))
      CALL TEST_FILTERX(MIJ(1,1,6))
      CALL TEST_FILTERX(ALP(1,1,1))
      CALL TEST_FILTERX(ALP(1,1,2))
      CALL TEST_FILTERX(ALP(1,1,3))
      CALL TEST_FILTERX(ALP(1,1,4))
      CALL TEST_FILTERX(ALP(1,1,5))
      CALL TEST_FILTERX(ALP(1,1,6))
      CALL TEST_FILTERX(ALP(1,1,7))
      CALL TEST_FILTERX(ALP(1,1,8))
      CALL TEST_FILTERX(ALP(1,1,9))

      DO K=1,N3M
      DO I=1,N1M
      SDXF2=(0.5*SDX(I-1)*(1.-FIXIU(I))+SDX(I)+0.5*SDX(I+1)*(1.-FIXIL(I)))**2.
      SDYF2=SDY(J)**2.
      SDZF2=SDZ(K)**2.
      B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,4)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,7)
      B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,5)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,8)
      B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,4)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,7)*ALP(I,K,9)
      B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)  &
          +SDYF2*ALP(I,K,5)*ALP(I,K,5)  &
          +SDZF2*ALP(I,K,8)*ALP(I,K,8)
      B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,5)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,8)*ALP(I,K,9)
      B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)  &
          +SDYF2*ALP(I,K,6)*ALP(I,K,6)  &
          +SDZF2*ALP(I,K,9)*ALP(I,K,9)
      STR(1)= ALP(I,K,1)
      STR(2)= 0.5*(ALP(I,K,2)+ALP(I,K,4))
      STR(3)= 0.5*(ALP(I,K,3)+ALP(I,K,7))
      STR(4)= ALP(I,K,5)
      STR(5)= 0.5*(ALP(I,K,6)+ALP(I,K,8))
      STR(6)= ALP(I,K,9)
      AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2.  &
         +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2.  &
         +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
      BBB=B(1)*B(4)+B(4)*B(6)+B(6)*B(1)-B(2)**2. -B(3)**2. -B(5)**2.
      BBB=AMAX1(BBB,1.0E-16)
      AMI=SQRT(BBB/AAA)
      MIJ(I,K,1)=AMI*STR(1)-MIJ(I,K,1)
      MIJ(I,K,2)=AMI*STR(2)-MIJ(I,K,2)
      MIJ(I,K,3)=AMI*STR(3)-MIJ(I,K,3)
      MIJ(I,K,4)=AMI*STR(4)-MIJ(I,K,4)
      MIJ(I,K,5)=AMI*STR(5)-MIJ(I,K,5)
      MIJ(I,K,6)=AMI*STR(6)-MIJ(I,K,6)
      ENDDO
      ENDDO

      DO K=1,N3M
      DO I=1,N1M
      UI (I,K,1)= 0.5*(U(I,J,K)+U(I+1,J,K))
      UI (I,K,2)= 0.5*(V(I,J,K)+V(I,J+1,K))
      UI (I,K,3)= 0.5*(W(I,J,K)+W(I,J,K+1))
      LIJ(I,K,1)= UI(I,K,1)*UI(I,K,1)
      LIJ(I,K,2)= UI(I,K,1)*UI(I,K,2)
      LIJ(I,K,3)= UI(I,K,1)*UI(I,K,3)
      LIJ(I,K,4)= UI(I,K,2)*UI(I,K,2)
      LIJ(I,K,5)= UI(I,K,2)*UI(I,K,3)
      LIJ(I,K,6)= UI(I,K,3)*UI(I,K,3)
      ENDDO
      ENDDO

      CALL TEST_FILTERX(UI(1,1,1))
      CALL TEST_FILTERX(UI(1,1,2))
      CALL TEST_FILTERX(UI(1,1,3))
      CALL TEST_FILTERX(LIJ(1,1,1))
      CALL TEST_FILTERX(LIJ(1,1,2))
      CALL TEST_FILTERX(LIJ(1,1,3))
      CALL TEST_FILTERX(LIJ(1,1,4))
      CALL TEST_FILTERX(LIJ(1,1,5))
      CALL TEST_FILTERX(LIJ(1,1,6))

      DO K=1,N3M
      DO I=1,N1M
      VOLUME=SDX(I)*SDY(J)*SDZ(K)
      LIJ(I,K,1)=LIJ(I,K,1)-UI(I,K,1)*UI(I,K,1)
      LIJ(I,K,2)=LIJ(I,K,2)-UI(I,K,1)*UI(I,K,2)
      LIJ(I,K,3)=LIJ(I,K,3)-UI(I,K,1)*UI(I,K,3)
      LIJ(I,K,4)=LIJ(I,K,4)-UI(I,K,2)*UI(I,K,2)
      LIJ(I,K,5)=LIJ(I,K,5)-UI(I,K,2)*UI(I,K,3)
      LIJ(I,K,6)=LIJ(I,K,6)-UI(I,K,3)*UI(I,K,3)
      LIJMIJ_V=LIJMIJ_V                                         &
              +(2.*LIJ(I,K,2)*MIJ(I,K,2)+LIJ(I,K,1)*MIJ(I,K,1)  &
               +2.*LIJ(I,K,3)*MIJ(I,K,3)+LIJ(I,K,4)*MIJ(I,K,4)  &
               +2.*LIJ(I,K,5)*MIJ(I,K,5)+LIJ(I,K,6)*MIJ(I,K,6)) &
               *FLOAT(NWALL_DVM(I,J,K))*VOLUME
      MIJMIJ_V=MIJMIJ_V                                 &
              +(2.*MIJ(I,K,2)**2.+MIJ(I,K,1)**2.        &
               +2.*MIJ(I,K,3)**2.+MIJ(I,K,4)**2.        &
               +2.*MIJ(I,K,5)**2.+MIJ(I,K,6)**2.)       &
               *FLOAT(NWALL_DVM(I,J,K))*VOLUME
      ENDDO
      ENDDO

      ENDDO
!$OMP END DO
      DEALLOCATE (B,STR)
      DEALLOCATE (SR,ALP)
      DEALLOCATE (LIJ,MIJ,UI)
!$OMP END PARALLEL

      CSGS=AMAX1(0.,-0.5*LIJMIJ_V/MIJMIJ_V)

!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
         TNU(I,J,K)=CSGS*TNU(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE SGS_DVMG_FZ !* Dymanic global model with filtering in z-
!=======================================================================
!
!     Calculate sub-grid scale viscosity, TNU
!       Park, Lee, Lee & Choi, 2006, A dynamic subgrid-scale eddy viscosity model with a global
!         model coefficient, Phys. Fluids, 18, 125109.
!       Lee, Choi & Park, 2010, Dynamic global model for large eddy simulation of transient
!         flow, Phys. Fluids, 22, 075106.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,TNU,NWALL_DVM
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: LIJMIJ_V,MIJMIJ_V
      REAL*8    :: SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME
      REAL*8, DIMENSION (:),     ALLOCATABLE :: B,STR
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: SR,ALP
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: LIJ,MIJ,UI

      LIJMIJ_V = 0.
      MIJMIJ_V = 0.
      CSGS  = 0.
      CSGSPS= 0.

!$OMP PARALLEL &
!$OMP private(B,STR,SR,ALP,LIJ,MIJ,UI)
      allocate (B(6),STR(6))
      allocate (SR(M1M,M3M,6),ALP(M1M,M3M,9))
      allocate (LIJ(M1M,M3M,6),MIJ(M1M,M3M,6),UI(M1M,M3M,3))
!$OMP DO &
!$OMP private(SDXF2,SDYF2,SDZF2,AAA,BBB,AMI,VOLUME) &
!$OMP reduction(+:LIJMIJ_V) &
!$OMP reduction(+:MIJMIJ_V)
      DO J=1,N2M

        CALL STRAIN(SR,ALP,J)

        DO K=1,N3M
        DO I=1,N1M
        SDXF2=SDX(I)**2.
        SDYF2=SDY(J)**2.
        SDZF2=SDZ(K)**2.
        B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,4)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,7)
        B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,5)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,8)
        B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,9)
        B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)    &
            +SDYF2*ALP(I,K,5)*ALP(I,K,5)    &
            +SDZF2*ALP(I,K,8)*ALP(I,K,8)
        B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,5)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,8)*ALP(I,K,9)
        B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,6)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,9)*ALP(I,K,9)
        AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2.    &
           +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2.    &
           +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
        BBB= B(1)*B(4)+B(4)*B(6)+B(6)*B(1)-B(2)**2. -B(3)**2. -B(5)**2.
        BBB= AMAX1(BBB,1.0E-16)
        TNU(I,J,K)=SQRT(BBB/AAA)
        MIJ(I,K,1)=TNU(I,J,K)*SR(I,K,1)
        MIJ(I,K,2)=TNU(I,J,K)*SR(I,K,2)
        MIJ(I,K,3)=TNU(I,J,K)*SR(I,K,3)
        MIJ(I,K,4)=TNU(I,J,K)*SR(I,K,4)
        MIJ(I,K,5)=TNU(I,J,K)*SR(I,K,5)
        MIJ(I,K,6)=TNU(I,J,K)*SR(I,K,6)
        ENDDO
        ENDDO

        CALL TEST_FILTERZ(MIJ(1,1,1))
        CALL TEST_FILTERZ(MIJ(1,1,2))
        CALL TEST_FILTERZ(MIJ(1,1,3))
        CALL TEST_FILTERZ(MIJ(1,1,4))
        CALL TEST_FILTERZ(MIJ(1,1,5))
        CALL TEST_FILTERZ(MIJ(1,1,6))
        CALL TEST_FILTERZ(ALP(1,1,1))
        CALL TEST_FILTERZ(ALP(1,1,2))
        CALL TEST_FILTERZ(ALP(1,1,3))
        CALL TEST_FILTERZ(ALP(1,1,4))
        CALL TEST_FILTERZ(ALP(1,1,5))
        CALL TEST_FILTERZ(ALP(1,1,6))
        CALL TEST_FILTERZ(ALP(1,1,7))
        CALL TEST_FILTERZ(ALP(1,1,8))
        CALL TEST_FILTERZ(ALP(1,1,9))

        DO K=1,N3M
        DO I=1,N1M
        SDXF2=SDX(I)**2.
        SDYF2=SDY(J)**2.
        SDZF2=(0.5*SDZ(K-1)*(1.-FIXKU(K))+SDZ(K)+0.5*SDZ(K+1)*(1.-FIXKL(K)))**2.
        B(1)=SDXF2*ALP(I,K,1)*ALP(I,K,1)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,4)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,7)
        B(2)=SDXF2*ALP(I,K,1)*ALP(I,K,2)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,5)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,8)
        B(3)=SDXF2*ALP(I,K,1)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,4)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,7)*ALP(I,K,9)
        B(4)=SDXF2*ALP(I,K,2)*ALP(I,K,2)    &
            +SDYF2*ALP(I,K,5)*ALP(I,K,5)    &
            +SDZF2*ALP(I,K,8)*ALP(I,K,8)
        B(5)=SDXF2*ALP(I,K,2)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,5)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,8)*ALP(I,K,9)
        B(6)=SDXF2*ALP(I,K,3)*ALP(I,K,3)    &
            +SDYF2*ALP(I,K,6)*ALP(I,K,6)    &
            +SDZF2*ALP(I,K,9)*ALP(I,K,9)
        STR(1)=ALP(I,K,1)
        STR(2)=0.5*(ALP(I,K,2)+ALP(I,K,4))
        STR(3)=0.5*(ALP(I,K,3)+ALP(I,K,7))
        STR(4)=ALP(I,K,5)
        STR(5)=0.5*(ALP(I,K,6)+ALP(I,K,8))
        STR(6)=ALP(I,K,9)
        AAA=ALP(I,K,1)**2.+ALP(I,K,2)**2.+ALP(I,K,3)**2.    &
           +ALP(I,K,4)**2.+ALP(I,K,5)**2.+ALP(I,K,6)**2.    &
           +ALP(I,K,7)**2.+ALP(I,K,8)**2.+ALP(I,K,9)**2.
        BBB= B(1)*B(4)+B(4)*B(6)+B(6)*B(1)-B(2)**2. -B(3)**2. -B(5)**2.
        BBB= AMAX1(BBB,1.0E-16)
        AMI= SQRT(BBB/AAA)
        MIJ(I,K,1)=AMI*STR(1)-MIJ(I,K,1)
        MIJ(I,K,2)=AMI*STR(2)-MIJ(I,K,2)
        MIJ(I,K,3)=AMI*STR(3)-MIJ(I,K,3)
        MIJ(I,K,4)=AMI*STR(4)-MIJ(I,K,4)
        MIJ(I,K,5)=AMI*STR(5)-MIJ(I,K,5)
        MIJ(I,K,6)=AMI*STR(6)-MIJ(I,K,6)
        ENDDO
        ENDDO

        DO K=1,N3M
        DO I=1,N1M
        UI (I,K,1)= 0.5*(U(I,J,K)+U(I+1,J,K))
        UI (I,K,2)= 0.5*(V(I,J,K)+V(I,J+1,K))
        UI (I,K,3)= 0.5*(W(I,J,K)+W(I,J,K+1))
        LIJ(I,K,1)= UI(I,K,1)*UI(I,K,1)
        LIJ(I,K,2)= UI(I,K,1)*UI(I,K,2)
        LIJ(I,K,3)= UI(I,K,1)*UI(I,K,3)
        LIJ(I,K,4)= UI(I,K,2)*UI(I,K,2)
        LIJ(I,K,5)= UI(I,K,2)*UI(I,K,3)
        LIJ(I,K,6)= UI(I,K,3)*UI(I,K,3)
        ENDDO
        ENDDO

        CALL TEST_FILTERZ(UI(1,1,1))
        CALL TEST_FILTERZ(UI(1,1,2))
        CALL TEST_FILTERZ(UI(1,1,3))
        CALL TEST_FILTERZ(LIJ(1,1,1))
        CALL TEST_FILTERZ(LIJ(1,1,2))
        CALL TEST_FILTERZ(LIJ(1,1,3))
        CALL TEST_FILTERZ(LIJ(1,1,4))
        CALL TEST_FILTERZ(LIJ(1,1,5))
        CALL TEST_FILTERZ(LIJ(1,1,6))

        DO K=1,N3M
        DO I=1,N1M
        VOLUME=SDX(I)*SDY(J)*SDZ(K)
        LIJ(I,K,1)=LIJ(I,K,1)-UI(I,K,1)*UI(I,K,1)
        LIJ(I,K,2)=LIJ(I,K,2)-UI(I,K,1)*UI(I,K,2)
        LIJ(I,K,3)=LIJ(I,K,3)-UI(I,K,1)*UI(I,K,3)
        LIJ(I,K,4)=LIJ(I,K,4)-UI(I,K,2)*UI(I,K,2)
        LIJ(I,K,5)=LIJ(I,K,5)-UI(I,K,2)*UI(I,K,3)
        LIJ(I,K,6)=LIJ(I,K,6)-UI(I,K,3)*UI(I,K,3)
        LIJMIJ_V=LIJMIJ_V                                           &
                +(2.*LIJ(I,K,2)*MIJ(I,K,2)+LIJ(I,K,1)*MIJ(I,K,1)    &
                 +2.*LIJ(I,K,3)*MIJ(I,K,3)+LIJ(I,K,4)*MIJ(I,K,4)    &
                 +2.*LIJ(I,K,5)*MIJ(I,K,5)+LIJ(I,K,6)*MIJ(I,K,6))   &
                 *FLOAT(NWALL_DVM(I,J,K))*VOLUME
        MIJMIJ_V=MIJMIJ_V                             &
                +(2.*MIJ(I,K,2)**2.+MIJ(I,K,1)**2.    &
                 +2.*MIJ(I,K,3)**2.+MIJ(I,K,4)**2.    &
                 +2.*MIJ(I,K,5)**2.+MIJ(I,K,6)**2.)   &
                 *FLOAT(NWALL_DVM(I,J,K))*VOLUME
        ENDDO
        ENDDO

      ENDDO
!$OMP END DO
      DEALLOCATE (B,STR)
      DEALLOCATE (SR,ALP)
      DEALLOCATE (LIJ,MIJ,UI)
!$OMP END PARALLEL

      CSGS= AMAX1(0.,-0.5*LIJMIJ_V/MIJMIJ_V)

!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        TNU(I,J,K)=CSGS*TNU(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE STRAIN(SR,ALP,J)
!=======================================================================
!
!     Calcultae velocity gradient and strain rate
!     SR : Strain rate = 0.5*(DUi/DXj+DUj/DXi)
!          SR(*,*,1)= VG11
!          SR(*,*,2)= 0.5*(VG12+VG21)
!          SR(*,*,3)= 0.5*(VG13+VG31)
!          SR(*,*,4)= VG22
!          SR(*,*,5)= 0.5*(VG23+VG32)
!          SR(*,*,6)= VG33
!     ALP: Velocity gradient = DUi/DXj
!          ALP(*,*,1)= VG11
!          ALP(*,*,2)= VG21
!          ALP(*,*,3)= VG31
!          ALP(*,*,4)= VG12
!          ALP(*,*,5)= VG22
!          ALP(*,*,6)= VG32
!          ALP(*,*,7)= VG13
!          ALP(*,*,8)= VG23
!          ALP(*,*,9)= VG33
!     J  : y-index
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8    :: SR(M1M,M3M,6),ALP(M1M,M3M,9)
      REAL*8    :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8    :: UP,UM,VP,VM,WP,WM

      DO K=1,N3M
      DO I=1,N1M

      VG11=SSDX(I)*(U(I+1,J,K)-U(I,J,K))

      UP=VVDY(J+1)*0.25                                     &
         *(SDY(J+1)*(U(I,J,K)+U(I+1,J,K))                   &
          +SDY(J)*(U(I,J+1,K)+U(I+1,J+1,K)))*(1.-FIXJU(J))  &
         +0.5*(U(I,N2,K)+U(I+1,N2,K))*FIXJU(J)
      UM=VVDY(J)*0.25                                       &
         *(SDY(J)*(U(I,J-1,K)+U(I+1,J-1,K))                 &
          +SDY(J-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXJL(J))    &
         +0.5*(U(I,0,K)+U(I+1,0,K))*FIXJL(J)
      VG12=SSDY(J)*(UP-UM)

      UP=VVDZ(K+1)*0.25                                     &
         *(SDZ(K+1)*(U(I,J,K)+U(I+1,J,K))                   &
          +SDZ(K)*(U(I,J,K+1)+U(I+1,J,K+1)))*(1.-FIXKU(K))  &
         +0.5*(U(I,J,N3)+U(I+1,J,N3))*FIXKU(K)
      UM=VVDZ(K)*0.25                                       &
         *(SDZ(K)*(U(I,J,K-1)+U(I+1,J,K-1))                 &
          +SDZ(K-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXKL(K))    &
         +0.5*(U(I,J,0)+U(I+1,J,0))*FIXKL(K)
      VG13=SSDZ(K)*(UP-UM)

      VP=VVDX(I+1)*0.25                                     &
         *(SDX(I+1)*(V(I,J,K)+V(I,J+1,K))                   &
          +SDX(I)*(V(I+1,J,K)+V(I+1,J+1,K)))*(1.-FIXIU(I))  &
         +0.5*(V(N1,J,K)+V(N1,J+1,K))*FIXIU(I)
      VM=VVDX(I)*0.25                                       &
         *(SDX(I)*(V(I-1,J,K)+V(I-1,J+1,K))                 &
          +SDX(I-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXIL(I))    &
         +0.5*(V(0,J,K)+V(0,J+1,K))*FIXIL(I)
      VG21=SSDX(I)*(VP-VM)

      VG22=SSDY(J)*(V(I,J+1,K)-V(I,J,K))

      VP=VVDZ(K+1)*0.25                                     &
         *(SDZ(K+1)*(V(I,J,K)+V(I,J+1,K))                   &
          +SDZ(K)*(V(I,J,K+1)+V(I,J+1,K+1)))*(1.-FIXKU(K))  &
         +0.5*(V(I,J,N3)+V(I,J+1,N3))*FIXKU(K)
      VM=VVDZ(K)*0.25                                       &
         *(SDZ(K)*(V(I,J,K-1)+V(I,J+1,K-1))                 &
          +SDZ(K-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXKL(K))    &
         +0.5*(V(I,J,0)+V(I,J+1,0))*FIXKL(K)
      VG23=SSDZ(K)*(VP-VM)

      WP=VVDX(I+1)*0.25                                     &
         *(SDX(I+1)*(W(I,J,K)+W(I,J,K+1))                   &
          +SDX(I)*(W(I+1,J,K)+W(I+1,J,K+1)))*(1.-FIXIU(I))  &
         +0.5*(W(N1,J,K)+W(N1,J,K+1))*FIXIU(I)
      WM=VVDX(I)*0.25                                       &
         *(SDX(I)*(W(I-1,J,K)+W(I-1,J,K+1))                 &
          +SDX(I-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXIL(I))    &
         +0.5*(W(0,J,K)+W(0,J,K+1))*FIXIL(I)
      VG31=SSDX(I)*(WP-WM)

      WP=VVDY(J+1)*0.25                                     &
         *(SDY(J+1)*(W(I,J,K)+W(I,J,K+1))                   &
          +SDY(J)*(W(I,J+1,K)+W(I,J+1,K+1)))*(1.-FIXJU(J))  &
         +0.5*(W(I,N2,K)+W(I,N2,K+1))*FIXJU(J)
      WM=VVDY(J)*0.25                                       &
         *(SDY(J)*(W(I,J-1,K)+W(I,J-1,K+1))                 &
          +SDY(J-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXJL(J))    &
         +0.5*(W(I,0,K)+W(I,0,K+1))*FIXJL(J)
      VG32=SSDY(J)*(WP-WM)

      VG33=SSDZ(K)*(W(I,J,K+1)-W(I,J,K))

      ALP(I,K,1)= VG11
      ALP(I,K,2)= VG21
      ALP(I,K,3)= VG31
      ALP(I,K,4)= VG12
      ALP(I,K,5)= VG22
      ALP(I,K,6)= VG32
      ALP(I,K,7)= VG13
      ALP(I,K,8)= VG23
      ALP(I,K,9)= VG33

      SR(I,K,1)= VG11
      SR(I,K,2)= 0.5*(VG12+VG21)
      SR(I,K,3)= 0.5*(VG13+VG31)
      SR(I,K,4)= VG22
      SR(I,K,5)= 0.5*(VG23+VG32)
      SR(I,K,6)= VG33

      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE TEST_FILTERXZ(A0)
!=======================================================================
!
!     Test filter in x- and z- for dynamic procedures
!
!     IPX,IPZ: 0 (not periodic), 1 (periodc)
!              These options are defined in the SIDEWALLBC subroutine.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFX1,CFZP1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: A0(M1M,M3M),TMP(M1M,M3M)

      IF (IPX.EQ.0) THEN
      DO K=1,N3M
      DO I=2,N1M-1
      TMP(I,K  )=CFX1(I,-1)*A0(I-1  ,K)  &
                +CFX1(I, 0)*A0(I    ,K)  &
                +CFX1(I, 1)*A0(I+1  ,K)
      ENDDO
      TMP(1  ,K)=A0(1  ,K)
      TMP(N1M,K)=A0(N1M,K)
      ENDDO
      ELSE IF (IPX.EQ.1) THEN
      DO K=1,N3M
      DO I=2,N1M-1
      TMP(I,K  )=CFX1(I,-1)*A0(I-1  ,K)  &
                +CFX1(I, 0)*A0(I    ,K)  &
                +CFX1(I, 1)*A0(I+1  ,K)
      ENDDO
      TMP(1  ,K)=CFX1(I,-1)*A0(N1M  ,K)  &
                +CFX1(I, 0)*A0(1    ,K)  &
                +CFX1(I, 1)*A0(2    ,K)
      TMP(N1M,K)=CFX1(I,-1)*A0(N1M-1,K)  &
                +CFX1(I, 0)*A0(N1M  ,K)  &
                +CFX1(I, 1)*A0(1    ,K)
      ENDDO
      ENDIF

      IF (IPZ.EQ.0) THEN
      DO I=1,N1M
      DO K=2,N3M-1
      A0(I,K  )=CFZP1(K  ,-1)*TMP(I,K-1)  &
               +CFZP1(K  , 0)*TMP(I,K  )  &
               +CFZP1(K  , 1)*TMP(I,K+1)
      ENDDO
      A0(I,1  )=TMP(I,1  )
      A0(I,N3M)=TMP(I,N3M)
      ENDDO
      ELSE IF (IPZ.EQ.1) THEN
      DO I=1,N1M
      DO K=2,N3M-1
      A0(I,K  )=CFZP1(K  ,-1)*TMP(I,K-1)  &
               +CFZP1(K  , 0)*TMP(I,K  )  &
               +CFZP1(K  , 1)*TMP(I,K+1)
      ENDDO
      A0(I,1  )=CFZP1(1  ,-1)*TMP(I,N3M)  &
               +CFZP1(1  , 0)*TMP(I,1  )  &
               +CFZP1(1  , 1)*TMP(I,2  )
      A0(I,N3M)=CFZP1(N3M,-1)*TMP(I,N3M-1) &
               +CFZP1(N3M, 0)*TMP(I,N3M  ) &
               +CFZP1(N3M, 1)*TMP(I,1    )
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE TEST_FILTERX(A0)
!=======================================================================
!
!     Test filter in x- for dynamic procedures
!
!     IPX,IPZ: 0 (not periodic), 1 (periodc)
!              These options are defined in the SIDEWALLBC subroutine.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFX1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    ::  A0(M1M,M3M),TMP(M1M,M3M)

      IF (IPX.EQ.0) THEN
      DO K=1,N3M
        DO I=2,N1M-1
        TMP(I,K  )=CFX1(I,-1)*A0(I-1  ,K)  &
                  +CFX1(I, 0)*A0(I    ,K)  &
                  +CFX1(I, 1)*A0(I+1  ,K)
        ENDDO
        TMP(1  ,K)=A0(1  ,K)
        TMP(N1M,K)=A0(N1M,K)
      ENDDO
      ELSE IF (IPX.EQ.1) THEN
        DO K=1,N3M
        DO I=2,N1M-1
        TMP(I,K  )=CFX1(I,-1)*A0(I-1  ,K) &
                  +CFX1(I, 0)*A0(I    ,K) &
                  +CFX1(I, 1)*A0(I+1  ,K)
        ENDDO
        TMP(1  ,K)=CFX1(I,-1)*A0(N1M  ,K)   &
                  +CFX1(I, 0)*A0(1    ,K)   &
                  +CFX1(I, 1)*A0(2    ,K)
        TMP(N1M,K)=CFX1(I,-1)*A0(N1M-1,K)   &
                  +CFX1(I, 0)*A0(N1M  ,K)   &
                  +CFX1(I, 1)*A0(1    ,K)
        ENDDO
      ENDIF

      A0=TMP

      RETURN
      END

!=======================================================================
      SUBROUTINE TEST_FILTERZ(A0)
!=======================================================================
!
!     Test filter in z- for dynamic procedures
!
!     IPX,IPZ: 0 (not periodic), 1 (periodc)
!              These options are defined in the SIDEWALLBC subroutine.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : CFZP1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    ::  A0(M1M,M3M),TMP(M1M,M3M)

      IF (IPZ.EQ.0) THEN
        DO I=1,N1M
          DO K=2,N3M-1
          TMP(I,K  )=CFZP1(K  ,-1)*A0(I,K-1)   &
                    +CFZP1(K  , 0)*A0(I,K  )   &
                    +CFZP1(K  , 1)*A0(I,K+1)
          ENDDO
          TMP(I,1  )=TMP(I,1  )
          TMP(I,N3M)=TMP(I,N3M)
        ENDDO
      ELSE IF (IPZ.EQ.1) THEN
        DO I=1,N1M
          DO K=2,N3M-1
          TMP(I,K  )=CFZP1(K  ,-1)*A0(I,K-1)  &
                    +CFZP1(K  , 0)*A0(I,K  )  &
                    +CFZP1(K  , 1)*A0(I,K+1)
          ENDDO
          TMP(I,1  )=CFZP1(1  ,-1)*A0(I,N3M)   &
                    +CFZP1(1  , 0)*A0(I,1  )   &
                    +CFZP1(1  , 1)*A0(I,2  )
          TMP(I,N3M)=CFZP1(N3M,-1)*A0(I,N3M-1) &
                    +CFZP1(N3M, 0)*A0(I,N3M  ) &
                    +CFZP1(N3M, 1)*A0(I,1    )
        ENDDO
      ENDIF

      A0=TMP

      RETURN
      END

!=======================================================================
      SUBROUTINE NUTZERO
!=======================================================================
!
!     Set TNU to zeo in the solid body,
!     Calculate avg. and max. TNU, and
!     Trace time history of sub-grid scale eddy viscosity (ftrles.dat).
!
!     NNZ: Number of TNU zero point
!     INZ, JNZ, KNZ: Index of TNU zero point
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : TNU,INZ,JNZ,KNZ
      IMPLICIT NONE
      INTEGER*8  :: I,J,K
      INTEGER*8  :: IMAX,JMAX,KMAX,L
      REAL*8     :: TNUAVG,TNUMAX

      IMAX = 0
      JMAX = 0
      KMAX = 0

      IF (IBMON.NE.0) THEN
!$OMP PARALLEL DO
        DO L=1,NNZ
          TNU(INZ(L),JNZ(L),KNZ(L)) = 0.
        ENDDO
      ENDIF

      TNUAVG = 0.
      TNUMAX = 0.
!$OMP PARALLEL DO &
!$OMP reduction(+:TNUAVG) &
!$OMP reduction(MAX:TNUMAX)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        TNUAVG = TNUAVG+TNU(I,J,K)*SDX(I)*SDY(J)*SDZ(K)
        !TNUMAX = AMAX1(TNUMAX,TNU(I,J,K))
        IF (TNU(I,J,K) .GT. TNUMAX) THEN
          TNUMAX=TNU(I,J,K)
          IMAX=I
          JMAX=J
          KMAX=K
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      TNUAVG=TNUAVG/XL/YL/ZL

      WRITE(5002,101) TIME,CSGS,TNUAVG*RE,TNUMAX*RE
  101 FORMAT(F13.5,3ES15.7)

      WRITE(*,99)  CSGS
      WRITE(*,100) TNUAVG*RE
      WRITE(*,110) TNUMAX*RE,XMP(IMAX),YMP(JMAX),ZMP(KMAX),IMAX,JMAX,KMAX
   99 FORMAT('C_VREMAN=',ES18.10)
  100 FORMAT('NUTAVG=  ',ES18.10)
  110 FORMAT('NUTMAX=  ',ES18.10,',',3F10.4,',',3I5)

      RETURN
      END

!=======================================================================
      SUBROUTINE NUTINTERPOL
!=======================================================================
!
!     INTERPOLATION OF NUT
!     MAKE FIELD OF TNUX,TNUY,TNUZ
!
!     B.C. FOR NUT (FOR GM MODEL)
!     AT THE NO SLIP WALL (J=0 WALL)         : NUT=0
!     INLET WALL                             : D(NUT)/DX=0
!     AT THE CONVECTIVE WALL(OUTLET WALL)    : D(NUT)/DX=0
!     FAR FIELD WALL(SIED WALL)              : D(NUT)/DX=0
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : TNU,TNUX,TNUY,TNUZ,INZ,JNZ,KNZ,NWALL_DVM
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      TNUX=0.
      TNUY=0.
      TNUZ=0.

      !* Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
        DO J=0,N2
        DO I=0,N1
           TNU(I,J,0) =TNU(I,J,N3M)
           TNU(I,J,N3)=TNU(I,J,1)
        ENDDO
        ENDDO
      ENDIF

      !* X PERIODICITY
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO
        DO K=0,N3
        DO J=0,N2
           TNU(0 ,J,K)=TNU(N1M,J,K)
           TNU(N1,J,K)=TNU(1  ,J,K)
        ENDDO
        ENDDO
      ENDIF

!$OMP PARALLEL DO
      DO 10 K=K_BGPZ,N3M
      DO 10 J=2,N2M
      DO 10 I=1,N1M
      TNUX(I,J,K)=VVDY(J)*VVDZ(K)*0.25                                  &
                *(SDY(J-1)*(SDZ(K-1)*TNU(I,J,K)+SDZ(K)*TNU(I,J,K-1))    &
                +SDY(J)*(SDZ(K-1)*TNU(I,J-1,K)+SDZ(K)*TNU(I,J-1,K-1)))
   10 CONTINUE

!$OMP PARALLEL DO
      DO 20 K=K_BGPZ,N3M
      DO 20 J=1,N2M
      DO 20 I=I_BGPX,N1M
      TNUY(I,J,K)=VVDZ(K)*VVDX(I)*0.25                                  &
                *(SDZ(K-1)*(SDX(I-1)*TNU(I,J,K)+SDX(I)*TNU(I-1,J,K))    &
                +SDZ(K)*(SDX(I-1)*TNU(I,J,K-1)+SDX(I)*TNU(I-1,J,K-1)))
   20 CONTINUE

!$OMP PARALLEL DO
      DO 30 K=1,N3M
      DO 30 J=2,N2M
      DO 30 I=I_BGPX,N1M
      TNUZ(I,J,K)=VVDX(I)*VVDY(J)*0.25                                  &
                *(SDX(I-1)*(SDY(J-1)*TNU(I,J,K)+SDY(J)*TNU(I,J-1,K))    &
                +SDX(I)*(SDY(J-1)*TNU(I-1,J,K)+SDY(J)*TNU(I-1,J-1,K)))
   30 CONTINUE

!     DEFINITION AT THE BOUNDARY
!
!      DO 40 J=2,N2M
!      JM=JMV(J)
!      DO 40 K=1,N3M
!      KM=KMV(K)
!      TNUY(1,J,K)=
!     &        VVDZ(K)*0.5*(SDZ(K)*TNU(2,J,KM)+SDZ(KM)*TNU(2,J,K))
!      TNUZ(1,J,K)=
!     &        VVDY(J)*0.5*(SDY(J)*TNU(2,JM,K)+SDY(JM)*TNU(2,J,K))
!      TNUY(N1,J,K)=
!     &        VVDZ(K)*0.5*(SDZ(K)*TNU(N1M,J,KM)+SDZ(KM)*TNU(N1M,J,K))
!      TNUZ(N1,J,K)=
!     &        VVDY(J)*0.5*(SDY(J)*TNU(N1M,JM,K)+SDY(JM)*TNU(N1M,J,K))
!   40 CONTINUE
!
!      DO 50 I=2,N1M
!      IM=IMV(I)
!      DO 50 K=1,N3M
!      KM=KMV(K)
!      TNUX(I,1,K)=
!     &        VVDZ(K)*0.5*(SDZ(K)*TNU(I,2,KM)+SDZ(KM)*TNU(I,2,K))
!      TNUZ(I,1,K)=
!     &        VVDX(I)*0.5*(SDX(I)*TNU(IM,2,K)+SDX(IM)*TNU(I,2,K))
!      TNUX(I,N2,K)=
!     &        VVDZ(K)*0.5*(SDZ(K)*TNU(I,N2M,KM)+SDZ(KM)*TNU(I,N2M,K))
!      TNUZ(I,N2,K)=
!     &        VVDX(I)*0.5*(SDX(I)*TNU(IM,N2M,K)+SDX(IM)*TNU(I,N2M,K))
!   50 CONTINUE
!
!      DO 60 I=2,N1M
!      IM=IMV(I)
!      DO 60 J=2,N2M
!      JM=JMV(J)
!      TNUX(I,J,1)=
!     &        VVDY(J)*0.5*(SDY(J)*TNU(I,JM,2)+SDY(JM)*TNU(I,J,2))
!      TNUY(I,J,1)=
!     &        VVDX(I)*0.5*(SDX(I)*TNU(IM,J,2)+SDX(IM)*TNU(I,J,2))
!      TNUX(I,J,N3)=
!     &        VVDY(J)*0.5*(SDY(J)*TNU(I,JM,N3M)+SDY(JM)*TNU(I,J,N3M))
!      TNUY(I,J,N3)=
!     &        VVDX(I)*0.5*(SDX(I)*TNU(IM,J,N3M)+SDX(IM)*TNU(I,J,N3M))
!   60 CONTINUE

      !* X NO PERIODICITY
      IF (IPX.NE.1) THEN
!$OMP PARALLEL DO
      DO K=K_BGPZ,N3M
      DO J=1,N2M
      !TNUY( 1,J,K)=VVDZ(K)*0.5*(SDZ(K)*TNU(  1,J,K-1)+SDZ(K-1)*TNU(  1,J,K))
      TNUY(N1,J,K)=VVDZ(K)*0.5*(SDZ(K)*TNU(N1M,J,K-1)+SDZ(K-1)*TNU(N1M,J,K))
      ENDDO
      ENDDO
!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=2,N2M
      !TNUZ( 1,J,K)=VVDY(J)*0.5*(SDY(J)*TNU(  1,J-1,K)+SDY(J-1)*TNU(  1,J,K))
      TNUZ(N1,J,K)=VVDY(J)*0.5*(SDY(J)*TNU(N1M,J-1,K)+SDY(J-1)*TNU(N1M,J,K))
      ENDDO
      ENDDO
      ENDIF

      !* Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=1,N2M
      DO I=1,N1M
         TNUX(I,J,N3) =TNUX(I,J,  1)
         TNUY(I,J,N3) =TNUY(I,J,  1)
         TNUZ(I,J,N3) =TNUZ(I,J,  1)
      ENDDO
      ENDDO
      ENDIF

      !* X PERIODICITY
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
         TNUX(N1,J,K)=TNUX(1  ,J,K)
         TNUY(N1,J,K)=TNUY(1  ,J,K)
         TNUZ(N1,J,K)=TNUZ(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE RHSSGS
!=======================================================================
!
!     Calculate D/DXj(NU_t*DUj/DXi)*(1-KroneckerDelta(ij)),
!     discretized terms using a RK3 method
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,TNU,TNUX,TNUY,TNUZ,RHS1
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      REAL*8    :: RHSU2,RHSU3,RHSV1,RHSV3,RHSW1,RHSW2

!$OMP PARALLEL DO &
!$OMP PRIVATE(RHSU2,RHSU3)
      DO 10 K=1,N3M
      DO 10 J=1,N2M
      DO 10 I=I_BGPX,N1M
      RHSU2=SSDY(J)*VVDX(I)                          &
           *(TNUZ(I,J+1,K)*(V(I,J+1,K)-V(I-1,J+1,K)) &
            -TNUZ(I,J,K)*(V(I,J,K)-V(I-1,J,K)))
      RHSU3=SSDZ(K)*VVDX(I)                          &
           *(TNUY(I,J,K+1)*(W(I,J,K+1)-W(I-1,J,K+1)) &
            -TNUY(I,J,K)*(W(I,J,K)-W(I-1,J,K)))
      RHS1(I,J,K,1)=RHSU2+RHSU3
   10 CONTINUE

!$OMP PARALLEL DO &
!$OMP PRIVATE(RHSV1,RHSV3)
      DO 20 K=1,N3M
      DO 20 J=2,N2M
      DO 20 I=1,N1M
      RHSV1=SSDX(I)*VVDY(J)                            &
            *(TNUZ(I+1,J,K)*(U(I+1,J,K)-U(I+1,J-1,K))  &
             -TNUZ(I,J,K)*(U(I,J,K)-U(I,J-1,K)))
      RHSV3=SSDZ(K)*VVDY(J)                            &
            *(TNUX(I,J,K+1)*(W(I,J,K+1)-W(I,J-1,K+1))  &
             -TNUX(I,J,K)*(W(I,J,K)-W(I,J-1,K)))
      RHS1(I,J,K,2)=RHSV1+RHSV3
   20 CONTINUE

!$OMP PARALLEL DO &
!$OMP PRIVATE(RHSW1,RHSW2)
      DO 30 K=K_BGPZ,N3M
      DO 30 J=1,N2M
      DO 30 I=1,N1M
      RHSW1=SSDX(I)*VVDZ(K)                            &
             *(TNUY(I+1,J,K)*(U(I+1,J,K)-U(I+1,J,K-1)) &
             -TNUY(I,J,K)*(U(I,J,K)-U(I,J,K-1)))
      RHSW2=SSDY(J)*VVDZ(K)                            &
             *(TNUX(I,J+1,K)*(V(I,J+1,K)-V(I,J+1,K-1)) &
             -TNUX(I,J,K)*(V(I,J,K)-V(I,J,K-1)))
      RHS1(I,J,K,3)=RHSW1+RHSW2
   30 CONTINUE

      RETURN
      END