!=======================================================================
!
!     LICA code
!     Large eddy simulation with an Immersed boundary method in the
!     CArtesian coordinates
!
!     Copyright(c) 2017 All rights reserved by
!     Haecheon Choi
!     Department of Mechanical & Aerospace Engineering
!     Seoul National University
!     E-mail: choi@snu.ac.kr
!     URL: http://tfc.snu.ac.kr
!
!     Noncommercial use with copyrighted mark
!
!     Source code is mainly contributed by
!       Prof. Haecheon Choi, Seoul National University
!       Prof. Dongju Kim,    Kumoh National Institute of Technology
!       Prof. Sungwon Kang,  Sogang University
!       Prof. Jungwoo Kim,   Seoul National University of Science and Technology
!       Prof. Jungil Lee,    Ajou University
!       Mr.   Jonghwan Park  Ph. D. Candidate in Seoul National University
!
!     Q&A: licacodeforum@gmail.com
!
!=======================================================================
      PROGRAM MAIN
      USE MOD_IBM_PRESSURE_NEUMANN
      USE MOD_STL_SDF, ONLY: LOAD_SDF
      USE MOD_ROTOR_MOTION, ONLY: LOAD_ROTOR, UPDATE_ROTOR, ROTOR_RE, REFERENCE_TIME, OMEGA_NOW

      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE
      INTEGER*8   ::  ICFL,JCFL,KCFL
      INTEGER*8   ::  I,J,K,N
      REAL*8      ::  DTCFL
      REAL*8 :: DEMO_END_TIME
      INTEGER*8 :: DEMO_UNIT
      LOGICAL :: DEMO_HAS_LIMIT, USER_STOP_REQUESTED
      INTEGER :: PROGRESS_UNIT

      CALL TIMESTAMP ( )
      CALL real_time(TOTAL_TIME_B)
      CALL LOAD_SDF('body.sdf')
      CALL READINPUTFILE                     ! Read input file (lica.in)
      CALL LOAD_ROTOR('rotor.in')
      RE=ROTOR_RE
      REI=1D0/RE
      CALL UPDATE_ROTOR(0D0)
      DEMO_END_TIME=HUGE(1D0)
      INQUIRE(FILE='run_limit.in',EXIST=DEMO_HAS_LIMIT)
      IF (DEMO_HAS_LIMIT) THEN
        OPEN(NEWUNIT=DEMO_UNIT,FILE='run_limit.in',STATUS='OLD')
        READ(DEMO_UNIT,*) DEMO_END_TIME
        CLOSE(DEMO_UNIT)
      ENDIF
      CALL SIDEWALLBC                        ! BOUNDARY CONDITION OPTION : lica_Ahmed.f90
      CALL GEOM                              ! GRID, MESH, INDEX ....
      CALL GETCDAREA                         ! GET FRONTAL AREA FOR DRAG COEFFI.
      CALL ALLO_ARRAY                        ! ALLOCATE FLOW ARRAY
      IF(IBMON.EQ.1 .AND. IMOVINGON.EQ.1 .AND. IFUNCorPTS.EQ.2) &
      CALL READCURVE                         ! READ THE LAGRANGIAN POINTS OF GEOMETRY

!================ TRACE, IBM, LES INITIALIZATION =======================
      IF (NTEST.NE.0) CALL TRACEINIT        ! TRACING POSITION INITIALIZE
      IF (IBMON.NE.0) THEN
        IF (IMOVINGON.EQ.1) THEN
          CALL FINDFORCING ! Moving STL geometry is built directly in memory
        ELSE
          CALL IBMINIT
          IF (ILES.EQ.1) CALL NUTZEROREAD
        ENDIF
      ENDIF
      IF (ILES.EQ.1)  CALL SGSFILTERINIT     ! CALCULATE FILTER COEFFICIENT

!==================== READ FIELD OR MAKE FIELD =========================
      IF (IREAD.NE.0) THEN
        CALL PREFLD                          ! READ INITIAL FIELD
        IF(ICH.EQ.1) CALL MEANPGR            ! CALCULATE MEAN PRESSURE GRADIENT FOR CONST. MASS FLOW RATE IN CHANNEL FLOW
      ELSE
        CALL MAKEFLD                         ! MAKE INITIAL FIELD
        IF(ICH.EQ.1) CALL MEANPGR
      ENDIF
      CALL CONVRGE1                          ! CONTINUITY CHECK

!========================= INITIALIZE OTHERS ===========================
      IF (IDTOPT .EQ. 0) DT=DTT
      IF ((IDTOPT .NE.0) .AND. (IREAD .EQ. 0)) DT=DTT
      IF (IRESET .EQ. 1) THEN
        IHIST = 0
        TIME  = 0.
      END IF
        NTIME = 0

      CALL UPDATE_ROTOR(TIME)
      IF (IREAD.NE.0 .AND. IMOVINGON.EQ.1) CALL FINDFORCING
      CALL INITIALIZE_PRESSURE_NEUMANN
      CALL EXPORT_IBM_GRID_SLICE
      CALL RK3COEF
      CALL SGSINIT                           ! TNU = 0
      CALL LHSINIT                           ! COEFFICIENTS OF LHS
      CALL POISINIT                          ! POISSON COEFFICIENTS

!======================= INITIALIZE FIELD AVG ==========================
      IF (IAVG .EQ. 1) CALL FIELD_AVG_INIT
      IF (IAVG .EQ. 2) CALL FIELD_AVG_INITZ

!============================== EXTRA ==================================
      CALL FTRFILES(0)                       ! OPEN FTR FILES
      CFLMAX=CFLFAC
      NV  = 101                              ! INITIAL FILE NUMBER(BINARY)
      NAV = 2001                             ! INITIAL AVG. FILE NUMBER(BINARY)

!================= START TIME DEPENDED CALCULATION =====================
      DO 2000 M=1,NTST

        CALL real_time(TIME_BEGIN)

        NTIME=NTIME+1
        IF (IBMON.NE.0) FCVAVG = 0.
        IF (ICH  .EQ.1) PMIAVG = 0.

        CALL CFL(CFLMAX,ICFL,JCFL,KCFL)         ! CALCULATE CFL NUMBER
        IF (IDTOPT.NE.0 .AND. CFLMAX.NE.0.) THEN
           DTCFL = AMIN1(DT*CFLFAC/CFLMAX, DT*(0.80+0.20*CFLFAC/CFLMAX))
           IF (IDTOPT.EQ.1) DT = DTCFL
           IF (IDTOPT.EQ.2) DT = AMIN1(DTT,DTCFL)
        ENDIF
        IF (CFLMAX.GT.(CFLFAC*1.1)) THEN
           PRINT*,' '
           WRITE(*,310) NTIME,CFLMAX,TIME
        ELSE
           PRINT*,' '
           WRITE(*,320) NTIME,TIME,DT
        ENDIF
  310   FORMAT(I15,'   CFL NUMBER EXCEED GIVEN CFL LIMIT :',ES18.5,' AT ',F12.5)
  320   FORMAT('--------------------------',I6,'  TIME=',F10.5,'  DT=',F12.8)

        DT=MIN(DT,DEMO_END_TIME-TIME)
        IF (DT.LE.1D-12) EXIT
        SUBDT = 0D0
!%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN SOLVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        DO 1000 MSUB= 1,3            ! SUB-ITERATION(K=1,2,3)

            ALPHA   = 0.5*(GAMMA(MSUB)+RO(MSUB))
            DTCONST = DT*(GAMMA(MSUB)+RO(MSUB))
            DTCONSTI= 1./DTCONST
            TEST1   = RESID1/DTCONST ! Poisson convergence criterion
            ACOEF   = ALPHA*DT/RE
            ACOEFI  = 1./ACOEF
            PMIAVG  = PMIAVG + 2.*ALPHA*PMI
            SUBDT   = SUBDT+DTCONST  ! Sub time for moving IBM

            IF(ICH.EQ.1) CALL QVOLCALC(QVOL1)
            IF (IMOVINGON.EQ.1) THEN
              CALL UPDATE_ROTOR(TIME+SUBDT)
              CALL FINDFORCING
            ENDIF
            CALL real_time(SGSTIME_B(MSUB))
            IF(ILES .EQ. 1) CALL SGSCALC                   ! CALCULATE COEFFICIENT FOR THE SGS MODLE (VREMAN COEF.)
            CALL real_time(SGSTIME_E(MSUB))
            IF(MSUB .EQ. 1) THEN
              WRITE(*,149) CFLMAX,XMP(ICFL),YMP(JCFL),ZMP(KCFL),ICFL,JCFL,KCFL
 149        FORMAT('CFLMAX=  ',F18.11,',',3F10.4,',',3I5)
            ENDIF
! Moving forcing points already updated before the SGS calculation
            CALL AMBIENT_GHOSTS
            CALL PREPARE_PRESSURE_NEUMANN
            CALL RHSNLHS                                   ! CALCULATE RHS & LHS
            CALL DIVGS(RHS1(1,1,1,2))                      ! CALCULATE RHS OF THE POISSON EQUATION. PHI: RHS1(1,1,1,1), DIVGSUM: RHS1(1,1,1,2)
            CALL real_time(POISSTIME_B(MSUB))
            CALL POISSON(RHS1(1,1,1,1),RHS1(1,1,1,2))      ! SOLVE THE POISSON EQUATION.
            CALL real_time(POISSTIME_E(MSUB))
            IF(ICH.EQ.1) CALL QVOLCALC(QVOL2)
            IF(ICH.EQ.1) PHCAP=(QVOL2-QVOL1)/(XL*YL*ZL)
            CALL UCALC(RHS1(1,1,1,1))                      ! UPDATE VELOCITY
            IF(ICH.EQ.1) CALL QVOLCALC(QVOL)
            CALL PCALC(RHS1(1,1,1,1),RHS1(1,1,1,2))        ! UPDATE PRESSURE
            CALL AMBIENT_GHOSTS
            IF(ICH.EQ.1) CALL MEANPGR
            CALL PRESSURE_STAGE_AUDIT
            IF(IBMON.NE.0) CALL LAGFORCE

 1000   CONTINUE
!%%%%%%%%%%%%%%%%%%%%%%%%%  END MAIN SOLVER  %%%%%%%%%%%%%%%%%%%%%%%%%%%

        TIME=TIME+DT

        CALL CONVRGE1       ! CONVERGENCY CHECK for FLOW
        IF (DVMAX.GT.RESID1) THEN
           WRITE(*,300) NTIME,DVMAX,TIME
        ENDIF
 300    FORMAT(I15,'   LOCAL RESIDUE FOR MASS FLOW CONSV.:',ES18.5,' AT ',F15.6)

        IF ((IBMON.NE.0) .AND. (MASSON.EQ.1)) CALL MASSCHECK ! CHECK MAGNITUDE OF MASS SOURCE/SINK
        IF (MOD(NTIME,NPIN).EQ.0) THEN
          IHIST=IHIST+1
          CALL WRITEHISTORY ! WRITE TIME HISTORY OF SIMULATION
        ENDIF
        IF ((NTEST.NE.0) .AND. (MOD(M,NTPRINT).EQ.0)) CALL TRACER
        IF (IBMON.NE.0) CALL DRAGLIFT    ! DRAG & LIFT OF A BODY
        IF (MOD(NTIME,NPRINT).EQ.0.OR.NTIME.EQ.5) THEN
          CALL WRITEFIELD                ! BINARY FIELD FILE
        ENDIF
        IF (IAVG .EQ. 1) CALL FIELD_AVG  ! AVERAGING
        IF (IAVG .EQ. 2) CALL FIELD_AVGZ
        CALL real_time(TIME_END)
        CALL WRITEFTRTIME
        OPEN(NEWUNIT=PROGRESS_UNIT,FILE='run_progress.dat',STATUS='REPLACE')
        WRITE(PROGRESS_UNIT,'(A)') 'step time_s signed_rpm dt_s CFL divergence_monitor'
        WRITE(PROGRESS_UNIT,'(I12,5ES24.15)') NTIME,TIME*REFERENCE_TIME, &
         OMEGA_NOW/REFERENCE_TIME*60D0/(2D0*ACOS(-1D0)),DT*REFERENCE_TIME,CFLMAX,DVMAX
        CLOSE(PROGRESS_UNIT)
        INQUIRE(FILE='STOP_AFTER_STEP',EXIST=USER_STOP_REQUESTED)
        IF(USER_STOP_REQUESTED)THEN
         WRITE(*,*) 'USER REQUEST: stopping after completed step; final field will be saved'
         EXIT
        ENDIF
        IF (TIME.GE.DEMO_END_TIME-1D-10) EXIT

 2000 CONTINUE
!=================== END TIME DEPENDED CALCULATION =====================
      CALL FTRFILES(1)                   ! CLOSE FTR FILES

      IF (MOD(NTIME,NPRINT).NE.0) THEN
        CALL WRITEFIELD                  ! BINARY FIELD FILE
      ENDIF

      CALL real_time(TOTAL_TIME_E)

      WRITE(*,*) ' '
      WRITE(*,210) TOTAL_TIME_E-TOTAL_TIME_B
 210  FORMAT('TOTAL TIME      : ',F12.2,' SECONDS')

      CALL TIMESTAMP ( )

      CALL DEALLO_ARRAY                  ! DEALLOCATE FLOW ARRAY

      STOP
      END
!=======================================================================
!=======================================================================
!=======================================================================
!=======================================================================
       SUBROUTINE READINPUTFILE
!=======================================================================
       USE MOD_COMMON
       IMPLICIT NONE
       CHARACTER(10) :: DUMMY

       OPEN(10,FILE='lica.in')
       READ(10,*) DUMMY
       READ(10,*) RE,IRESET,IREAD,IAVG,IPZERO,EPS_PTR
       READ(10,*) DUMMY
       READ(10,*) NTST,NPRINT,NPRIAVG,NPIN
       READ(10,*) DUMMY
       READ(10,*) IDTOPT,DTT,CFLFAC
       READ(10,*) DUMMY
       READ(10,*) RESID1,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR
       READ(10,*) DUMMY
       READ(10,*) ILES,ISGSMOD,IDYNMON,CSGS,CSGSPS,IFILTER
       READ(10,*) DUMMY
       READ(10,*) IBMON,MASSON,IINTP,IMOVINGON,IFUNCorPTS
       READ(10,*) DUMMY
       READ(10,'(a)') gridfile
       READ(10,'(a)') fileprevel
       READ(10,*) DUMMY
       READ(10,*) NTEST,NTPRINT

       WRITE(*,101) RE,IRESET,IREAD,IAVG,IPZERO,EPS_PTR
       WRITE(*,102) NTST,NPRINT,NPRIAVG,NPIN
       WRITE(*,103) IDTOPT,DTT,CFLFAC
       WRITE(*,104) RESID1,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR
       WRITE(*,105) ILES,ISGSMOD,IDYNMON,CSGS,CSGSPS,IFILTER
       WRITE(*,106) IBMON,MASSON,IINTP,IMOVINGON,IFUNCorPTS
       WRITE(*,107) gridfile
       WRITE(*,108) fileprevel


 101  FORMAT('RE=',ES11.3,'  IRESET=',I5,'  IREAD=',I5,'  IAVG=',I5,'  IPZERO=',I5,'  EPS_PRT=',F7.3)
 102  FORMAT('NTST=',I10,'  NPRINT=',I8,'  NPRIAVG=',I8,'  NPIN=',I5)
 103  FORMAT('IDTOPT=',I5,'  DT=',ES13.5,'  CFLFAC=',F11.3)
 104  FORMAT('RESID=',ES12.4,'  NLEV=',I5,'  IWC=',I5,'  NBLI=',I5,'  MGITR=',I5,'  IOLDV=',I5,'  IMGSOR=',I5,'  WWSOR=',F7.3)
 105  FORMAT('ILES=',I5,'  ISGSMOD=',I5,'  IDYNMON=',I5,'  CSGS=',F12.4,'  CSGSPS=',F12.4,'  IFILTER=',I5)
 106  FORMAT('IBMON=',I5,'  MASSON=',I5,'  IINTP=',I5,'  IMOVINGON=',I3,'  IFUNCorPTS=',I2)
 107  FORMAT('gridfile  = ', A20)
 108  FORMAT('fileprevel= ', A20)

       REI  = 1./RE

       IF (ISGSMOD.EQ.1) THEN
        WRITE(*,*) 'SMAGORINSKY MODEL IS NOT INCLUDED YET! STOP~~'
        STOP
       ENDIF

       RETURN
       END
!=======================================================================
      SUBROUTINE GEOM
!=======================================================================
!
!     Read 'grid.dat' and define various variables related to grid
!
!     GEOM
!       K_BGPZ,I_BGPX: intial index in x-,z-directions
!     READGRID
!       N1,N2,N3   : number of grid in x-,y-,z-directions
!       N1M,N2M,N3M: values obtained by subtracting 1 from N1,N2,N3
!       XL,YL,ZL   : domain size
!       X,Y,Z      : cell face positions in each direction
!     INDICES
!       IPV,JPV,KPV: I,J,K => X,Y,Z and P => +1, IPV(I)=I+1
!       IMV,JMV,KMV:                    M => -1, IMV(I)=I-1
!     INDXFIX
!       FIXIL,FIXJL,FIXKL: lower boundary indicators, which are only 1
!                          at the lower (I,J,K=1) boundary
!       FIXIU,FIXJU,FIXKU: upper boundary indicators, which are only 1
!                          at the upper (I,J,K=Max) boundary
!     MESHES
!       SDX,SDY,SDZ: distance between cell face,   SDX(I)=X(I+1)-X(I)
!       VDX,VDY,VDZ: distance between cell center, VDX(I)=(SDX(I)+SDX(I-1))/2
!       SSDX,SSDY,SSDZ: reciprocal of SDX,SDY,SDZ for spatial discretization
!       VVDX,VVDY,VVDZ: reciprocal of VDX,VDY,VDZ for spatial discretization
!     PHYSPOS
!       XMP,YMP,ZMP: cell center positions,        XMP(I)=X(I)+SDX(I)/2
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      CALL READGRID !* Read grid.dat
      CALL INDICES  !* Define index variables
      CALL INDXFIX  !* Define index variables for boundary conditions
      CALL MESHES   !* Define mesh sizes
      CALL PHYSPOS  !* Define cell centers

      !* Z periodicity
      IF (IPZ .EQ. 1) THEN
        FIXKL(1)  = 0.
        FIXKU(N3M)= 0.
        SDZ(0)    = SDZ(N3M)
        SDZ(N3)   = SDZ(1)
        SSDZ(0)   = SSDZ(N3M)
        SSDZ(N3)  = SSDZ(1)
        VDZ(1)    = 0.5*(SDZ(1)+SDZ(N3M))
        VDZ(N3)   = 0.5*(SDZ(1)+SDZ(N3M))
        DO K=1,N3
         VVDZ(K)= 1./VDZ(K)
        ENDDO
        KMV(1)  = N3M
        KPV(N3M)= 1
      ENDIF
      IF (IPZ.EQ.0) THEN
        K_BGPZ= 2
      ELSE IF (IPZ.EQ.1) THEN
        K_BGPZ= 1
      ENDIF

      !* X periodicity
      IF (IPX .EQ. 1) THEN
        FIXIL(1)  = 0.
        FIXIU(N1M)= 0.
        SDX(0)    = SDX(N1M)
        SDX(N1)   = SDX(1)
        SSDX(0)   = SSDX(N1M)
        SSDX(N1)  = SSDX(1)
        VDX(1)    = 0.5*(SDX(1)+SDX(N1M))
        VDX(N1)   = 0.5*(SDX(1)+SDX(N1M))
        DO I=1,N1
        VVDX(I)= 1./VDX(I)
        ENDDO
        IMV(1)  = N1M
        IPV(N1M)= 1
      ENDIF
      IF (IPX.EQ.0) THEN
        I_BGPX= 2
      ELSE IF (IPX.EQ.1) THEN
        I_BGPX= 1
      ENDIF

      RETURN
      END

!=======================================================================
     SUBROUTINE READGRID
!=======================================================================
!
!     N1,N2,N3   : number of grid in x-,y-,z-directions
!     N1M,N2M,N3M: values obtained by subtracting 1 from N1,N2,N3
!     XL,YL,ZL   : domain size
!     X,Y,Z      : cell face positions in each direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      OPEN(11,FILE=gridfile)
      READ(11,*) N1,N2,N3
      WRITE(*,101) N1,N2,N3
      READ(11,*) XL,YL,ZL
      WRITE(*,102) XL,YL,ZL
 101  FORMAT('NX=',I12,'  NY=',I12,'  NZ=',I12)
 102  FORMAT('XL=',F12.4,'  YL=',F12.4,'  ZL=',F12.4)

      N1M=N1-1
      N2M=N2-1
      N3M=N3-1

      IF ((N1.GT.M1) .OR. (N2.GT.M2) .OR. (N3.GT.M3)) THEN
         PRINT*, 'ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF

      READ(11,*) (X(I),I=1,N1)
      READ(11,*) (Y(J),J=1,N2)
      READ(11,*) (Z(K),K=1,N3)
      CLOSE(11)

      RETURN
      END

!=======================================================================
      SUBROUTINE INDICES
!=======================================================================
!
!     IPV,JPV,KPV: I,J,K => X,Y,Z and P => +1, IPV(I)=I+1
!     IMV,JMV,KMV:                    M => -1, IMV(I)=I-1
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

!---- STREAMWISE DIRECTION
      DO I=1,N1M
      IPV(I)=I+1
      IMV(I)=I-1
      ENDDO

!---- NORMAL DIRECTION
      DO J=1,N2M
      JPV(J)=J+1
      JMV(J)=J-1
      ENDDO

!---- SPANWISE DIRECTION
      DO K=1,N3M
      KPV(K)=K+1
      KMV(K)=K-1
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE INDXFIX
!=======================================================================
!
!     FIXIL,FIXJL,FIXKL: Lower boundary indicators, which are only 1
!                        at the lower (I,J,K=1) boundary
!     FIXIU,FIXJU,FIXKU: Upper boundary indicators, which are only 1
!                        at the upper (I,J,K=Max) boundary
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      DO I=1,N1M
      FIXIL(I)  = 0.
      FIXIU(I)  = 0.
      ENDDO
      FIXIL(1)  = 1.
      FIXIU(N1M)= 1.

      DO J=1,N2M
      FIXJL(J)  = 0.
      FIXJU(J)  = 0.
      ENDDO
      FIXJL(1)  = 1.
      FIXJU(N2M)= 1.

      DO K=1,N3M
      FIXKL(K)  = 0.
      FIXKU(K)  = 0.
      ENDDO
      FIXKL(1)  = 1.
      FIXKU(N3M)= 1.

      RETURN
      END

!=======================================================================
      SUBROUTINE MESHES
!=======================================================================
!
!     SDX,SDY,SDZ: distance between cell face,   SDX(I)=X(I+1)-X(I)
!     VDX,VDY,VDZ: distance between cell center, VDX(I)=(SDX(I)+SDX(I-1))/2
!     SSDX,SSDY,SSDZ: reciprocal of SDX,SDY,SDZ for spatial discretization
!     VVDX,VVDY,VVDZ: reciprocal of VDX,VDY,VDZ for spatial discretization
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      DO I=1,N1M
      SDX(I)=X(I+1)-X(I)
      ENDDO

      DO I=2,N1M
      VDX(I) = 0.5*(SDX(I)+SDX(I-1))
      ENDDO
      VDX(1) = 0.5*SDX(1)             ! LEFT BOUNDED COND.
      VDX(N1)= 0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

      DO I=1,N1M
      SSDX(I) = 1./SDX(I)
      ENDDO
      DO I=1,N1
      VVDX(I) = 1./VDX(I)
      ENDDO

      DO J=1,N2M
      SDY(J) = Y(J+1)-Y(J)
      ENDDO

      DO J=2,N2M
      VDY(J) = 0.5*(SDY(J)+SDY(J-1))
      ENDDO
      VDY(1) = 0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
      VDY(N2)= 0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

      DO J=1,N2M
      SSDY(J)= 1./SDY(J)
      ENDDO
      DO J=1,N2
      VVDY(J)= 1./VDY(J)
      ENDDO

      IF (N3M.EQ.1) THEN
        ZL = 1.
        SDZ(1) = ZL/N3M
        SSDZ(1)= 1./SDZ(1)
        VDZ(1) = SDZ(1)
        VVDZ(1)= SSDZ(1)

      ELSE
        DO K=1,N3M
        SDZ(K)=Z(K+1)-Z(K)
        ENDDO
        DO K=2,N3M
        VDZ(K)=0.5*(SDZ(K)+SDZ(K-1))
        ENDDO
        VDZ(1)=0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
        VDZ(N3)=0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.
        DO K=1,N3M
        SSDZ(K)=1./SDZ(K)
        ENDDO
        DO K=1,N3
        VVDZ(K)=1./VDZ(K)
        ENDDO
      ENDIF

!-----set by zero
      SDX(0) = 0.
      SDX(M1)= 0.
      SDY(0) = 0.
      SDY(M2)= 0.
      SDZ(0) = 0.
      SDZ(M3)= 0.

      RETURN
      END

!=======================================================================
      SUBROUTINE PHYSPOS
!=======================================================================
!
!     XMP,YMP,ZMP: cell center positions,        XMP(I)=X(I)+SDX(I)/2
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      DO I=1,N1M
      XMP(I)=X(I)+0.5*SDX(I)
      ENDDO
      XMP(N1)=X(N1)
      XMP(0)=X(1)

      DO J=1,N2M
      YMP(J)=Y(J)+0.5*SDY(J)
      ENDDO
      YMP(N2)=Y(N2)
      YMP(0)=Y(1)

      DO K=1,N3M
      ZMP(K)=Z(K)+0.5*SDZ(K)
      ENDDO
      ZMP(N3)=Z(N3)
      ZMP(0)=Z(1)

      RETURN
      END

!=======================================================================
      SUBROUTINE RK3COEF
!=======================================================================
!
!     Coefficients used in 3rd-Order Runge-Kutta method (subroutine RHS)
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE

      GAMMA(1) = 8./15.
      GAMMA(2) = 5./12.
      GAMMA(3) = 3./4.
      RO(1) = 0.
      RO(2) = -17./60.
      RO(3) = -5./12.

      RETURN
      END

!=======================================================================
      SUBROUTINE LHSINIT
!=======================================================================
!
!     CALCULATE COEFFICIENTS FOR LHS
!       These coefficients are used in LHSU,LHSV,LHSW subroutines to
!       calculate intermediate velocity.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: IC,JC,KC
      INTEGER*8 :: IM,IP,JM,JP,KM,KP

      DO 10 IC=I_BGPX,N1M
      IM=IMV(IC)
      AIU(IC)=-VVDX(IC)*SSDX(IM)                 !   :i-1
      CIU(IC)=-VVDX(IC)*SSDX(IC)                 !   :i+1
   10 CONTINUE

      DO 20 IC=1,N1M
      IP=IPV(IC)
      AIVW(IC)=-VVDX(IC)*SSDX(IC)
      CIVW(IC)=-VVDX(IP)*SSDX(IC)
   20 CONTINUE

      DO 30 JC=2,N2M
      JM=JMV(JC)
      AJV(JC)=-VVDY(JC)*SSDY(JM)
      CJV(JC)=-VVDY(JC)*SSDY(JC)
   30 CONTINUE

      DO 40 JC=1,N2M
      JP=JPV(JC)
      AJUW(JC)=-VVDY(JC)*SSDY(JC)
      CJUW(JC)=-VVDY(JP)*SSDY(JC)
   40 CONTINUE

      DO 50 KC=K_BGPZ,N3M
      KM=KMV(KC)
      AKW(KC)=-VVDZ(KC)*SSDZ(KM)
      CKW(KC)=-VVDZ(KC)*SSDZ(KC)
   50 CONTINUE

      DO 60 KC=1,N3M
      KP=KPV(KC)
      AKUV(KC)=-VVDZ(KC)*SSDZ(KC)
      CKUV(KC)=-VVDZ(KP)*SSDZ(KC)
   60 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRACEINIT
!=======================================================================
!
!     Initializing the time-tracing position.
!
!     NTEST: number of tracing positions in lica.in
!     ITR(N), JTR(N), KTR(N): tracing indices listed under the NTEST in lica.in
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: N

      DO N=1,NTEST
        READ(10,*) ITR(N),JTR(N),KTR(N)
      ENDDO

      WRITE(*,*) '========= TRACE POSITION ========='
      DO  N=1,NTEST
        WRITE(*,16) N,ITR(N),JTR(N),KTR(N),XMP(ITR(N)),YMP(JTR(N)),ZMP(KTR(N))
      ENDDO

   16 FORMAT(I5,3I6,3F13.5)

      RETURN
      END

!=======================================================================
      SUBROUTINE TRACER
!=======================================================================
!
!     This subroutine is called when MOD(M,NTPRINT)=0.
!     Velocity components and the pressure at the cell center of the
!     tracing cells are saved.
!
!     NTPRINT : determining the time-tracing interval
!
!      5003->FILE='ftru.dat'
!      5004->FILE='ftrv.dat'
!      5005->FILE='ftrw.dat'
!      5006->FILE='ftrp.dat'
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P
      IMPLICIT NONE
      INTEGER*8 :: N

      WRITE(5003,101)TIME,(0.5*(U(ITR(N),JTR(N),KTR(N))              &
                             +U(ITR(N)+1,JTR(N),KTR(N))),N=1,NTEST)
      WRITE(5004,101)TIME,(0.5*(V(ITR(N),JTR(N),KTR(N))              &
                             +V(ITR(N),JTR(N)+1,KTR(N))),N=1,NTEST)
      WRITE(5005,101)TIME,(0.5*(W(ITR(N),JTR(N),KTR(N))              &
                             +W(ITR(N),JTR(N),KPV(KTR(N)))),N=1,NTEST)
      WRITE(5006,101)TIME,(P(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
  101 FORMAT(F15.7,10000ES15.7)

      RETURN
      END

!=======================================================================
      SUBROUTINE IBMINIT
!=======================================================================
!
!     Read IBM pre-processing data which is calculated by using
!       ibmpre_cyl.f90
!
!     Option
!       IINTP: IBM interpolation ON(1)/OFF(0)
!
!     Variables
!       NFC_INTP():  Number of momentum forcing points (parameter: u,v,w)
!                    which conduct interpolation
!       NFC_INNER(): Number of momentum forcing points (parameter: u,v,w)
!                    which do not conduct interpolation
!       IFC(),JFC(),KFC(): I,J,K index of momentum forcing points
!       MPI(): Index increase or decrease in order to know neighboring
!              outside points for interpolation
!              (named after Minus Plus for Interpolation)
!              (2nd parameter: u v w, 3rd parameter: x y z)
!              (+1: increase, -1: decrease)
!              (0: NO interpolation for that direction)
!       GFI(): Geometric factor for interpolation
!       FCV(): Momentum forcing to satisfy no-slip boundary codition
!              on the body surface
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : IFC,JFC,KFC,MPI,GFI,FCV
      IMPLICIT NONE
      INTEGER*8  :: N,L
      INTEGER*8  :: I,J,K,ITMP

      !* Stationary body (if IMOVINGON=1, goto FINDFORCING subroutines.)
      IF(IMOVINGON.EQ.0) THEN
        OPEN(11,FILE='ibmpre_intp.bin',FORM='UNFORMATTED')
        READ(11) NFC_INTP(1),NFC_INTP(2),NFC_INTP(3)
        READ(11) NFC_INNER(1),NFC_INNER(2),NFC_INNER(3)
        READ(11) (IFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
        READ(11) (IFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
        READ(11) (IFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
        READ(11) (JFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
        READ(11) (JFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
        READ(11) (JFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
        READ(11) (KFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
        READ(11) (KFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
        READ(11) (KFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
        READ(11) ((MPI(N,1,L),N=1,NFC_INTP(1)),L=1,3)
        READ(11) ((MPI(N,2,L),N=1,NFC_INTP(2)),L=1,3)
        READ(11) ((MPI(N,3,L),N=1,NFC_INTP(3)),L=1,3)
        READ(11) ((((GFI(N,1,I,J,K),N=1,NFC_INTP(1)),I=0,2),J=0,2),K=0,2)
        READ(11) ((((GFI(N,2,I,J,K),N=1,NFC_INTP(2)),I=0,2),J=0,2),K=0,2)
        READ(11) ((((GFI(N,3,I,J,K),N=1,NFC_INTP(3)),I=0,2),J=0,2),K=0,2)
        CLOSE(11)

        WRITE(*,35) NFC_INTP(1),NFC_INTP(2),NFC_INTP(3)
        WRITE(*,36) NFC_INNER(1),NFC_INNER(2),NFC_INNER(3)
 35     FORMAT('FU_intp :',I12,'  FV_intp :',I12,'  FW_intp :',I12)
 36     FORMAT('FU_inner:',I12,'  FV_inner:',I12,'  FW_inner:',I12)

        IF (IINTP.EQ.0) THEN
          GFI = 0.
          MPI = 0
          DO L=1,3
            ITMP = NFC_INTP(L)+NFC_INNER(L)
            NFC_INTP(L) = 0
            NFC_INNER(L)= ITMP
          ENDDO

          WRITE(*,* ) 'IINTP = 0, RESET OF THE NUMBER OF MTM FC POINTS'
          WRITE(*,35) NFC_INTP(1),NFC_INTP(2),NFC_INTP(3)
          WRITE(*,36) NFC_INNER(1),NFC_INNER(2),NFC_INNER(3)
        ENDIF
      ENDIF

      !* Initialization of FCV to zero
      FCV = 0.

      RETURN
      END

!=======================================================================
      SUBROUTINE RHS_IBM
!=======================================================================
!
!     Calculate momentum forcing
!
!     Option
!       IMOVINGON = 0, Stationary body => UBODY,VBODY,WBODY for translational vel.
!       IMOVINGON = 1, Moving body     => UBD,VBD,WBD in lica_cylinder.f90
!
!     Variables
!       UTARG,VTARG,WTARG: Target velocity to satisfy no-slip b.c.
!       FCV   : Momentum forcing from the target velocity
!       FCVAVG: Averaging forcing values    to calculate the force on a body
!       DUDTR : Time derivative of velocity to calculate the force on a body
!               Ref. Lee et al., 2011, Sources of spurious force
!               oscillations from an immersed boundary method for moving
!               -body problems, J. Comp. Phys., 230, 2677-2695.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : RHS1,U,V,W,IFC,JFC,KFC,MPI,GFI,FCV,FCVAVG,DUDTR
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: SURFACE_VELOCITY,EXTENSION_VELOCITY
      USE MOD_IBM_COUPLED_TARGETS, ONLY: COMPUTE_COUPLED_TARGETS,COUPLED_TARGET
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,N,L,II,JJ,KK,IP,JP,KP,IPP,JPP,KPP
      REAL*8    :: UBODY,VBODY,WBODY,DTI
      REAL*8    :: UTARG,VTARG,WTARG
      REAL*8    :: RHSTMP(M1M,M2M,M3M,3)
      REAL*8    :: UBD,VBD,WBD

!---- Compute target velocities & forcing values at forcing points
      UBODY = 0.
      VBODY = 0.
      WBODY = 0.
      DTI   = 1./DT

!$OMP PARALLEL DO
      DO K= 1,N3M
      DO J= 1,N2M
      DO I= 1,N1M
      DO L= 1,3
        RHSTMP(I,J,K,L)=RHS1(I,J,K,L)
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      CALL COMPUTE_COUPLED_TARGETS(RHSTMP)

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,UTARG) &
!$OMP firstprivate(UBODY)
      DO N=1,NFC_INTP(1)
         II =IFC(N,1)
         JJ =JFC(N,1)
         KK =KFC(N,1)
         IP =IFC(N,1)+MPI(N,1,1)
         JP =JFC(N,1)+MPI(N,1,2)
         KP =KFC(N,1)+MPI(N,1,3)
         IPP=IFC(N,1)+MPI(N,1,1)*2
         JPP=JFC(N,1)+MPI(N,1,2)*2
         KPP=KFC(N,1)+MPI(N,1,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         IF (IMOVINGON.EQ.1) UBODY=SURFACE_VELOCITY(N,1)

         UTARG=COUPLED_TARGET(N,1)
         FCV(N,1)=(UTARG-U(II,JJ,KK)-RHSTMP(II,JJ,KK,1))*DTI
         DUDTR(N,1)=(UTARG-U(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,1)=UTARG-U(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,VTARG) &
!$OMP firstprivate(VBODY)
      DO N=1,NFC_INTP(2)
         II =IFC(N,2)
         JJ =JFC(N,2)
         KK =KFC(N,2)
         IP =IFC(N,2)+MPI(N,2,1)
         JP =JFC(N,2)+MPI(N,2,2)
         KP =KFC(N,2)+MPI(N,2,3)
         IPP=IFC(N,2)+MPI(N,2,1)*2
         JPP=JFC(N,2)+MPI(N,2,2)*2
         KPP=KFC(N,2)+MPI(N,2,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         IF (IMOVINGON.EQ.1) VBODY=SURFACE_VELOCITY(N,2)

         VTARG=COUPLED_TARGET(N,2)
         FCV(N,2)=(VTARG-V(II,JJ,KK)-RHSTMP(II,JJ,KK,2))*DTI
         DUDTR(N,2)=(VTARG-V(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,2)=VTARG-V(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,IP,JP,KP,IPP,JPP,KPP,WTARG) &
!$OMP firstprivate(WBODY)
      DO N=1,NFC_INTP(3)
         II =IFC(N,3)
         JJ =JFC(N,3)
         KK =KFC(N,3)
         IP =IFC(N,3)+MPI(N,3,1)
         JP =JFC(N,3)+MPI(N,3,2)
         KP =KFC(N,3)+MPI(N,3,3)
         IPP=IFC(N,3)+MPI(N,3,1)*2
         JPP=JFC(N,3)+MPI(N,3,2)*2
         KPP=KFC(N,3)+MPI(N,3,3)*2

         IF ((IPP.GE.N1).OR.(IPP.LE.0)) CALL REINDEX_I(IP,IPP)
         IF ((KPP.GE.N3).OR.(KPP.LE.0)) CALL REINDEX_K(KP,KPP)
         IF (IMOVINGON.EQ.1) WBODY=SURFACE_VELOCITY(N,3)

         WTARG=COUPLED_TARGET(N,3)
         FCV(N,3)=(WTARG-W(II,JJ,KK)-RHSTMP(II,JJ,KK,3))*DTI
         DUDTR(N,3)=(WTARG-W(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,3)=WTARG-W(II,JJ,KK)
      ENDDO

!*************************       CAUTION       *************************
!     FOLLOWING 3 DO-LOOPS REQUIRES THE OPTION OF
!     '-Wf "-pvctl vwork=stack"'
!     WHEN THE CODE IS COMPILED ON NEC MACHINE.
!     FOR DETAILED EXPLANATION OF THE REASON,
!     YOU CAN CONSULT "NEC PORTING GUIDE"
!     PROVIDED BY NEC SUPERCOMPUTIONG CENTER.
!***********************************************************************

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,UTARG) &
!$OMP firstprivate(UBODY)
      DO N=NFC_INTP(1)+1,NFC_INTP(1)+NFC_INNER(1)
         II=IFC(N,1)
         JJ=JFC(N,1)
         KK=KFC(N,1)
         IF (IMOVINGON.EQ.1) UBODY=EXTENSION_VELOCITY(N,1)
         UTARG=UBODY
         FCV(N,1)=(UTARG-U(II,JJ,KK)-RHS1(II,JJ,KK,1))*DTI
         DUDTR(N,1)=(UTARG-U(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,1)=UTARG-U(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,VTARG) &
!$OMP firstprivate(VBODY)
      DO N=NFC_INTP(2)+1,NFC_INTP(2)+NFC_INNER(2)
         II=IFC(N,2)
         JJ=JFC(N,2)
         KK=KFC(N,2)
         IF (IMOVINGON.EQ.1) VBODY=EXTENSION_VELOCITY(N,2)
         VTARG=VBODY
         FCV(N,2)=(VTARG-V(II,JJ,KK)-RHS1(II,JJ,KK,2))*DTI
         DUDTR(N,2)=(VTARG-V(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,2)=VTARG-V(II,JJ,KK)
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(II,JJ,KK,WTARG) &
!$OMP firstprivate(WBODY)
      DO N=NFC_INTP(3)+1,NFC_INTP(3)+NFC_INNER(3)
         II=IFC(N,3)
         JJ=JFC(N,3)
         KK=KFC(N,3)
         IF (IMOVINGON.EQ.1) WBODY=EXTENSION_VELOCITY(N,3)
         WTARG=WBODY
         FCV(N,3)=(WTARG-W(II,JJ,KK)-RHS1(II,JJ,KK,3))*DTI
         DUDTR(N,3)=(WTARG-W(II,JJ,KK))*DTI
         RHS1(II,JJ,KK,3)=WTARG-W(II,JJ,KK)
      ENDDO

!-----compute average forcing values during RK3 steps
!$OMP PARALLEL DO
      DO L=1,3
      DO N=1,NFC_INTP(L)+NFC_INNER(L)
         FCVAVG(N,L)=FCVAVG(N,L)+FCV(N,L)
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE REINDEX_I(IP,IPP)
!=======================================================================
!
!     Assign appropriate indices to the variables used in IBM interpolation
!     in order to satisfy the periodic BC in x-direction.
!      IP =I +1
!      IPP=IP+1
!     *Generally, (IPP.GE.N1)~ situations occur in case of periodic BC
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: IP,IPP

      IF (IPP.EQ.N1) THEN
       IPP = 1
      ELSEIF (IPP.EQ.N1+1) THEN
       IP  = 1
       IPP = 2
      ELSEIF (IPP.EQ.0) THEN
       IPP = N1M
      ELSEIF (IPP.EQ.-1) THEN
       IP  = N1M
       IPP = N1M-1
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE REINDEX_K(KP,KPP)
!=======================================================================
!
!     Same as the SUBROUTINE REINDEX_I(IP,IPP)
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8     :: KP,KPP

      IF (KPP.EQ.N3) THEN
       KPP = 1
      ELSEIF (KPP.EQ.N3+1) THEN
       KP  = 1
       KPP = 2
      ELSEIF (KPP.EQ.0) THEN
       KPP = N3M
      ELSEIF (KPP.EQ.-1) THEN
       KP  = N3M
       KPP = N3M-1
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE MASSCHECK
!=======================================================================
!
!     Calculate maximum mass source/sink in IBM
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : QMASS
      IMPLICIT NONE
      INTEGER*8     :: I,J,K

      QMMAX = 0.

!$OMP PARALLEL DO &
!$OMP reduction(MAX:QMMAX)
      DO 100 K=1,N3M
      DO 100 J=1,N2M
      DO 100 I=1,N1M
      QMMAX = AMAX1(QMMAX,ABS(QMASS(I,J,K)))
 100  CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE CFL(CFLM,ICFL,JCFL,KCFL)
!=======================================================================
!
!     Calculate the maximum CFL number of flow field
!       CFL number is a necessary condition for convergence while solving certain
!       PDE numerically. (http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition)
!
!     Using a semi-implicit scheme, CFL# has to be within the required
!      CFL limit provided by the user in "lica.in". The CFL# limit varies with the
!      numerical scheme and the kind of problem.
!
!     CFL#=U*DT/DX
!     CFLI=(Uc/DX+Vc/DY+Wc/DZ)*DT  (c for cell center)
!     I,J,KCFL: index where the CFL# is maximum
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8   ::  ICFL,JCFL,KCFL
      REAL*8      ::  CFLM,CFLI1,CFLI2,CFLI3,CFLI

      CFLM = 0.
      ICFL=1;JCFL=1;KCFL=1

!$OMP PARALLEL DO &
!$OMP private(CFLI1,CFLI2,CFLI3,CFLI) &
!$OMP reduction(MAX:CFLM)
      DO 13 K=1,N3M
      DO 13 J=1,N2M
      DO 13 I=1,N1M
         CFLI1=ABS(U(I,J,K)+U(I+1,J,K))*SSDX(I)
         CFLI2=ABS(V(I,J,K)+V(I,J+1,K))*SSDY(J)
         CFLI3=ABS(W(I,J,K)+W(I,J,K+1))*SSDZ(K)
         CFLI=(CFLI1+CFLI2+CFLI3)*0.5*DT
!         CFLM=AMAX1(CFLI,CFLM)
         IF (CFLI .GE. CFLM) THEN
              CFLM=CFLI
              ICFL=I
              JCFL=J
              KCFL=K
         ENDIF
   13 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE CONVRGE1
!=======================================================================
!
!     Calculate maximum of DIVERGENCE U in order to check convergence
!     QMASS is subtracted when the IBM mass-source/sink is activated (MASSON=1)
!
!     If maximum divervence is greater than the residue -Poisson convergence
!     criterion (RESID1)- set by the user, a line informing the local DVMAX will be printed.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,QMASS
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8 :: DVG11,DVG12,DVG13,DVG1

      DVMAX = 0.

!$OMP PARALLEL DO                        &
!$OMP private(DVG11,DVG12,DVG13,DVG1)    &
!$OMP reduction(MAX:DVMAX)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
         DVG11=(U(I+1,J,K)-U(I,J,K))*SSDX(I)
         DVG12=(V(I,J+1,K)-V(I,J,K))*SSDY(J)
         DVG13=(W(I,J,K+1)-W(I,J,K))*SSDZ(K)
         IF (MASSON.EQ.1) THEN
          DVG1=DVG11+DVG12+DVG13-QMASS(I,J,K)
         ELSE
          DVG1=DVG11+DVG12+DVG13
         ENDIF
         DVMAX=AMAX1(ABS(DVG1),DVMAX)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
        SUBROUTINE RHSNLHS
!=======================================================================
!
!     Solving N-S equation to get intermediate velocity, u_i hat
!     Intermediate velocity, u_i hat: {u_i}^(k-1) -> u_i hat -> {u_i}^k
!           is based on the pressure at timestep k (p^k).
!           Indeed, u_i hat is not satisfied with continuity equation.
!           For details about u_i hat, see 'Fractional step method'.
!
!     RHS sub-grid scale part (LES) -> RHS -> RHS IB(Immersed boundary)
!     -> Computing convective boundary condition (exit condition)
!     -> Other boundary conditions -> LHS (Get u_i hat)
!
!     Schemes used in this code
!       Time : 3rd order Runge-Kutta (Convection/Non-linear term) +
!              2nd order Crank-Nicolson (Diffusion/Linear term)
!       Space: 2nd order Central-difference
!
!     Non-linear term: N(u_i)=  d(u_i*u_j)/d(x_j)
!     Linear term    : L(u_i)=( d(d(u_i))/(d(x_j)d(x_j)) )/Re
!
!
!     See the following papers for details:
!
!     Papers for Fractional step method (computation procedure):
!      Kim, J., Moin, P. & Moser, R. 1987 Turbulence statistics in
!        fully developed channel flow at low Reynolds number.
!        J. Fluid Mech. 177, 133.
!      Choi, H., Moin, P. & Kim, J. 1994 Active turbulence control for
!        drag reduction in wall-bounded flows. J. Fluid Mech. 262,
!        75-110.
!
!     Paper for Large eddy simulation (Dynamic global model):
!      Lee, J., Choi, H. & Park, N. 2010 Dynamic global model for
!        large eddy simulation of transient flow. Phys. Fluids 22,
!        075106.
!
!     Paper for Immersed Boundary(IB) method:
!      Kim, J., Kim, D. & Choi, H. 2001 An Immersed-Boundary Finite
!        Volume Method for Simulations of Flow in Complex Geometries.
!        J. Comput. Phys. 171 132-150.
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       IMPLICIT NONE

       CALL real_time(RHSTIME_B(MSUB)) ! stamping present real time
       IF (ILES .EQ. 1) CALL RHSSGS    ! SubGrid-Scale computations (LES case)
       CALL RHS
       IF (IBMON.NE.0) CALL RHS_IBM    ! IB method computation
       CALL CONVBC                     ! defined in lica_[bodyname, e.g. cylinder].f90 file
       IF (ICH.NE.1) CALL RHSINCORPBC
       CALL real_time(RHSTIME_E(MSUB))

       CALL real_time(LHSTIME_B(MSUB))
       CALL LHSU
       CALL LHSV
       IF (N3M.NE.1) CALL LHSW
       CALL real_time(LHSTIME_E(MSUB))

       ! WHEN DO BLOWING/SUCTION RETRV_UVW SHOULD BE TURNED ON
       CALL RETRV_UVW

       RETURN
       END

!=======================================================================
       SUBROUTINE RHS
!=======================================================================
!
!     Computing intermediate velocity, u hat, step in delta form
!
!     variables in common:
!           x, y, z         : coordinate direction
!           u, v, w         : velocity for x, y, z direction
!           n, s, e, w, c, f: + & - for each x, y, z direction
!           AN, AL, TAL     : Non-linear, linear, turbulent (SGS)
!
!     RHS1(x,y,z,L):
!           RHS term consists of Non-linear/Linear/SGS terms
!           Components for IB method will be added in RHS_IBM subroutine
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,P,TNU,TNUX,TNUY,TNUZ,RHS1   &
                                ,RK3XO,RK3YO,RK3ZO,RK3XOO,RK3YOO,RK3ZOO
       USE MOD_IBM_PRESSURE_NEUMANN, ONLY: PRESSURE_GRADIENT
      IMPLICIT NONE
       INTEGER*8   :: I,J,K

!------------ variables for U (x, streamwise) velocity
       REAL*8      :: UN,US,UC,UF,VNX,VSX,WCX,WFX
       REAL*8      :: ANXE,ANXW,ANXX,ANX1,ANX2,RK3X
       REAL*8      :: ALX1,ALX2,ALX3,ALX4,ALX5,ALX6,ALXX,ALXY,ALXZ,ALX
       REAL*8      :: TALXX,TALXY,TALXZ
!------------ variables for V (y, traseverse) velocity
       REAL*8      :: VE,VW,VC,VF,UEY,UWY,WCY,WFY
       REAL*8      :: ANYN,ANYS,ANYY,ANY1,ANY2,RK3Y
       REAL*8      :: ALY1,ALY2,ALY3,ALY4,ALY5,ALY6,ALYX,ALYY,ALYZ,ALY
       REAL*8      :: TALYX,TALYY,TALYZ
!------------ variables for W (z, spanwise) velocity
       REAL*8      :: WE,WW,WN,WS,UEZ,UWZ,VNZ,VSZ
       REAL*8      :: ANZC,ANZF,ANZZ,ANZ1,ANZ2,RK3Z
       REAL*8      :: ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZX,ALZY,ALZZ,ALZ
       REAL*8      :: TALZX,TALZY,TALZZ

!-----RHS1 calculation for u momentum -----------------
!$OMP PARALLEL DO  &
!$OMP private(UN,US,UC,UF,VNX,VSX,WCX,WFX)   &
!$OMP private(ANXE,ANXW,ANXX,ANX1,ANX2,RK3X) &
!$OMP private(ALX1,ALX2,ALX3,ALX4,ALX5,ALX6) &
!$OMP private(ALXX,ALXY,ALXZ,ALX,TALXX,TALXY,TALXZ)
       DO K=1,N3M
       DO J=1,N2M
       DO I=I_BGPX,N1M

       UN=0.5*(SDY(J+1)*U(I,J,K)+SDY(J)*U(I,J+1,K))*VVDY(J+1)    &
           *(1.-FIXJU(J))+U(I,N2,K)*FIXJU(J)
       US=0.5*(SDY(J)*U(I,J-1,K)+SDY(J-1)*U(I,J,K))*VVDY(J)      &
           *(1.-FIXJL(J))+U(I,0,K)*FIXJL(J)
       UC=0.5*(SDZ(K+1)*U(I,J,K)+SDZ(K)*U(I,J,K+1))*VVDZ(K+1)    &
           *(1.-FIXKU(K))+U(I,J,N3)*FIXKU(K)
       UF=0.5*(SDZ(K)*U(I,J,K-1)+SDZ(K-1)*U(I,J,K))*VVDZ(K)      &
           *(1.-FIXKL(K))+U(I,J,0)*FIXKL(K)
       VNX=0.5*(SDX(I-1)*V(I,J+1,K)+SDX(I)*V(I-1,J+1,K))*VVDX(I)
       VSX=0.5*(SDX(I-1)*V(I,J,K)+SDX(I)*V(I-1,J,K))*VVDX(I)
       WCX=0.5*(SDX(I-1)*W(I,J,K+1)+SDX(I)*W(I-1,J,K+1))*VVDX(I)
       WFX=0.5*(SDX(I-1)*W(I,J,K)+SDX(I)*W(I-1,J,K))*VVDX(I)

       ANXE=0.5*(U(I+1,J,K)+U(I,J,K))
       ANXW=0.5*(U(I,J,K)+U(I-1,J,K))
       ANXX=(ANXE**2-ANXW**2)*VVDX(I)
       ANX1=(UN*VNX-US*VSX)*SSDY(J)
       ANX2=(UC*WCX-UF*WFX)*SSDZ(K)
       RK3X=-ANXX-ANX1-ANX2               ! Non-linear term at k-substep

       ALX1=(U(I+1,J,K)-U(I,J,K))*SSDX(I)
       ALX2=(U(I,J,K)-U(I-1,J,K))*SSDX(I-1)
       ALX3=(U(I,J+1,K)-U(I,J,K))*VVDY(J+1)*(1.-FIXJU(J)*FLOAT(JUT))
       ALX4=(U(I,J,K)-U(I,J-1,K))*VVDY(J)  *(1.-FIXJL(J)*FLOAT(JUB))
       ALX5=(U(I,J,K+1)-U(I,J,K))*VVDZ(K+1)*(1.-FIXKU(K)*FLOAT(KUT))
       ALX6=(U(I,J,K)-U(I,J,K-1))*VVDZ(K)  *(1.-FIXKL(K)*FLOAT(KUB))
       ALXX=(ALX1-ALX2)*VVDX(I)
       ALXY=(ALX3-ALX4)*SSDY(J)
       ALXZ=(ALX5-ALX6)*SSDZ(K)
       ALX=REI*(ALXX+ALXY+ALXZ)           ! Linear terms at k-substep

!-----LES
      IF (ILES.EQ.1) THEN
       TALXX=(TNU(I,J,K)*ALX1-TNU(I-1,J,K)*ALX2)*VVDX(I)
       TALXY=(TNUZ(I,J+1,K)*ALX3-TNUZ(I,J,K)*ALX4)*SSDY(J)
       TALXZ=(TNUY(I,J,K+1)*ALX5-TNUY(I,J,K)*ALX6)*SSDZ(K)
       ALX=ALX+FLOAT(ILES)*(2.*TALXX+TALXY+TALXZ)
       RK3X=RK3X+FLOAT(ILES)*RHS1(I,J,K,1)!RHS1 in here: RHSSGS(lica_sgs.f90)
      ENDIF
!-----LES

       RHS1(I,J,K,1)=DT                                    &
               *( GAMMA(MSUB)*RK3X+RO(MSUB)*RK3XO(I,J,K)   &
               +2.*ALPHA*ALX                               &
               -2.*ALPHA*PRESSURE_GRADIENT(I,J,K,1)     &
               -2.*ALPHA*PMI*FLOAT(ICH) )
       RK3XOO(I,J,K)=RK3XO(I,J,K) ! added term for Fx
       RK3XO(I,J,K)=RK3X
       ENDDO
       ENDDO
       ENDDO

!-----RHS1 calculation for v momentum -----------------
!$OMP PARALLEL DO  &
!$OMP private(VE,VW,VC,VF,UEY,UWY,WCY,WFY)   &
!$OMP private(ANYN,ANYS,ANYY,ANY1,ANY2,RK3Y) &
!$OMP private(ALY1,ALY2,ALY3,ALY4,ALY5,ALY6) &
!$OMP private(ALYX,ALYY,ALYZ,ALY,TALYX,TALYY,TALYZ)
       DO K=1,N3M
       DO J=2,N2M
       DO I=1,N1M

       VE=(0.5*(SDX(I+1)*V(I,J,K)+SDX(I)*V(I+1,J,K))*VVDX(I+1))   &
          *(1.-FIXIU(I))+V(N1,J,K)*FIXIU(I)
       VW=(0.5*(SDX(I)*V(I-1,J,K)+SDX(I-1)*V(I,J,K))*VVDX(I))     &
          *(1.-FIXIL(I))+V(0,J,K)*FIXIL(I)
       VC=(0.5*(SDZ(K+1)*V(I,J,K)+SDZ(K)*V(I,J,K+1))*VVDZ(K+1))   &
          *(1.-FIXKU(K))+V(I,J,N3)*FIXKU(K)
       VF=(0.5*(SDZ(K)*V(I,J,K-1)+SDZ(K-1)*V(I,J,K))*VVDZ(K))     &
          *(1.-FIXKL(K))+V(I,J,0)*FIXKL(K)
       UEY=0.5*(SDY(J-1)*U(I+1,J,K)+SDY(J)*U(I+1,J-1,K))*VVDY(J)
       UWY=0.5*(SDY(J-1)*U(I,J,K)+SDY(J)*U(I,J-1,K))*VVDY(J)
       WCY=0.5*(SDY(J-1)*W(I,J,K+1)+SDY(J)*W(I,J-1,K+1))*VVDY(J)
       WFY=0.5*(SDY(J-1)*W(I,J,K)+SDY(J)*W(I,J-1,K))*VVDY(J)

       ANYN=0.5*(V(I,J+1,K)+V(I,J,K))
       ANYS=0.5*(V(I,J,K)+V(I,J-1,K))
       ANYY=VVDY(J)*(ANYN**2-ANYS**2)
       ANY1=(UEY*VE-UWY*VW)*SSDX(I)
       ANY2=(VC*WCY-VF*WFY)*SSDZ(K)
       RK3Y=-ANYY-ANY1-ANY2                ! Ny term

       ALY1=(V(I+1,J,K)-V(I,J,K))*VVDX(I+1)
       ALY2=(V(I,J,K)-V(I-1,J,K))*VVDX(I)
       ALY3=(V(I,J+1,K)-V(I,J,K))*SSDY(J)
       ALY4=(V(I,J,K)-V(I,J-1,K))*SSDY(J-1)
       ALY5=(V(I,J,K+1)-V(I,J,K))*VVDZ(K+1)*(1.-FIXKU(K)*FLOAT(KVT))
       ALY6=(V(I,J,K)-V(I,J,K-1))*VVDZ(K) *(1.-FIXKL(K)*FLOAT(KVB))
       ALYX=(ALY1-ALY2)*SSDX(I)
       ALYY=(ALY3-ALY4)*VVDY(J)
       ALYZ=(ALY5-ALY6)*SSDZ(K)
       ALY=REI*(ALYX+ALYY+ALYZ)

!-----LES
      IF (ILES.EQ.1) THEN
       TALYX=(TNUZ(I+1,J,K)*ALY1-TNUZ(I,J,K)*ALY2)*SSDX(I)
       TALYY=(TNU(I,J,K)*ALY3-TNU(I,J-1,K)*ALY4)*VVDY(J)
       TALYZ=(TNUX(I,J,K+1)*ALY5-TNUX(I,J,K)*ALY6)*SSDZ(K)
       ALY=ALY+FLOAT(ILES)*(TALYX+2.*TALYY+TALYZ)
       RK3Y=RK3Y+FLOAT(ILES)*RHS1(I,J,K,2)
      ENDIF
!-----LES

       RHS1(I,J,K,2)=DT                                    &
              *( GAMMA(MSUB)*RK3Y+RO(MSUB)*RK3YO(I,J,K)    &
                +2.*ALPHA*ALY                              &
                -2.*ALPHA*PRESSURE_GRADIENT(I,J,K,2))
       RK3YOO(I,J,K)=RK3YO(I,J,K)  ! added term for Fy
       RK3YO(I,J,K)=RK3Y

       ENDDO
       ENDDO
       ENDDO

       IF (N3M.EQ.1) GOTO 100
!-----RHS1 calculation for w momentum -----------------
!$OMP PARALLEL DO &
!$OMP private(WE,WW,WN,WS,UEZ,UWZ,VNZ,VSZ)    &
!$OMP private(ANZC,ANZF,ANZZ,ANZ1,ANZ2,RK3Z)  &
!$OMP private(ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6)  &
!$OMP private(ALZX,ALZY,ALZZ,ALZ,TALZX,TALZY,TALZZ)
       DO K=K_BGPZ,N3M
       DO J=1,N2M
       DO I=1,N1M

       WE=0.5*(SDX(I+1)*W(I,J,K)+SDX(I)*W(I+1,J,K))*VVDX(I+1)    &
           *(1.-FIXIU(I))+FIXIU(I)*W(N1,J,K)
       WW=0.5*(SDX(I)*W(I-1,J,K)+SDX(I-1)*W(I,J,K))*VVDX(I)      &
           *(1.-FIXIL(I))+FIXIL(I)*W(0,J,K)
       WN=(0.5*(SDY(J+1)*W(I,J,K)+SDY(J)*W(I,J+1,K))*VVDY(J+1))  &
           *(1.-FIXJU(J))+FIXJU(J)*W(I,N2,K)
       WS=(0.5*(SDY(J)*W(I,J-1,K)+SDY(J-1)*W(I,J,K))*VVDY(J))    &
           *(1.-FIXJL(J))+FIXJL(J)*W(I,0,K)
       UEZ=0.5*(SDZ(K-1)*U(I+1,J,K)+SDZ(K)*U(I+1,J,K-1))*VVDZ(K)
       UWZ=0.5*(SDZ(K-1)*U(I,J,K)+SDZ(K)*U(I,J,K-1))*VVDZ(K)
       VNZ=0.5*(SDZ(K-1)*V(I,J+1,K)+SDZ(K)*V(I,J+1,K-1))*VVDZ(K)
       VSZ=0.5*(SDZ(K-1)*V(I,J,K)+SDZ(K)*V(I,J,K-1))*VVDZ(K)

       ANZC=0.5*(W(I,J,K+1)+W(I,J,K))
       ANZF=0.5*(W(I,J,K)+W(I,J,K-1))
       ANZZ=(ANZC**2-ANZF**2)*VVDZ(K)
       ANZ1=(UEZ*WE-UWZ*WW)*SSDX(I)
       ANZ2=(VNZ*WN-VSZ*WS)*SSDY(J)
       RK3Z=-ANZZ-ANZ1-ANZ2                        ! Nz term

       ALZ1=(W(I+1,J,K)-W(I,J,K))*VVDX(I+1)
       ALZ2=(W(I,J,K)-W(I-1,J,K))*VVDX(I)
       ALZ3=(W(I,J+1,K)-W(I,J,K))*VVDY(J+1)*(1.-FIXJU(J)*FLOAT(JWT))
       ALZ4=(W(I,J,K)-W(I,J-1,K))*VVDY(J)  *(1.-FIXJL(J)*FLOAT(JWB))
       ALZ5=(W(I,J,K+1)-W(I,J,K))*SSDZ(K)
       ALZ6=(W(I,J,K)-W(I,J,K-1))*SSDZ(K-1)
       ALZX=(ALZ1-ALZ2)*SSDX(I)
       ALZY=(ALZ3-ALZ4)*SSDY(J)
       ALZZ=(ALZ5-ALZ6)*VVDZ(K)
       ALZ=REI*(ALZX+ALZY+ALZZ)

!-----LES
      IF (ILES.EQ.1) THEN
       TALZX=(TNUY(I+1,J,K)*ALZ1-TNUY(I,J,K)*ALZ2)*SSDX(I)
       TALZY=(TNUX(I,J+1,K)*ALZ3-TNUX(I,J,K)*ALZ4)*SSDY(J)
       TALZZ=(TNU(I,J,K)*ALZ5-TNU(I,J,K-1)*ALZ6)*VVDZ(K)
       ALZ=ALZ+FLOAT(ILES)*(TALZX+TALZY+2.*TALZZ)
       RK3Z=RK3Z+FLOAT(ILES)*RHS1(I,J,K,3)
      ENDIF
!-----LES

       RHS1(I,J,K,3)=DT                                   &
              *( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)   &
                +2.*ALPHA*ALZ                             &
                -2.*ALPHA*PRESSURE_GRADIENT(I,J,K,3))
       RK3ZOO(I,J,K)=RK3ZO(I,J,K)      ! added term for Fz
       RK3ZO(I,J,K)=RK3Z

       ENDDO
       ENDDO
       ENDDO

   100 RETURN
       END

!=======================================================================
      SUBROUTINE LHSU
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute streamwise velocity (u hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,TNU,TNUX,TNUY,TNUZ
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE=FLOAT(ILES)*RE
       CRE2=2.*FLOAT(ILES)*RE

!=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=1,N2M
      DO 141 K=1,N3M
      DO 141 I=I_BGPX,N1M
       AK(I,K)=AKUV(K)*(1.+CRE*TNUY(I,J,K)) *(1.-FIXKL(K)*FLOAT(KUB))
       CK(I,K)=CKUV(K)*(1.+CRE*TNUY(I,J,K+1))*(1.-FIXKU(K)*FLOAT(KUT))
       BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
       GK(I,K)=ACOEFI*RHS1(I,J,K,1)
  141 CONTINUE

      IF (IPZ .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,I_BGPX,N1M)
      ELSE IF (IPZ .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,I_BGPX,N1M)  !z periodicity
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=I_BGPX,N1M
      RHS1(I,J,K,1)=GK(I,K)
  151 CONTINUE
  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL  &
!$OMP private(AJ,CJ,BJ,GJ)  &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=1,N3M

!-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=I_BGPX,N1M
      AJ(I,J)=AJUW(J)*(1.+CRE*TNUZ(I,J,K))  *(1.-FIXJL(J)*FLOAT(JUB))
      CJ(I,J)=CJUW(J)*(1.+CRE*TNUZ(I,J+1,K))*(1.-FIXJU(J)*FLOAT(JUT))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,1)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,I_BGPX,N1M)

!-----X-DIRECTION
      DO 51 I=I_BGPX,N1M
      DO 51 J=1,N2M
      AI(J,I)=AIU(I)*(1.+CRE2*TNU(I-1,J,K))
      CI(J,I)=CIU(I)*(1.+CRE2*TNU(I,J,K))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (IPX .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,I_BGPX,N1M,1,N2M)
      ELSE IF (IPX .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,I_BGPX,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=I_BGPX,N1M
        U(I,J,K)=GI(J,I)+U(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE LHSV
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute transverse velocity (v hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,TNU,TNUX,TNUY,TNUZ
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE = FLOAT(ILES)*RE
       CRE2= 2.*FLOAT(ILES)*RE

!=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=2,N2M

      DO 141 K=1,N3M
      DO 141 I=1,N1M
      AK(I,K)=AKUV(K)*(1.+CRE*TNUX(I,J,K)) *(1.-FIXKL(K)*FLOAT(KVB))
      CK(I,K)=CKUV(K)*(1.+CRE*TNUX(I,J,K+1))*(1.-FIXKU(K)*FLOAT(KVT))
      BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
      GK(I,K)=ACOEFI*RHS1(I,J,K,2)
  141 CONTINUE

      IF (IPZ .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ELSE IF (IPZ .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,2)=GK(I,K)
  151 CONTINUE
  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL &
!$OMP private(AJ,CJ,BJ,GJ) &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=1,N3M

!-----Y-DIRECTION
      DO 91 J=2,N2M
      DO 91 I=1,N1M
      AJ(I,J)=AJV(J)*(1.+CRE2*TNU(I,J-1,K))
      CJ(I,J)=CJV(J)*(1.+CRE2*TNU(I,J,K))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,2)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,2,N2M,1,N1M)

!-----X-DIRECTION
      DO 51 I=1,N1M
      DO 51 J=2,N2M
      AI(J,I)=AIVW(I)*(1.+CRE*TNUZ(I,J,K))
      CI(J,I)=CIVW(I)*(1.+CRE*TNUZ(I+1,J,K))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (IPX .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,2,N2M)
      ELSE IF (IPX .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,2,N2M)
      ENDIF

      DO 61 J=2,N2M
      DO 61 I=1,N1M
      V(I,J,K)=GI(J,I)+V(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE LHSW
!=======================================================================
!
!     Calculate intermediate velocity, u_i hat through TDMA
!           In this routine, compute spanwise velocity (w hat)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,TNU,TNUX,TNUY,TNUZ
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AI,BI,CI,GI
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AJ,BJ,CJ,GJ
       REAL*8, DIMENSION (:,:), ALLOCATABLE :: AK,BK,CK,GK

       CRE = FLOAT(ILES)*RE
       CRE2= 2.*FLOAT(ILES)*RE

!=====ADI STARTS

!-----Z-DIRECTION
!$OMP PARALLEL &
!$OMP private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=1,N2M

      DO 141 K=K_BGPZ,N3M
      DO 141 I=1,N1M
      AK(I,K)=AKW(K)*(1.+CRE2*TNU(I,J,K-1))
      CK(I,K)=CKW(K)*(1.+CRE2*TNU(I,J,K))
      BK(I,K)=ACOEFI-AK(I,K)-CK(I,K)
      GK(I,K)=ACOEFI*RHS1(I,J,K,3)
  141 CONTINUE

      IF (IPZ .EQ. 0) THEN
         CALL TRDIAG3(AK,BK,CK,GK,GK,2,N3M,1,N1M)
      ELSE IF (IPZ .EQ. 1) THEN
         CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ENDIF

      DO 151 K=K_BGPZ,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,3)=GK(I,K)
  151 CONTINUE

  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL


!$OMP PARALLEL  &
!$OMP private(AJ,CJ,BJ,GJ)  &
!$OMP private(AI,CI,BI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=K_BGPZ,N3M

!-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=1,N1M
      AJ(I,J)=AJUW(J)*(1.+CRE*TNUX(I,J,K))  *(1.-FIXJL(J)*FLOAT(JWB))
      CJ(I,J)=CJUW(J)*(1.+CRE*TNUX(I,J+1,K))*(1.-FIXJU(J)*FLOAT(JWT))
      BJ(I,J)=ACOEFI-AJ(I,J)-CJ(I,J)
      GJ(I,J)=ACOEFI*RHS1(I,J,K,3)
   91 CONTINUE

      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,1,N1M)

!-----X-DIRECTION
      DO 51 I=1,N1M
      DO 51 J=1,N2M
      AI(J,I)=AIVW(I)*(1.+CRE*TNUY(I,J,K))
      CI(J,I)=CIVW(I)*(1.+CRE*TNUY(I+1,J,K))
      BI(J,I)=ACOEFI-AI(J,I)-CI(J,I)
      GI(J,I)=ACOEFI*GJ(I,J)
   51 CONTINUE

      IF (IPX .EQ. 0) THEN
         CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,1,N2M)
      ELSE IF (IPX .EQ. 1) THEN
         CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=1,N1M
      W(I,J,K)=GI(J,I)+W(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

!=======================================================================
      SUBROUTINE RETRV_UVW
      CALL AMBIENT_PREDICT
      RETURN
      END
!=======================================================================
       SUBROUTINE RHSINCORPBC
!=======================================================================
!
!     Applying exit boundary condition & top boundary condition to RHS
!           [RHS incorporate B.C.]
!
!     ACOEF: main solver in this file
!     CIU, CIVW, CJV: subroutine LHSINIT in this file
!     DUOUT, DVOUT, DVTOP, DWOUT: subroutine CONVBC (lica_cylinder.f90)
!     TNU, TNUZ, TNUY: subroutine NUTINTERPOL (lica_sgs.f90)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,RHS1,TNU,TNUX,TNUY,TNUZ
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: CRE,CRE2

       CRE=FLOAT(ILES)*RE
       CRE2=2.*FLOAT(ILES)*RE

!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=1,N2M
       RHS1(N1M,J,K,1)=RHS1(N1M,J,K,1)                                  &
                      -ACOEF*CIU(N1M)*(1.+CRE2*TNU(N1M,J,K))*DUOUT(J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=1,N3M
       DO J=2,N2M
       RHS1(N1M,J,K,2)=RHS1(N1M,J,K,2)                                  &
                      -ACOEF*CIVW(N1M)*(1.+CRE*TNUZ(N1,J,K))*DVOUT(J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=K_BGPZ,N3M
       DO J=1,N2M
       RHS1(N1M,J,K,3)=RHS1(N1M,J,K,3)                                  &
                      -ACOEF*CIVW(N1M)*(1.+CRE*TNUY(N1,J,K))*DWOUT(J,K)
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO K=1,N3M
       DO I=1,N1M
       RHS1(I,N2M,K,2)=RHS1(I,N2M,K,2)                                  &
                      -ACOEF*CJV(N2M)*(1.+CRE2*TNU(I,N2M,K))*DVTOP(I,K)
       ENDDO
       ENDDO

       RETURN
       END

!=======================================================================
      SUBROUTINE DIVGS(DIVGSUM)
!=======================================================================
!
!     RHS OF POISSON EQUATION
!
!     Options
!       MASSON = 0, IMOVINGON = 0: w/o IBM body
!       MASSON = 1, IMOVINGON = 0: w/  stationary body
!       MASSON = 1, IMOVINGON = 1: w/  moving body
!
!     Reference for IBM
!     Kim, Kim & Choi, 2001, An immersed-boundary finite-volume method
!       for simulations of flow in complex geometries, J. Comp. Phys.,
!       171, 132-150.
!
!     Variables
!       DIVGSUM: divergence of velocity
!       QMASS: mass source/sink for IBM body
!       UBD,VBD,WBD: translational velocity of body (IMOVINGON=1)
!                    These velocities are defined in lica_movingIBM.f90.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,QMASS,IFC,JFC,KFC
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: EXTENSION_VELOCITY
      IMPLICIT NONE
      INTEGER*8    :: I,J,K
      INTEGER*8    :: II,JJ,KK,IMM,JMM,KMM,N
      REAL*8       :: DIVG1,DIVG2,DIVG3
      REAL*8       :: DIVGSUM(0:M1,0:M2,0:M3)
      REAL*8       :: USTMP,VSTMP,WSTMP
      REAL*8       :: UBD,VBD,WBD

!############################################################### W/O IBM
!$OMP PARALLEL DO  &
!$OMP private(DIVG1,DIVG2,DIVG3)
      DO 73 K=1,N3M
      DO 73 J=1,N2M
      DO 73 I=1,N1M
        DIVG1=(U(I+1,J,K)-U(I,J,K))*SSDX(I)
        DIVG2=(V(I,J+1,K)-V(I,J,K))*SSDY(J)
        DIVG3=(W(I,J,K+1)-W(I,J,K))*SSDZ(K)
        DIVGSUM(I,J,K)=(DIVG1+DIVG2+DIVG3)
   73 CONTINUE

!############################################### W/  STATIONALY IBM body
      IF (MASSON.EQ.1 .AND. IMOVINGON.EQ.0) THEN

      QMASS = 0.

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        QMASS(II,JJ,KK) = QMASS(II,JJ,KK)-U(II,JJ,KK)*SSDX(II)
        DIVGSUM(II,JJ,KK) = DIVGSUM(II,JJ,KK)+U(II,JJ,KK)*SSDX(II)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,IMM)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        IMM=IMV(II)
        QMASS(IMM,JJ,KK)=QMASS(IMM,JJ,KK)+U(II,JJ,KK)*SSDX(IMM)
        DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-U(II,JJ,KK)*SSDX(IMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-V(II,JJ,KK)*SSDY(JJ)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+V(II,JJ,KK)*SSDY(JJ)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,JMM)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        JMM=JMV(JJ)
        QMASS(II,JMM,KK)=QMASS(II,JMM,KK)+V(II,JJ,KK)*SSDY(JMM)
        DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-V(II,JJ,KK)*SSDY(JMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-W(II,JJ,KK)*SSDZ(KK)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+W(II,JJ,KK)*SSDZ(KK)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,KMM)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        KMM=KMV(KK)
        QMASS(II,JJ,KMM)=QMASS(II,JJ,KMM)+W(II,JJ,KK)*SSDZ(KMM)
        DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-W(II,JJ,KK)*SSDZ(KMM)
      ENDDO

!################################################### W/  MOVING IBM body
      ELSEIF (MASSON.EQ.1 .AND. IMOVINGON.EQ.1) THEN

      QMASS = 0.

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,USTMP)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)

        USTMP=U(II,JJ,KK)-EXTENSION_VELOCITY(N,1)

        QMASS (II,JJ,KK) = QMASS  (II,JJ,KK)-USTMP*SSDX(II)
        DIVGSUM(II,JJ,KK)= DIVGSUM(II,JJ,KK)+USTMP*SSDX(II)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,IMM,USTMP)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        IMM=IMV(II)

        USTMP=U(II,JJ,KK)-EXTENSION_VELOCITY(N,1)

        QMASS  (IMM,JJ,KK)=QMASS  (IMM,JJ,KK)+USTMP*SSDX(IMM)
        DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-USTMP*SSDX(IMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,VSTMP)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)

        VSTMP=V(II,JJ,KK)-EXTENSION_VELOCITY(N,2)

        QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-VSTMP*SSDY(JJ)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+VSTMP*SSDY(JJ)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,JMM,VSTMP)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        JMM=JMV(JJ)

        VSTMP=V(II,JJ,KK)-EXTENSION_VELOCITY(N,2)

        QMASS  (II,JMM,KK)=QMASS  (II,JMM,KK)+VSTMP*SSDY(JMM)
        DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-VSTMP*SSDY(JMM)
      ENDDO

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,WSTMP)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)

        WSTMP=W(II,JJ,KK)-EXTENSION_VELOCITY(N,3)

        QMASS  (II,JJ,KK)=QMASS  (II,JJ,KK)-WSTMP*SSDZ(KK)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+WSTMP*SSDZ(KK)
      ENDDO
!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK,KMM,WSTMP)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        KMM=KMV(KK)

        WSTMP=W(II,JJ,KK)-EXTENSION_VELOCITY(N,3)

        QMASS  (II,JJ,KMM)=QMASS  (II,JJ,KMM)+WSTMP*SSDZ(KMM)
        DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-WSTMP*SSDZ(KMM)
      ENDDO

      ENDIF

!$OMP PARALLEL DO
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
      DIVGSUM(I,J,K)=DIVGSUM(I,J,K)*DTCONSTI
   80 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE UCALC(PHI)
!=======================================================================
!
!     Calculate velocity (u_i) from u_i hat
!     u_i hat is derived from pseudo-pressure, phi.
!
!     Definition of pseudo-pressure, phi, is as follows:
!     ({u_i}^k - u_i hat)/(2.*alpha_k*dt) = - d({phi}^k)/d(x_i)
!
!     From the definition of phi, following equation is derived:
!     {u_i}^k = u_i hat - 2.*a_k*DT*d({phi}^k)/d(x_i)
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       INTEGER*8     :: IM,KM
       REAL*8        :: PHI(0:M1,0:M2,0:M3)

!$OMP PARALLEL DO &
!$OMP private(IM)
      DO 21 K=1,N3M
      DO 21 J=1,N2M
      DO 21 I=I_BGPX,N1M                  ! I=1,N1 => BOUNDARY
         IM=IMV(I)
         U(I,J,K)=U(I,J,K)                                   &
                 -DTCONST*(PHI(I,J,K)-PHI(IM,J,K))*VVDX(I)   &
                 -PHCAP*FLOAT(ICH)
   21 CONTINUE

!$OMP PARALLEL DO
      DO 31 K=1,N3M
      DO 31 J=2,N2M                       ! J=1,N2 => BOUNDARY
      DO 31 I=1,N1M
         V(I,J,K)=V(I,J,K)-DTCONST*(PHI(I,J,K)-PHI(I,J-1,K))*VVDY(J)
   31 CONTINUE

!$OMP PARALLEL DO &
!$OMP private(KM)
      DO 41 K=K_BGPZ,N3M                  ! K=0,N3 => BOUNDARY
         KM=KMV(K)
      DO 41 J=1,N2M
      DO 41 I=1,N1M
         W(I,J,K)=W(I,J,K)-DTCONST*(PHI(I,J,K)-PHI(I,J,KM))*VVDZ(K)
   41 CONTINUE

      CALL AMBIENT_PROJECT(PHI)
!!!!!!     NEUMANN & DIRICHLET BOUNDARY CONDITION
      IF (JUT.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=0,N3
         DO I=1,N1
            U(I,N2,K)=U(I,N2M,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JWT.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=1,N3
         DO I=0,N1
            W(I,N2,K)=W(I,N2M,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JUB.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=0,N3
         DO I=1,N1
            U(I,0,K)=U(I,1,K)
         ENDDO
         ENDDO
      ENDIF
      IF (JWB.EQ.1) THEN
!$OMP PARALLEL DO
         DO K=1,N3
         DO I=0,N1
            W(I,0,K)=W(I,1,K)
         ENDDO
         ENDDO
      ENDIF
      IF (KUT.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=0,N2
         DO I=1,N1
            U(I,J,N3)=U(I,J,N3M)
         ENDDO
         ENDDO
      ENDIF
      IF (KVT.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=1,N2
         DO I=0,N1
            V(I,J,N3)=V(I,J,N3M)
         ENDDO
         ENDDO
      ENDIF
      IF (KUB.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=0,N2
         DO I=1,N1
            U(I,J,0)=U(I,J,1)
         ENDDO
         ENDDO
      ENDIF
      IF (KVB.EQ.1) THEN
!$OMP PARALLEL DO
         DO J=1,N2
         DO I=0,N1
            V(I,J,0)=V(I,J,1)
         ENDDO
         ENDDO
      ENDIF

!     Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=0,N2
      DO I=1,N1
         U(I,J,0) =U(I,J,N3M)
         U(I,J,N3)=U(I,J,1)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO J=1,N2
      DO I=0,N1
         V(I,J,0) =V(I,J,N3M)
         V(I,J,N3)=V(I,J,1)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO J=0,N2
      DO I=0,N1
         W(I,J,0) =W(I,J,N3M)
         W(I,J,N3)=W(I,J,1)
      ENDDO
      ENDDO
      ENDIF

!     X PERIODICITY
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=0,N3
      DO J=0,N2
         U(0 ,J,K)=U(N1M,J,K)
         U(N1,J,K)=U(1  ,J,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO K=0,N3
      DO J=1,N2
         V(0 ,J,K)=V(N1M,J,K)
         V(N1,J,K)=V(1  ,J,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO
      DO K=1,N3
      DO J=0,N2
         W(0 ,J,K)=W(N1M,J,K)
         W(N1,J,K)=W(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE PCALC(PHI,DIVGSUM)
!=======================================================================
!
!     Calculate pressure from DIVGSUM & PHI
!
!     ({u_i}^k-u_i hat)/(2.*a_k*DT)
!                             = - [ d(p^k)/d(x_i) - d(p^(k-1))/d(x_i) ]
!                               + 0.5*[L({u_i}^k) - L(u_i hat)]
!                             = - d({phi}^k)/d(x_i)
!
!     By the definition of phi, final equation is derived as follows:
!       p^k = p^(k-1) + {phi}^k - (a_k*DT/Re)*(d2({phi^k})/d(x_j)d(x_j))
!
!-----------------------------------------------------------------------
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,P
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8        :: PHIREF
       REAL*8        :: PHI(0:M1,0:M2,0:M3),DIVGSUM(0:M1,0:M2,0:M3)


!$OMP PARALLEL DO
      DO 92 K=1,N3M
      DO 92 J=1,N2M
      DO 92 I=1,N1M
      P(I,J,K)=P(I,J,K)+PHI(I,J,K)-0.5*DTCONST*DIVGSUM(I,J,K)*REI
   92 CONTINUE

      IF(ICH.EQ.1) THEN
!     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
      PHIREF = 0.
!$OMP PARALLEL DO &
!$OMP reduction(+:PHIREF)
      DO 80 K=1,N3M
      DO 80 I=1,N1M
      PHIREF=PHI(I,N2M,K)*SDX(I)*SDZ(1)+PHIREF
   80 CONTINUE
      PHIREF=PHIREF/XL/ZL

!$OMP PARALLEL DO
      DO 93 K=1,N3M
      DO 93 J=1,N2M
      DO 93 I=1,N1M
      P(I,J,K)=P(I,J,K)-PHIREF
   93 CONTINUE
      ENDIF

!     Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO
      DO J=1,N2M
      DO I=1,N1M
         P(I,J,0) =P(I,J,N3M)
         P(I,J,N3)=P(I,J,1)
      ENDDO
      ENDDO
      ENDIF

!     X PERIODICITY
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO
      DO K=1,N3M
      DO J=1,N2M
         P(0 ,J,K)=P(N1M,J,K)
         P(N1,J,K)=P(1  ,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG1(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x1-direction (TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!       At each i,
!       |B_{L1}   C_{L1}   0        0               | = |R_{L1  }|
!       |A_{L1+1} B_{L1+1} C_{L1+1} 0               | = |R_{L1+1}|
!       |0        A_{L1+2} B_{L1+2} C_{L1+2}        | = |R_{L1+2}|
!       |0                                          | = |  :     |
!       |0                                          | = |  :     |
!       |0                             A_{L2} B_{L2}| = |R_{L2}  |
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,L1,L2,LL1,LL2
      REAL*8    :: GAM(M2,M1),A(M2,M1),B(M2,M1),C(M2,M1),R(M2,M1),UU(M2,M1),BET(M2)

      DO 10 I=LL1,LL2
      BET(I)  = 1./B(I,L1)
      UU(I,L1)= R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)= C(I,J-1)*BET(I)
      BET(I)  = 1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J) = (R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J) = UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG2(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-2 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x2-direction (TDMA)
!       LL2,LL2: 1st index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J
      INTEGER*8 :: L1,L2,LL1,LL2
      REAL*8    :: GAM(M1,M2),A(M1,M2),B(M1,M2),C(M1,M2),R(M1,M2),UU(M1,M2),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG3(A,B,C,R,UU,L1,L2,LL1,LL2)
!=======================================================================
!
!     SOLVE THE TRIDIAGONAL MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x3-direction (TDMA)
!       LL2,LL2: 1st index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J
      INTEGER*8 :: L1,L2,LL1,LL2
      REAL*8    :: GAM(M1,M3),A(M1,M3),B(M1,M3),C(M1,M3),R(M1,M3),UU(M1,M3),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE

      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG1P(A,B,C,F,J1,J2,L1,L2)
!=======================================================================
!
!     If IPX = 1 (X periodicity),
!     INVERT A PERIODIC MATRIX (X-1 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x1-direction (periodic TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K,J1,J2,L1,L2,JA,JJ
      REAL*8    :: A(M2,M1),B(M2,M1),C(M2,M1),F(M2,M1)
      REAL*8    :: Q(M2,M1),S(M2,M1),QE(M2,M1),FN(M2),PN(M2)
      REAL*8    :: BINV

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV= 1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

!     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

!     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2) = 1.
      QE(K,J2)= 0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1)) &
             /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

!     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE TRDIAG3P(A,B,C,F,J1,J2,L1,L2)
!=======================================================================
!
!     If IPZ = 1 (Z periodicity),
!     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
!     DIFFERENT COEFFICIENTS in LHSU,LHSV,LHSW subroutines
!
!     Variables
!       A(),B(),C()  : Coefficients of tridiagonal matrix
!       R()          : RHS
!       UU()         : Solution
!       L1,L2  : 2nd index indicating x3-direction (periodic TDMA)
!       LL2,LL2: 1st index indicating x2-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      INTEGER*8 :: J1,J2,L1,L2,JA,JJ
      REAL*8    :: A(M1,M3),B(M1,M3),C(M1,M3),F(M1,M3)
      REAL*8    :: Q(M1,M3),S(M1,M3),QE(M1,M3),FN(M1),PN(M1)
      REAL*8    :: BINV

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

!     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

!     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))  &
             /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

!     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE PTDIAG1A(A,B,C,D,E,F,UU,I1,IN,J1,JN)
!=======================================================================
!
!     SOLVE THE PENTADIAGONAL MATRIX (X-1 PENTADIAGONAL)
!
!     Variables
!       A(),B(),C(),D(),E(),F(): Coefficients of tridiagonal matrix
!       UU()   : RHS and solution
!       I1,IN  : 1st index indicating x2-direction
!       J1,JN  : 2nd index indicating x1-direction
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8 :: I,J,K
      INTEGER*8 :: I1,IN,J1,JN,J2,J3,JNM,JM,JMM
      REAL*8    ::  A(M2,M1),B(M2,M1),C(M2,M1),D(M2,M1),E(M2,M1)
      REAL*8    ::  O(M2,M1),Q(M2,M1),R(M2,M1),F(M2,M1),UU(M2,M1)
      REAL*8    ::  PDENO2(M2),PDENOJ(M2)

      J2=J1+1
      J3=J2+1
      JNM=JN-1
      DO 1 I=I1,IN
      O(I,J1)=-B(I,J1)/A(I,J1)
      Q(I,J1)=-C(I,J1)/A(I,J1)
      R(I,J1)=F(I,J1)/A(I,J1)
    1 CONTINUE

      DO 2 I=I1,IN
      PDENO2(I)=1./(A(I,J2)+D(I,J2)*O(I,J1))
      O(I,J2)=-(B(I,J2)+D(I,J2)*Q(I,J1))*PDENO2(I)
      Q(I,J2)=-C(I,J2)*PDENO2(I)
      R(I,J2)=(F(I,J2)-D(I,J2)*R(I,J1))*PDENO2(I)
    2 CONTINUE

      DO 10 J=J3,JN
      JM=J-1
      JMM=JM-1
      DO 10 I=I1,IN
      PDENOJ(I)=1./(A(I,J)+E(I,J)*Q(I,JMM)                            &
                  +(D(I,J)+E(I,J)*O(I,JMM))*O(I,JM))
      O(I,J)=-(B(I,J)+(D(I,J)+E(I,J)*O(I,JMM))*Q(I,JM))*PDENOJ(I)
      Q(I,J)=-C(I,J)*PDENOJ(I)
      R(I,J)=(F(I,J)-E(I,J)*R(I,JMM)-(D(I,J)+E(I,J)*O(I,JMM))*R(I,JM)) &
            *PDENOJ(I)
   10 CONTINUE

      DO 11 I=I1,IN
      UU(I,JN)=R(I,JN)
   11 CONTINUE

      DO 12 I=I1,IN
      UU(I,JNM)=O(I,JNM)*UU(I,JN)+R(I,JNM)
   12 CONTINUE

      DO 20 J=JNM-1,J1,-1
      DO 20 I=I1,IN
      UU(I,J)=O(I,J)*UU(I,J+1)+Q(I,J)*UU(I,J+2)+R(I,J)
   20 CONTINUE

      RETURN
      END

!=======================================================================
      SUBROUTINE FIELD_AVG_INIT
!=======================================================================
!
!     Subroutine for initializing the variables for time-averaging.
!     ...most of which are required to obtain turbulence statistics
!     Refer to mod_common.f90 for meaning of each variable
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : UAVG,VAVG,WAVG,UIUJAVG,PAVG,P2AVG,VORAVG,WRMS,SSAVG
      IMPLICIT NONE

      UAVG   = 0.
      VAVG   = 0.
      WAVG   = 0.
      UIUJAVG= 0.
      PAVG   = 0.
      P2AVG  = 0.
      SSAVG  = 0.
      VORAVG = 0.
      WRMS   = 0.

      TIMEINIT = TIME
      IHISTINIT= IHIST

      RETURN
      END

!=======================================================================
      SUBROUTINE FIELD_AVG
!=======================================================================
!     SUBROUTINE FOR FIELD AVERAGING
!     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
!     VARIABLES ARE LINEARLY INTERPOLATED
!     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES
!
!     No spatial averaging.
!
!     Print an averaged flow-field file, when MOD(NTIME,NPRIAVG)=0,
!     satisfying the averaged flow-field file printing interval (NPRIAVG set in lica.in).
!     tfn1='fav': prefix for average flow field
!     tname     : average field file name  ex) fav100000-110000
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P,UAVG,VAVG,WAVG,UIUJAVG,PAVG,P2AVG,VORAVG,WRMS,SSAVG
      IMPLICIT NONE
      INTEGER*8     :: I,J,K,L
      REAL*8        :: UCC,VCC,WCC
      REAL*8        :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8        :: UP,UM,VP,VM,WP,WM,SR(6),SRSR
      REAL*8        :: WX,WY,WZ
      REAL*8        :: TIMEEND
      INTEGER*8     :: IHISTEND
      INTEGER*8     :: idg1,idg2,idg3,idg4,idg5,idg6
      CHARACTER*16  :: tname
      CHARACTER*3   :: tfn1
      CHARACTER*6   :: tfn2,tfn3
      CHARACTER*1   :: tfnh

      tfnh='-'

!$OMP PARALLEL DO &
!$OMP private(UCC,VCC,WCC)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       UCC=0.5*(U(I,J,K)+U(I+1,J,K))
       VCC=0.5*(V(I,J,K)+V(I,J+1,K))
       WCC=0.5*(W(I,J,K)+W(I,J,K+1))
         UAVG(I,J,K)     =UAVG(I,J,K)     +UCC         *DT
         VAVG(I,J,K)     =VAVG(I,J,K)     +VCC         *DT
         WAVG(I,J,K)     =WAVG(I,J,K)     +WCC         *DT
         UIUJAVG(I,J,K,1)=UIUJAVG(I,J,K,1)+UCC**2.     *DT
         UIUJAVG(I,J,K,2)=UIUJAVG(I,J,K,2)+UCC*VCC     *DT
         UIUJAVG(I,J,K,3)=UIUJAVG(I,J,K,3)+UCC*WCC     *DT
         UIUJAVG(I,J,K,4)=UIUJAVG(I,J,K,4)+VCC**2.     *DT
         UIUJAVG(I,J,K,5)=UIUJAVG(I,J,K,5)+VCC*WCC     *DT
         UIUJAVG(I,J,K,6)=UIUJAVG(I,J,K,6)+WCC**2.     *DT
         PAVG(I,J,K)     =PAVG(I,J,K)     +P(I,J,K)    *DT
         P2AVG(I,J,K)    =P2AVG(I,J,K)    +P(I,J,K)**2.*DT
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33) &
!$OMP private(UP,UM,VP,VM,WP,WM,SR,SRSR) &
!$OMP private(WX,WY,WZ)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M

       VG11=SSDX(I)*(U(I+1,J,K)-U(I,J,K))
       UP=VVDY(J+1)*0.25                                      &
          *(SDY(J+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +SDY(J)*(U(I,J+1,K)+U(I+1,J+1,K)))*(1.-FIXJU(J))   &
          +0.5*(U(I,N2,K)+U(I+1,N2,K))*FIXJU(J)
       UM=VVDY(J)*0.25                                        &
          *(SDY(J)*(U(I,J-1,K)+U(I+1,J-1,K))                  &
           +SDY(J-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXJL(J))     &
          +0.5*(U(I,0,K)+U(I+1,0,K))*FIXJL(J)
       VG12=SSDY(J)*(UP-UM)

       UP=VVDZ(K+1)*0.25                                      &
          *(SDZ(K+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +SDZ(K)*(U(I,J,K+1)+U(I+1,J,K+1)))*(1.-FIXKU(K))   &
          +0.5*(U(I,J,N3)+U(I+1,J,N3))*FIXKU(K)
       UM=VVDZ(K)*0.25                                        &
          *(SDZ(K)*(U(I,J,K-1)+U(I+1,J,K-1))                  &
           +SDZ(K-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXKL(K))     &
          +0.5*(U(I,J,0)+U(I+1,J,0))*FIXKL(K)
       VG13=SSDZ(K)*(UP-UM)

       VP=VVDX(I+1)*0.25                                      &
          *(SDX(I+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +SDX(I)*(V(I+1,J,K)+V(I+1,J+1,K)))*(1.-FIXIU(I))   &
          +0.5*(V(N1,J,K)+V(N1,J+1,K))*FIXIU(I)
       VM=VVDX(I)*0.25                                        &
          *(SDX(I)*(V(I-1,J,K)+V(I-1,J+1,K))                  &
           +SDX(I-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXIL(I))     &
          +0.5*(V(0,J,K)+V(0,J+1,K))*FIXIL(I)
       VG21=SSDX(I)*(VP-VM)

       VG22=SSDY(J)*(V(I,J+1,K)-V(I,J,K))

       VP=VVDZ(K+1)*0.25                                      &
          *(SDZ(K+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +SDZ(K)*(V(I,J,K+1)+V(I,J+1,K+1)))*(1.-FIXKU(K))   &
          +0.5*(V(I,J,N3)+V(I,J+1,N3))*FIXKU(K)
       VM=VVDZ(K)*0.25                                        &
          *(SDZ(K)*(V(I,J,K-1)+V(I,J+1,K-1))                  &
           +SDZ(K-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXKL(K))     &
          +0.5*(V(I,J,0)+V(I,J+1,0))*FIXKL(K)
       VG23=SSDZ(K)*(VP-VM)

       WP=VVDX(I+1)*0.25                                      &
          *(SDX(I+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +SDX(I)*(W(I+1,J,K)+W(I+1,J,K+1)))*(1.-FIXIU(I))   &
          +0.5*(W(N1,J,K)+W(N1,J,K+1))*FIXIU(I)
       WM=VVDX(I)*0.25                                        &
          *(SDX(I)*(W(I-1,J,K)+W(I-1,J,K+1))                  &
           +SDX(I-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXIL(I))     &
          +0.5*(W(0,J,K)+W(0,J,K+1))*FIXIL(I)
       VG31=SSDX(I)*(WP-WM)

       WP=VVDY(J+1)*0.25                                      &
          *(SDY(J+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +SDY(J)*(W(I,J+1,K)+W(I,J+1,K+1)))*(1.-FIXJU(J))   &
          +0.5*(W(I,N2,K)+W(I,N2,K+1))*FIXJU(J)
       WM=VVDY(J)*0.25                                        &
          *(SDY(J)*(W(I,J-1,K)+W(I,J-1,K+1))                  &
           +SDY(J-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXJL(J))     &
          +0.5*(W(I,0,K)+W(I,0,K+1))*FIXJL(J)
       VG32=SSDY(J)*(WP-WM)

       VG33=SSDZ(K)*(W(I,J,K+1)-W(I,J,K))

       SR(1)= VG11
       SR(2)= 0.5*(VG12+VG21)
       SR(3)= 0.5*(VG13+VG31)
       SR(4)= VG22
       SR(5)= 0.5*(VG23+VG32)
       SR(6)= VG33

       WX=VG32-VG23
       WY=VG13-VG31
       WZ=VG21-VG12

       VORAVG(I,J,K,1)=VORAVG(I,J,K,1)+WX*DT
       VORAVG(I,J,K,2)=VORAVG(I,J,K,2)+WY*DT
       VORAVG(I,J,K,3)=VORAVG(I,J,K,3)+WZ*DT

       WRMS(I,J,K,1)=WRMS(I,J,K,1)+WX*WX*DT
       WRMS(I,J,K,2)=WRMS(I,J,K,2)+WX*WY*DT
       WRMS(I,J,K,3)=WRMS(I,J,K,3)+WX*WZ*DT
       WRMS(I,J,K,4)=WRMS(I,J,K,4)+WY*WY*DT
       WRMS(I,J,K,5)=WRMS(I,J,K,5)+WY*WZ*DT
       WRMS(I,J,K,6)=WRMS(I,J,K,6)+WZ*WZ*DT

       SRSR= 2.*SR(2)**2.+SR(1)**2.+2.*SR(3)**2.+SR(4)**2.+2.*SR(5)**2.+SR(6)**2.
       SSAVG(I,J,K)= SSAVG(I,J,K)+SRSR*DT

      ENDDO
      ENDDO
      ENDDO

      IF (MOD(NTIME,NPRIAVG).EQ.0) THEN

       TIMEEND =TIME
       IHISTEND=IHIST

       tfn1='fav'
       idg1=ihistinit/100000
       idg2=(ihistinit-idg1*100000)/10000
       idg3=(ihistinit-idg1*100000-idg2*10000)/1000
       idg4=(ihistinit-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn2=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       idg1=ihist/100000
       idg2=(ihist-idg1*100000)/10000
       idg3=(ihist-idg1*100000-idg2*10000)/1000
       idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn3=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       tname=tfn1//tfn2//tfnh//tfn3

       OPEN(NAV,FILE=tname,FORM='UNFORMATTED')
       WRITE(NAV)N1M,N2M,N3M,RE
       WRITE(NAV)TIMEINIT,TIMEEND,IHISTINIT,IHISTEND
       WRITE(NAV)((( UAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( VAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( WAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((((UIUJAVG(I,J,K,L),I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
       WRITE(NAV)((( PAVG(I,J,K)     ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((( P2AVG(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
       WRITE(NAV)((((VORAVG(I,J,K,L) ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,3)
       WRITE(NAV)((((WRMS(I,J,K,L)   ,I=1,N1M),J=1,N2M),K=1,N3M),L=1,6)
       WRITE(NAV)((( SSAVG(I,J,K)    ,I=1,N1M),J=1,N2M),K=1,N3M)
       CLOSE(NAV)

       NAV      = NAV+1
       UAVG     = 0.
       VAVG     = 0.
       WAVG     = 0.
       UIUJAVG  = 0.
       PAVG     = 0.
       P2AVG    = 0.
       VORAVG   = 0.
       WRMS     = 0.
       SSAVG    = 0.
       TIMEINIT = TIME
       IHISTINIT= IHIST

      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE FIELD_AVG_INITZ
!=======================================================================
!     Subroutine for initializing the variables for time-averaging.
!     ...most of which are required to obtain turbulence statistics
!     Refer to mod_common.f90 for meaning of each variable
!
!     Used if the flow is homogeneous in the z-direction
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : UAVGZ,VAVGZ,WAVGZ,UIUJAVGZ,PAVGZ    &
                               ,P2AVGZ,VORAVGZ,WRMSZ,SSAVGZ,SIJAVGZ
      IMPLICIT NONE

      UAVGZ     = 0.
      VAVGZ     = 0.
      WAVGZ     = 0.
      UIUJAVGZ  = 0.
      PAVGZ     = 0.
      P2AVGZ    = 0.
      VORAVGZ   = 0.
      WRMSZ     = 0.
      SSAVGZ    = 0.
      SIJAVGZ   = 0.

      TIMEINITZ = TIME
      IHISTINITZ= IHIST

      RETURN
      END
!=======================================================================
      SUBROUTINE FIELD_AVGZ
!=======================================================================
!     SUBROUTINE FOR FIELD AVERAGING
!     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
!     VARIABLES ARE LINEARLY INTERPOLATED
!     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES

!     Spatial averaging in z-direction.
!
!     Print an averaged flow-field file, when MOD(NTIME,NPRIAVG)=0,
!     satisfying the averaged flow-field file printing interval (NPRIAVG set in lica.in).
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P,UAVGZ,VAVGZ,WAVGZ,UIUJAVGZ,PAVGZ    &
                               ,P2AVGZ,VORAVGZ,WRMSZ,SSAVGZ,SIJAVGZ
      IMPLICIT NONE
      INTEGER*8     :: I,J,K,L
      REAL*8        :: UCC,VCC,WCC
      REAL*8        :: VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33
      REAL*8        :: UP,UM,VP,VM,WP,WM,SR(6),SRSR
      REAL*8        :: WX,WY,WZ
      REAL*8        :: TIMEENDZ
      INTEGER*8     :: IHISTENDZ
      INTEGER*8     :: idg1,idg2,idg3,idg4,idg5,idg6
      CHARACTER*17  :: tname
      CHARACTER*4   :: tfn1
      CHARACTER*6   :: tfn2,tfn3
      CHARACTER*1   :: tfnh

      tfnh='-'

!$OMP PARALLEL DO &
!$OMP private(UCC,VCC,WCC)
      DO J=1,N2M
      DO K=1,N3M
      DO I=1,N1M
      UCC = 0.5*(U(I,J,K)+U(I+1,J,K))
      VCC = 0.5*(V(I,J,K)+V(I,J+1,K))
      WCC = 0.5*(W(I,J,K)+W(I,J,K+1))
         UAVGZ(I,J)     = UAVGZ(I,J)     + UCC         *DT
         VAVGZ(I,J)     = VAVGZ(I,J)     + VCC         *DT
         WAVGZ(I,J)     = WAVGZ(I,J)     + WCC         *DT
         UIUJAVGZ(I,J,1)= UIUJAVGZ(I,J,1)+ UCC**2.     *DT
         UIUJAVGZ(I,J,2)= UIUJAVGZ(I,J,2)+ UCC*VCC     *DT
         UIUJAVGZ(I,J,3)= UIUJAVGZ(I,J,3)+ UCC*WCC     *DT
         UIUJAVGZ(I,J,4)= UIUJAVGZ(I,J,4)+ VCC**2.     *DT
         UIUJAVGZ(I,J,5)= UIUJAVGZ(I,J,5)+ VCC*WCC     *DT
         UIUJAVGZ(I,J,6)= UIUJAVGZ(I,J,6)+ WCC**2.     *DT
         PAVGZ(I,J)     = PAVGZ(I,J)     + P(I,J,K)    *DT
         P2AVGZ(I,J)    = P2AVGZ(I,J)    + P(I,J,K)**2.*DT
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO &
!$OMP private(VG11,VG12,VG13,VG21,VG22,VG23,VG31,VG32,VG33) &
!$OMP private(UP,UM,VP,VM,WP,WM,SR,SRSR) &
!$OMP private(WX,WY,WZ)
      DO J=1,N2M
      DO K=1,N3M
      DO I=1,N1M

       VG11=SSDX(I)*(U(I+1,J,K)-U(I,J,K))
       UP=VVDY(J+1)*0.25                                      &
          *(SDY(J+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +SDY(J)*(U(I,J+1,K)+U(I+1,J+1,K)))*(1.-FIXJU(J))   &
          +0.5*(U(I,N2,K)+U(I+1,N2,K))*FIXJU(J)
       UM=VVDY(J)*0.25                                        &
          *(SDY(J)*(U(I,J-1,K)+U(I+1,J-1,K))                  &
           +SDY(J-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXJL(J))     &
          +0.5*(U(I,0,K)+U(I+1,0,K))*FIXJL(J)
       VG12=SSDY(J)*(UP-UM)

       UP=VVDZ(K+1)*0.25                                      &
          *(SDZ(K+1)*(U(I,J,K)+U(I+1,J,K))                    &
           +SDZ(K)*(U(I,J,K+1)+U(I+1,J,K+1)))*(1.-FIXKU(K))   &
          +0.5*(U(I,J,N3)+U(I+1,J,N3))*FIXKU(K)
       UM=VVDZ(K)*0.25                                        &
          *(SDZ(K)*(U(I,J,K-1)+U(I+1,J,K-1))                  &
           +SDZ(K-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXKL(K))     &
          +0.5*(U(I,J,0)+U(I+1,J,0))*FIXKL(K)
       VG13=SSDZ(K)*(UP-UM)

       VP=VVDX(I+1)*0.25                                      &
          *(SDX(I+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +SDX(I)*(V(I+1,J,K)+V(I+1,J+1,K)))*(1.-FIXIU(I))   &
          +0.5*(V(N1,J,K)+V(N1,J+1,K))*FIXIU(I)
       VM=VVDX(I)*0.25                                        &
          *(SDX(I)*(V(I-1,J,K)+V(I-1,J+1,K))                  &
           +SDX(I-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXIL(I))     &
          +0.5*(V(0,J,K)+V(0,J+1,K))*FIXIL(I)
       VG21=SSDX(I)*(VP-VM)

       VG22=SSDY(J)*(V(I,J+1,K)-V(I,J,K))

       VP=VVDZ(K+1)*0.25                                      &
          *(SDZ(K+1)*(V(I,J,K)+V(I,J+1,K))                    &
           +SDZ(K)*(V(I,J,K+1)+V(I,J+1,K+1)))*(1.-FIXKU(K))   &
          +0.5*(V(I,J,N3)+V(I,J+1,N3))*FIXKU(K)
       VM=VVDZ(K)*0.25                                        &
          *(SDZ(K)*(V(I,J,K-1)+V(I,J+1,K-1))                  &
           +SDZ(K-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXKL(K))     &
          +0.5*(V(I,J,0)+V(I,J+1,0))*FIXKL(K)
       VG23=SSDZ(K)*(VP-VM)

       WP=VVDX(I+1)*0.25                                      &
          *(SDX(I+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +SDX(I)*(W(I+1,J,K)+W(I+1,J,K+1)))*(1.-FIXIU(I))   &
          +0.5*(W(N1,J,K)+W(N1,J,K+1))*FIXIU(I)
       WM=VVDX(I)*0.25                                        &
          *(SDX(I)*(W(I-1,J,K)+W(I-1,J,K+1))                  &
           +SDX(I-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXIL(I))     &
          +0.5*(W(0,J,K)+W(0,J,K+1))*FIXIL(I)
       VG31=SSDX(I)*(WP-WM)

       WP=VVDY(J+1)*0.25                                      &
          *(SDY(J+1)*(W(I,J,K)+W(I,J,K+1))                    &
           +SDY(J)*(W(I,J+1,K)+W(I,J+1,K+1)))*(1.-FIXJU(J))   &
          +0.5*(W(I,N2,K)+W(I,N2,K+1))*FIXJU(J)
       WM=VVDY(J)*0.25                                        &
          *(SDY(J)*(W(I,J-1,K)+W(I,J-1,K+1))                  &
           +SDY(J-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXJL(J))     &
          +0.5*(W(I,0,K)+W(I,0,K+1))*FIXJL(J)
       VG32=SSDY(J)*(WP-WM)

       VG33=SSDZ(K)*(W(I,J,K+1)-W(I,J,K))

       SR(1)=VG11
       SR(2)=0.5*(VG12+VG21)
       SR(3)=0.5*(VG13+VG31)
       SR(4)=VG22
       SR(5)=0.5*(VG23+VG32)
       SR(6)=VG33

       WX=VG32-VG23
       WY=VG13-VG31
       WZ=VG21-VG12

       SIJAVGZ(I,J,1)=SIJAVGZ(I,J,1)+SR(1)*DT
       SIJAVGZ(I,J,2)=SIJAVGZ(I,J,2)+SR(2)*DT
       SIJAVGZ(I,J,3)=SIJAVGZ(I,J,3)+SR(3)*DT
       SIJAVGZ(I,J,4)=SIJAVGZ(I,J,4)+SR(4)*DT
       SIJAVGZ(I,J,5)=SIJAVGZ(I,J,5)+SR(5)*DT
       SIJAVGZ(I,J,6)=SIJAVGZ(I,J,6)+SR(6)*DT

       VORAVGZ(I,J,1)=VORAVGZ(I,J,1)+WX*DT
       VORAVGZ(I,J,2)=VORAVGZ(I,J,2)+WY*DT
       VORAVGZ(I,J,3)=VORAVGZ(I,J,3)+WZ*DT

       WRMSZ(I,J,1)=WRMSZ(I,J,1)+WX*WX*DT
       WRMSZ(I,J,2)=WRMSZ(I,J,2)+WX*WY*DT
       WRMSZ(I,J,3)=WRMSZ(I,J,3)+WX*WZ*DT
       WRMSZ(I,J,4)=WRMSZ(I,J,4)+WY*WY*DT
       WRMSZ(I,J,5)=WRMSZ(I,J,5)+WY*WZ*DT
       WRMSZ(I,J,6)=WRMSZ(I,J,6)+WZ*WZ*DT

       SRSR= 2.*SR(2)**2.+SR(1)**2.+2.*SR(3)**2.+SR(4)**2.+2.*SR(5)**2.+SR(6)**2.
       SSAVGZ(I,J)=SSAVGZ(I,J)+SRSR*DT

      ENDDO
      ENDDO
      ENDDO

      IF (MOD(NTIME,NPRIAVG).EQ.0) THEN

       TIMEENDZ=TIME
       IHISTENDZ=IHIST

       tfn1='favZ'
       idg1=ihistinitz/100000
       idg2=(ihistinitz-idg1*100000)/10000
       idg3=(ihistinitz-idg1*100000-idg2*10000)/1000
       idg4=(ihistinitz-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihistinitz-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihistinitz-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn2=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       idg1=ihist/100000
       idg2=(ihist-idg1*100000)/10000
       idg3=(ihist-idg1*100000-idg2*10000)/1000
       idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
       idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
       idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
       tfn3=char(idg1+48)//char(idg2+48)//char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
       tname=tfn1//tfn2//tfnh//tfn3

       OPEN(NAV,FILE=tname,FORM='UNFORMATTED')
       WRITE(NAV)N1M,N2M,N3M,RE
       WRITE(NAV)TIMEINITZ,TIMEENDZ,IHISTINITZ,IHISTENDZ
       WRITE(NAV)(( UAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
       WRITE(NAV)(( VAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
       WRITE(NAV)(( WAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
       WRITE(NAV)(((UIUJAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,6)
       WRITE(NAV)(( PAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
       WRITE(NAV)(( P2AVGZ(I,J)    ,I=1,N1M),J=1,N2M)
       WRITE(NAV)(((SIJAVGZ(I,J,L) ,I=1,N1M),J=1,N2M),L=1,6)
       WRITE(NAV)(((VORAVGZ(I,J,L) ,I=1,N1M),J=1,N2M),L=1,3)
       WRITE(NAV)(((WRMSZ(I,J,L)   ,I=1,N1M),J=1,N2M),L=1,6)
       WRITE(NAV)(( SSAVGZ(I,J)    ,I=1,N1M),J=1,N2M)
       CLOSE(NAV)

       NAV       = NAV+1
       UAVGZ    = 0.
       VAVGZ    = 0.
       WAVGZ    = 0.
       UIUJAVGZ = 0.
       PAVGZ    = 0.
       P2AVGZ   = 0.
       SIJAVGZ  = 0.
       VORAVGZ  = 0.
       SSAVGZ   = 0.
       WRMSZ    = 0.
       TIMEINITZ = TIME
       IHISTINITZ= IHIST

      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE PREFLD
!=======================================================================
!
!     Initial field is read from previous flow field
!     (fileprevel in lica.in, read in subroutine READINPUTFILE).
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8 :: NN1,NN2,NN3,IPSS,IIPX,IIPZ,IDUM
      REAL*8    :: RRE,PRAA,DUM

      U = 0.
      V = 0.
      W = 0.
      P = 0.

!     dum for future use
      OPEN(12,FILE=fileprevel,FORM='UNFORMATTED')
      READ(12) NN1,NN2,NN3,RRE
      READ(12) IHIST,M,TIME,DT
      READ(12) IPSS,IIPX,idum,IIPZ,idum,idum,idum,idum,idum,idum
      READ(12) PRAA, dum, dum, dum, dum, dum, dum, dum, dum, dum
      READ(12) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      READ(12) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      READ(12) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      READ(12) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      CLOSE(12)

      IF (IPZERO .EQ. 1) P = 0.

      WRITE(*,*)' '
      WRITE(*,100)
      WRITE(*,101)
      WRITE(*,102) fileprevel
      WRITE(*,103) RRE
      WRITE(*,104) NN1,NN2,NN3
      WRITE(*,105) IHIST,M,TIME,DT
      WRITE(*,106) IPSS,PRAA

  100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
  101 FORMAT('INITIAL FIELD      : READING DONE')
  102 FORMAT('INITIAL FIELD NAME : ',A30)
  103 FORMAT('RE=',ES12.3)
  104 FORMAT('N1=',I12,'  N2=',I12,'  N3=',I12)
  105 FORMAT('IHIST=',I9,'  M=',I9,'  TIME=',F10.5,'  DT=',F12.8)
  106 FORMAT('IPS=',I11,'  PRA=',F11.3)

!     Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
       DO J=0,N2
       DO I=1,N1
          U(I,J,0) =U(I,J,N3M)
          U(I,J,N3)=U(I,J,1)
       ENDDO
       ENDDO

       DO J=1,N2
       DO I=0,N1
          V(I,J,0) =V(I,J,N3M)
          V(I,J,N3)=V(I,J,1)
       ENDDO
       ENDDO

       DO J=0,N2
       DO I=0,N1
          W(I,J,0) =W(I,J,N3M)
          W(I,J,N3)=W(I,J,1)
       ENDDO
       ENDDO
      ENDIF

!     X PERIODICITY
      IF (IPX .EQ. 1) THEN
       DO K=0,N3
       DO J=0,N2
          U(0 ,J,K)=U(N1M,J,K)
          U(N1,J,K)=U(1  ,J,K)
       ENDDO
       ENDDO

       DO K=0,N3
       DO J=1,N2
          V(0 ,J,K)=V(N1M,J,K)
          V(N1,J,K)=V(1  ,J,K)
       ENDDO
       ENDDO

       DO K=1,N3
       DO J=0,N2
          W(0 ,J,K)=W(N1M,J,K)
          W(N1,J,K)=W(1  ,J,K)
       ENDDO
       ENDDO
      ENDIF

!     Z PERIODICITY
      IF (IPZ .EQ. 1) THEN
       DO J=1,N2M
       DO I=1,N1M
          P(I,J,0) =P(I,J,N3M)
          P(I,J,N3)=P(I,J,1)
       ENDDO
       ENDDO
      ENDIF

!     X PERIODICITY
      IF (IPX .EQ. 1) THEN
       DO K=1,N3M
       DO J=1,N2M
          P(0 ,J,K)=P(N1M,J,K)
          P(N1,J,K)=P(1  ,J,K)
       ENDDO
       ENDDO
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE WRITEFIELD
!=======================================================================
!
!     Make an output file of instantaneous field, when MOD(NTIME,NPRINT)=0,
!     NPRINT : instantaneous field file printing interval
!     tfn1='fld': prefix for instantaneous flow field
!     tname     : instantaneous field file name ex) fld006100
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,V,W,P
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8     :: IDUM,idg1,idg2,idg3,idg4,idg5,idg6
      REAL*8        :: DUM
      CHARACTER*9   :: tname
      CHARACTER*3   :: tfn1

      idum = 0
      dum  = 0.

      tfn1='fld'
      idg1=ihist/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tname=tfn1//char(idg1+48)//char(idg2+48)//                       &
           char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)

      OPEN(NV,FILE=tname,FORM='UNFORMATTED')
      WRITE(NV) N1,N2,N3,RE
      WRITE(NV) IHIST,M,TIME,DT
      WRITE(NV) IPS,IPX,idum,IPZ,idum,idum,idum,idum,idum,idum
      WRITE(NV) PRA, dum, dum, dum, dum, dum, dum, dum, dum, dum
      WRITE(NV) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      WRITE(NV) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      WRITE(NV) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      WRITE(NV) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
      CLOSE(NV)

      NV=NV+1

      RETURN
      END

!=======================================================================
       SUBROUTINE MEANPGR
!=======================================================================
!
!     CALCULATE MEAN PRESSURE GRADIENT FOR CONST. MASS FLOW RATE IN CHANNEL FLOW
!
!     The mean pressure gradient (dP/dx) required to keep the mass flow rate constant
!     is determined by integrating the wall-shear stresses at the channel walls.
!
!     When there is an IBM body inside the channel, the IBM forcing has to be included
!     -> FCV(N,1) is considered (forcing in the streamwise direction)
!
!     PMIL: wall-shear stress at the lower-channel wall
!     PMIU: wall-shear stress at the upper-channel wall
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U,FCV,IFC,JFC,KFC
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      INTEGER*8     :: N
      REAL*8        :: UGU,UGL,TMP1,TMP2,TMP3

         PMI  = 0.
         PMIU = 0.
         PMIL = 0.

!$OMP PARALLEL DO       &
!$OMP private(UGU,UGL)  &
!$OMP reduction(+:PMI)  &
!$OMP reduction(+:PMIU) &
!$OMP reduction(+:PMIL)
      DO K=1,N3M
      DO I=1,N1M
         UGU  = (U(I,N2,K)-U(I,N2M,K))*VVDY(N2)*VDX(I)*SDZ(K)
         UGL  = (U(I,1,K)-U(I,0,K))*VVDY(1)*VDX(I)*SDZ(K)
         PMI  = PMI  + (UGU-UGL)*REI
         PMIU = PMIU +  UGU*REI
         PMIL = PMIL -  UGL*REI
      ENDDO
      ENDDO

         PMI=PMI/(XL*YL*ZL)
         PMIL=PMIL/(XL*YL*ZL)
         PMIU=PMIU/(XL*YL*ZL)

      IF (IDPDX.EQ.1) THEN
         PMI=-1.
      ENDIF


      IF (IBMON.NE.0) THEN

         TMP1 = 0.
         TMP2 = 0.
         TMP3 = 0.

!$OMP PARALLEL DO        &
!$OMP reduction(-:TMP1)  &
!$OMP reduction(-:TMP2)  &
!$OMP reduction(-:TMP3)
        DO N=1,NFC_INTP(1)+NFC_INNER(1)
          TMP1=TMP1+FCV(N,1)*VDX(IFC(N,1))*SDY(JFC(N,1))*SDZ(KFC(N,1))
          IF(YMP(JFC(N,1)) .LE. 0) THEN
           TMP2=TMP2+FCV(N,1)*VDX(IFC(N,1))*SDY(JFC(N,1))*SDZ(KFC(N,1))
          ELSE
           TMP3=TMP3+FCV(N,1)*VDX(IFC(N,1))*SDY(JFC(N,1))*SDZ(KFC(N,1))
          ENDIF
        ENDDO

         PMI =PMI +TMP1/(XL*YL*ZL)
         PMIL=PMIL+TMP2/(XL*YL*ZL)
         PMIU=PMIU+TMP3/(XL*YL*ZL)

      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE QCALC
!=======================================================================
!
!     QREF  : 5600./4200. REb(=Ub*2H/nu)=5600, REc(=Uc(lam)*H/nu)=4200.
!     QFLUX : check flow rate. QREF is multiplied to check the flow rate is constant
!             corresponding to Rec=4200
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8        :: PI,QREF

      QFLUX = 0.
      QREF  = 4./3.

!$OMP PARALLEL DO &
!$OMP reduction(+:QFLUX)
      DO K=1,N3M
      DO J=1,N2M
      QFLUX = QFLUX+U(1,J,K)*SDY(J)*SDZ(K)
      ENDDO
      ENDDO
      QFLUX=QFLUX/ZL/QREF

      RETURN
      END

!=======================================================================
      SUBROUTINE QVOLCALC(QM)
!=======================================================================
!
!     Calculate the volume flow rate at certain steps to keep the flow rate constant
!     QVOL : volume flow rate
!     QVOL1: volume flow rate of previous step
!     QVOL2: volume flow rate of present step
!     -> PHCAP=QVOL1-QVOL2 : correction term for constant volume flow rate
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : U
      IMPLICIT NONE
      INTEGER*8     :: I,J,K
      REAL*8      :: QM

      QM = 0.

!$OMP PARALLEL DO &
!$OMP reduction(+:QM)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        QM = QM+U(I,J,K)*VDX(I)*SDY(J)*SDZ(K)
      ENDDO
      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE LAGFORCE
!=======================================================================
!
!     Material derivate of velocity, one of major contribution of the
!       force on a body
!
!     Reference
!       Lee et al., 2011, Sources of spurious force oscillations
!       from an immersed boundary method for moving-body problems,
!       J. Comp. Phys., 230, 2677-2695.
!
!     Variables
!       DUDTA,DVDTA,DWDTA: Volume integration of material derivative of
!                          velocity over a body
!       DUDTR: Time derivative of velocity
!       RK3XO,RK3XOO: Convective terms at each RK3 sub-step
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY :  RK3XO,RK3YO,RK3ZO,DUDTR  &
                                ,RK3XOO,RK3YOO,RK3ZOO,IFC,JFC,KFC
      IMPLICIT NONE
      INTEGER*8     :: I,J,K,N,L,II,JJ,KK

      IF(MSUB .EQ. 1) THEN
        DUDTA = 0.
        DVDTA = 0.
        DWDTA = 0.
      ENDIF

!$OMP PARALLEL DO  &
!$OMP private(II,JJ,KK)
      DO L=1,3
      DO N=1,NFC_INTP(L)+NFC_INNER(L)
        II=IFC(N,L)
        JJ=JFC(N,L)
        KK=KFC(N,L)
        !Intermediate information
        IF( L .EQ. 1) THEN
          DUDTA=DUDTA+VDX(II)*SDY(JJ)*SDZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3XO(II,JJ,KK)+RO(MSUB)*RK3XOO(II,JJ,KK)))
        ELSEIF( L .EQ. 2) THEN
          DVDTA=DVDTA+SDX(II)*VDY(JJ)*SDZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3YO(II,JJ,KK)+RO(MSUB)*RK3YOO(II,JJ,KK)))
        ELSEIF( L .EQ. 3) THEN
          DWDTA=DWDTA+SDX(II)*SDY(JJ)*VDZ(KK)*(DUDTR(N,L)-              &
                (GAMMA(MSUB)*RK3ZO(II,JJ,KK)+RO(MSUB)*RK3ZOO(II,JJ,KK)))
        ENDIF

      ENDDO
      ENDDO

      RETURN
      END

!=======================================================================
      SUBROUTINE DRAGLIFT
!=======================================================================
!
!     Force coefficients on a body are traced when IBMON = 1.
!
!     Variables
!       CDAREA: Frontral area for drag/lift coefficients
!               It is defined in GETCDAREA subroutine in lica_cylinder.f90
!
!     Output
!       5001->FILE='ftr_cdcl.dat' only intergration of momentum forcing
!       5010->FILE='ftrcdcl2.dat' integration of both mometum forcing and
!                                 material derivative of velocity
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : IFC,JFC,KFC,FCVAVG
      IMPLICIT NONE
      INTEGER*8     :: N,L
      REAL*8        :: TMP1,TMP2,TMP3

      FORCESUM  = 0.
      FORCESUMA = 0.
      TMP1 = 0.
      TMP2 = 0.
      TMP3 = 0.

!-----sum forcing values and compute Cd and Cl
!-----sphere
!$OMP PARALLEL DO &
!$OMP reduction(-:TMP1)
      DO 40 N=1,NFC_INTP(1)+NFC_INNER(1)
      TMP1=TMP1-FCVAVG(N,1)*VDX(IFC(N,1))*SDY(JFC(N,1))*SDZ(KFC(N,1))
 40   CONTINUE

!$OMP PARALLEL DO &
!$OMP reduction(-:TMP2)
      DO 42 N=1,NFC_INTP(2)+NFC_INNER(2)
      TMP2=TMP2-FCVAVG(N,2)*SDX(IFC(N,2))*VDY(JFC(N,2))*SDZ(KFC(N,2))
 42   CONTINUE

!$OMP PARALLEL DO &
!$OMP reduction(-:TMP3)
      DO 44 N=1,NFC_INTP(3)+NFC_INNER(3)
      TMP3=TMP3-FCVAVG(N,3)*SDX(IFC(N,3))*SDY(JFC(N,3))*VDZ(KFC(N,3))
 44   CONTINUE

      FORCESUM(1)=TMP1
      FORCESUM(2)=TMP2
      FORCESUM(3)=TMP3

      FORCESUMA(1)=TMP1+DUDTA
      FORCESUMA(2)=TMP2+DVDTA
      FORCESUMA(3)=TMP3+DWDTA

      DO L= 1,3
       FORCESUM(L)=FORCESUM(L)*2./CDAREA
       FORCESUMA(L)=FORCESUMA(L)*2./CDAREA
      ENDDO

      WRITE(5001,101) TIME,(FORCESUM(L),L=1,3)
      WRITE(5010,101) TIME,(FORCESUMA(L),L=1,3)
  101 FORMAT(F13.5,3ES15.7)

      RETURN
      END

!=======================================================================
      SUBROUTINE WRITEHISTORY
!=======================================================================
!
!     Trace the maximum CFL#, divergenceU, and QMASS(IBM).
!     The tracing interval is set in lica.in (NPIN)
!
!     *In case of channel flow,
!     - the flow rate and the mean pressure gradient is traced
!     - QFLUX should be constant, while PMI, PMIL, PMIU fluctuate
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE

      IF (ICH.EQ.1) CALL QCALC

      WRITE(5000,130) TIME,DT,CFLMAX,DVMAX,QMMAX
      FLUSH(5000)
 130  FORMAT(F13.5,4ES15.7)

      IF (ICH.EQ.1) THEN
      WRITE(5007,140) TIME,QFLUX,-PMIU,-PMIL,-PMI
 140  FORMAT(F13.5,4ES20.12)
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE WRITEFTRTIME
!=======================================================================
!
!     Real time measurements for individual processes
!     Printed at every computational step
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      REAL*8     ::   FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIME5

      FTRTIME1=(SGSTIME_E(3)-SGSTIME_B(3)+SGSTIME_E(2)-SGSTIME_B(2)+SGSTIME_E(1)-SGSTIME_B(1))
      FTRTIME2=(RHSTIME_E(3)-RHSTIME_B(3)+RHSTIME_E(2)-RHSTIME_B(2)+RHSTIME_E(1)-RHSTIME_B(1))
      FTRTIME3=(LHSTIME_E(3)-LHSTIME_B(3)+LHSTIME_E(2)-LHSTIME_B(2)+LHSTIME_E(1)-LHSTIME_B(1))
      FTRTIME4=(POISSTIME_E(3)-POISSTIME_B(3)+POISSTIME_E(2)-POISSTIME_B(2)+POISSTIME_E(1)-POISSTIME_B(1))
      FTRTIME5=TIME_END-TIME_BEGIN

      WRITE(*,201) FTRTIME1
      WRITE(*,202) FTRTIME2
      WRITE(*,203) FTRTIME3
      WRITE(*,204) FTRTIME4
      WRITE(*,205) FTRTIME5
 201  FORMAT('TIME FOR SGS    : ',F12.3,' SECONDS')
 202  FORMAT('TIME FOR RHS    : ',F12.3,' SECONDS')
 203  FORMAT('TIME FOR LHS    : ',F12.3,' SECONDS')
 204  FORMAT('TIME FOR POISSON: ',F12.3,' SECONDS')
 205  FORMAT('TIME OF OPERTION: ',F12.3,' SECONDS')
      WRITE(5009,206)TIME,FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIME5
 206  FORMAT(F13.5,5F12.4)

      RETURN
      END

!=======================================================================
      SUBROUTINE FTRFILES(IOPT)
!=======================================================================
!
!     IOPT=0, open the files for recording the time-traced variables.
!     When the initial field is read from a previous flow-field (IREAD=1)
!     tracing is continued on from the existing files.
!
!     IOPT=1 when the calculation is finished (M=NTST),
!     and the time-tracing files are closed.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      IMPLICIT NONE
      INTEGER*8    :: IOPT

      IF (IOPT.EQ.0) THEN
       IF (IREAD.EQ.0) THEN
        OPEN(5000,FILE='ftrhist.dat')
        IF (IBMON .NE. 0) OPEN(5001,FILE='ftrcdcl.dat')
        IF (ILES  .EQ. 1) OPEN(5002,FILE='ftrles.dat')
        IF (NTEST .NE. 0) OPEN(5003,FILE='ftru.dat')
        IF (NTEST .NE. 0) OPEN(5004,FILE='ftrv.dat')
        IF (NTEST .NE. 0) OPEN(5005,FILE='ftrw.dat')
        IF (NTEST .NE. 0) OPEN(5006,FILE='ftrp.dat')
        IF (ICH   .EQ. 1) OPEN(5007,FILE='ftrch.dat')
        OPEN(5009,FILE='ftrtime.dat')
        IF (IBMON .NE. 0) OPEN(5010,FILE='ftrcdcl2.dat')
       ELSE IF (IREAD.EQ.1) THEN
        OPEN(5000,FILE='ftrhist.dat',POSITION='APPEND')
        IF (IBMON .NE. 0) OPEN(5001,FILE='ftrcdcl.dat',POSITION='APPEND')
        IF (ILES  .EQ. 1) OPEN(5002,FILE='ftrles.dat',POSITION='APPEND')
        IF (NTEST .NE. 0) OPEN(5003,FILE='ftru.dat',POSITION='APPEND')
        IF (NTEST .NE. 0) OPEN(5004,FILE='ftrv.dat',POSITION='APPEND')
        IF (NTEST .NE. 0) OPEN(5005,FILE='ftrw.dat',POSITION='APPEND')
        IF (NTEST .NE. 0) OPEN(5006,FILE='ftrp.dat',POSITION='APPEND')
        IF (ICH   .EQ. 1) OPEN(5007,FILE='ftrch.dat',POSITION='APPEND')
        OPEN(5009,FILE='ftrtime.dat',POSITION='APPEND')
        IF (IBMON .NE. 0) OPEN(5010,FILE='ftrcdcl2.dat',POSITION='APPEND')
       ENDIF


      ELSE
         CLOSE(5000)
         IF (IBMON .NE. 0) CLOSE(5001)
         IF (ILES  .EQ. 1) CLOSE(5002)
         IF (NTEST .NE. 0) CLOSE(5003)
         IF (NTEST .NE. 0) CLOSE(5004)
         IF (NTEST .NE. 0) CLOSE(5005)
         IF (NTEST .NE. 0) CLOSE(5006)
         IF (ICH   .EQ. 1) CLOSE(5007)
         CLOSE(5009)
         IF (IBMON .NE. 0) CLOSE(5010)
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE ALLO_ARRAY
!=======================================================================
!
!     Allocate flow array
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE

      ALLOCATE( U(0:M1,0:M2,0:M3))
      ALLOCATE( V(0:M1,0:M2,0:M3))
      ALLOCATE( W(0:M1,0:M2,0:M3))
      ALLOCATE( P(0:M1,0:M2,0:M3))
      ALLOCATE( RHS1(0:M1,0:M2,0:M3,3))
      ALLOCATE( TNU (0:M1,0:M2,0:M3))
      ALLOCATE( TNUX(0:M1,0:M2,0:M3))
      ALLOCATE( TNUY(0:M1,0:M2,0:M3))
      ALLOCATE( TNUZ(0:M1,0:M2,0:M3))
      ALLOCATE( RK3XO(M1M,M2M,M3M))
      ALLOCATE( RK3YO(M1M,M2M,M3M))
      ALLOCATE( RK3ZO(M1M,M2M,M3M))
      ALLOCATE(RK3XOO(M1M,M2M,M3M))
      ALLOCATE(RK3YOO(M1M,M2M,M3M))
      ALLOCATE(RK3ZOO(M1M,M2M,M3M))

      IF (IBMON.EQ.1) THEN
        ALLOCATE(IFC(MBODY,3),JFC(MBODY,3),KFC(MBODY,3))
        ALLOCATE(FCV(MBODY,3),FCVAVG(MBODY,3))
        ALLOCATE(MPI(MBODY0,3,3))
        ALLOCATE(GFI(MBODY0,3,0:2,0:2,0:2))
        ALLOCATE(QMASS(M1M,M2M,M3M))
        ALLOCATE(DUDTR(MBODY,3))
      ENDIF

      IF (IMOVINGON.EQ.1) THEN
        ALLOCATE(JTYPE(MBODY0,3))
        ALLOCATE(INOUT(0:M1,0:M2,0:M3,3))
        ALLOCATE(ROTOR_PHI_CENTER(0:M1,0:M2,0:M3))
      ENDIF

      IF (ILES .EQ. 1) THEN
        ALLOCATE(INZ(MBODY),JNZ(MBODY),KNZ(MBODY))
        ALLOCATE(NWALL_DVM(M1M,M2M,M3M))
        ALLOCATE(CFX1 (M1M,-1:1))
        ALLOCATE(CFX2 (M1M,-2:2))
        ALLOCATE(CFZP1(M3M,-1:1))
        ALLOCATE(CFZP2(M3M,-2:2))
      ENDIF

      IF (IAVG .EQ. 1) THEN
        ALLOCATE(UAVG(M1M,M2M,M3M))
        ALLOCATE(VAVG(M1M,M2M,M3M))
        ALLOCATE(WAVG(M1M,M2M,M3M))
        ALLOCATE(PAVG(M1M,M2M,M3M))
        ALLOCATE(P2AVG(M1M,M2M,M3M))
        ALLOCATE(SSAVG(M1M,M2M,M3M))
        ALLOCATE(UIUJAVG(M1M,M2M,M3M,6))
        ALLOCATE(VORAVG(M1M,M2M,M3M,3))
        ALLOCATE(WRMS(M1M,M2M,M3M,6))
      ELSE IF (IAVG .EQ. 2) THEN
        ALLOCATE(UAVGZ(M1M,M2M))
        ALLOCATE(VAVGZ(M1M,M2M))
        ALLOCATE(WAVGZ(M1M,M2M))
        ALLOCATE(PAVGZ(M1M,M2M))
        ALLOCATE(P2AVGZ(M1M,M2M))
        ALLOCATE(SSAVGZ(M1M,M2M))
        ALLOCATE(VORAVGZ(M1M,M2M,3))
        ALLOCATE(WRMSZ(M1M,M2M,6))
        ALLOCATE(UIUJAVGZ(M1M,M2M,6))
        ALLOCATE(SIJAVGZ(M1M,M2M,6))
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE DEALLO_ARRAY
!=======================================================================
!
!     Deallocate flow array
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      IMPLICIT NONE

      DEALLOCATE( U,V,W,P)
      DEALLOCATE( RHS1)
      DEALLOCATE(TNU ,TNUX,TNUY,TNUZ)
      DEALLOCATE(RK3XO,RK3YO,RK3ZO)
      DEALLOCATE(RK3XOO,RK3YOO,RK3ZOO)

      IF (IBMON.EQ.1) THEN
        DEALLOCATE(IFC,JFC,KFC)
        DEALLOCATE(FCV,FCVAVG)
        DEALLOCATE(DUDTR)
        DEALLOCATE(MPI)
        DEALLOCATE(GFI)
        DEALLOCATE(QMASS)
      ENDIF

      IF (IMOVINGON.EQ.1) THEN
        DEALLOCATE(JTYPE)
        DEALLOCATE(INOUT)
        DEALLOCATE(ROTOR_PHI_CENTER)
      ENDIF

      IF (ILES .EQ. 1) THEN
        DEALLOCATE(INZ,JNZ,KNZ)
        DEALLOCATE(NWALL_DVM)
        DEALLOCATE(CFX1,CFX2,CFZP1,CFZP2)
      ENDIF

      IF (IAVG .EQ. 1) THEN
        DEALLOCATE(UAVG,VAVG,WAVG,PAVG)
        DEALLOCATE(P2AVG,SSAVG)
        DEALLOCATE(UIUJAVG,VORAVG,WRMS)
      ELSE IF (IAVG .EQ. 2) THEN
        DEALLOCATE(UAVGZ,VAVGZ,WAVGZ,PAVGZ)
        DEALLOCATE(P2AVGZ,SSAVGZ)
        DEALLOCATE(VORAVGZ,WRMSZ,UIUJAVGZ,SIJAVGZ)
      ENDIF

      RETURN
      END

!=======================================================================
! Export actual native staggered point classifications, not a voxel surrogate.
subroutine export_ibm_grid_slice
 use mod_common
 use mod_flowarray,only:INOUT,IFC,JFC,KFC,JTYPE
 use mod_rotor_motion,only:component_functions,rotating_owner,length_mm,axis_x_mm
 implicit none
 integer::l,i,j,k,n,unit,planeid
 integer::forcing(0:M2,0:M3)
 real*8::xx,yy,zz,fs,fr,fh,dmin
 open(newunit=unit,file='native_ibm_grid_slice.dat',status='replace')
 write(unit,'(A)') '# component i j k x_mm y_mm z_mm inout interpolation_type rotating_owner'
 do l=1,3
  dmin=huge(1D0);planeid=1
  do i=2,N1M-1
   xx=XMP(i)
   if(l==1)xx=X(i)
   if(abs(xx*length_mm+axis_x_mm+62D0)<dmin)then
    dmin=abs(xx*length_mm+axis_x_mm+62D0);planeid=i
   endif
  enddo
  forcing=0
  do n=1,NFC_INTP(l)
   if(IFC(n,l)==planeid)forcing(JFC(n,l),KFC(n,l))=JTYPE(n,l)
  enddo
  i=planeid;xx=XMP(i)
  if(l==1)xx=X(i)
  do k=2,N3M-1
   zz=ZMP(k)
   if(l==3)zz=Z(k)
   if(zz*length_mm<9.3D0.or.zz*length_mm>11.5D0)cycle
   do j=2,N2M-1
    yy=YMP(j)
    if(l==2)yy=Y(j)
    if(abs(yy*length_mm)>3.5D0)cycle
    write(unit,'(4I8,3ES24.15,3I8)')l,i,j,k,xx*length_mm+axis_x_mm,yy*length_mm,zz*length_mm, &
      INOUT(i,j,k,l),forcing(j,k),merge(1,0,rotating_owner(xx,yy,zz))
   enddo
  enddo
 enddo
 close(unit)
end subroutine
