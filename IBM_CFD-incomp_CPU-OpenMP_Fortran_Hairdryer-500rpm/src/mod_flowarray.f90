MODULE MOD_FLOWARRAY

!======================== BASIC VARIABLES ==============================
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: U,V,W,P
      REAL*8, ALLOCATABLE :: ROTOR_PHI_CENTER(:,:,:)
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: RHS1
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: RK3XO,RK3YO,RK3ZO
!     U,V,W,P           : Velocity and Pressure
!     RHS1              : Explicitly treated terms in Navier-Stokes equation
!     RK3XO,RK3YO,RK3ZO : Previous step information for Runge-Kutta method (RK3)

!======================== IBM VARIABLES ================================
      INTEGER*8, DIMENSION (:,:),    ALLOCATABLE :: IFC,JFC,KFC
      INTEGER*8, DIMENSION (:,:,:),  ALLOCATABLE :: MPI
      INTEGER*8, DIMENSION (:,:,:,:),ALLOCATABLE :: INOUT
      INTEGER*8, DIMENSION (:,:),    ALLOCATABLE :: JTYPE
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: FCV,FCVAVG
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: DUDTR
      REAL*8, DIMENSION (:,:,:,:,:), ALLOCATABLE :: GFI
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: QMASS
!     IFC,JFC,KFC : Index of interpolation point for IBM
!     MPI         : Distinguish the solid-to-fluid direction
!     INOUT       : 1 => Outside point, 0 => Inside point
!     JTYPE       : The type of interpolation (1: 1D, 2: 2D, 3: 3D)
!     FCV,FCVAVG  : Instantaneous and average forces obtained from the momentum forcing in IBM
!     DUDTR       : Inertia contribution in IBM forcing
!     GFI         : Geoetric factor for IBM interpolation
!     QMASS       : Mass source and sinke

!======================== LES VARIABLES ================================
      INTEGER*8, DIMENSION (:),      ALLOCATABLE :: INZ,JNZ,KNZ
      INTEGER*8, DIMENSION (:,:,:),  ALLOCATABLE :: NWALL_DVM
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: TNU,TNUX,TNUY,TNUZ
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: CFX1,CFX2,CFZP1,CFZP2
!     INZ,JNZ,KNZ           : Index for zero turbulent viscosity point
!     NWALL_DVM             : Zero turbulent viscosity in the IB body
!     TNU,TNUX,TNUY,TNUZ    : Turbulent viscosity at each velocity control volume
!     CFX1,CFX2,CFZP1,CFZP2 : Coefficient for LES filtering

!======================== SEMI_IMPLICIT ================================
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: RK3XOO,RK3YOO,RK3ZOO
!     RK3XOO,RK3YOO,RK3ZOO : Convection term from inertia contribution in IBM forcing

!====================== FILED AVG. VARIABLES ===========================
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: UAVG,VAVG,WAVG,PAVG
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: P2AVG,SSAVG
      REAL*8, DIMENSION (:,:,:,:),   ALLOCATABLE :: UIUJAVG,VORAVG,WRMS

      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: UAVGZ,VAVGZ,WAVGZ,PAVGZ
      REAL*8, DIMENSION (:,:),       ALLOCATABLE :: P2AVGZ,SSAVGZ
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: VORAVGZ,WRMSZ
      REAL*8, DIMENSION (:,:,:),     ALLOCATABLE :: UIUJAVGZ,SIJAVGZ
!     UAVG,VAVG,WAVG,PAVG : AVERAGED VELOCITY AND PRESSURE
!     P2AVG,SSAVG         : AVERAGED PRESSURE FLUCTUATIONS, AVERAGED MAGNITUDE OF STRAIN
!     UIUJAVG,VORAVG,WRMS : AVERAGED REYNOLDS STRESS, AVERAGED ROOT-MEAN-SQUARE OF VORTICITY
!
!     UAVGZ,VAVGZ,WAVGZ,PAVGZ : AVERAGED VELOCITY AND PRESSURE OVER Z-DIRECTION
!     P2AVGZ,SSAVGZ           : AVERAGED PRESSURE FLUCTUATIONS OVER Z-DIRECTION,
!                               AVERAGED MAGNITUDE OF STRAIN OVER Z-DIRECTION
!     VORAVGZ,WRMSZ           : AVERAGED VORTICITY OVER Z-DIRECTION,
!                             : AVERAGED ROOT-MEAN-SQUARE OF VORTICITY OVER Z-DIRECTION
!     UIUJAVGZ,SIJAVGZ        : AVERAGED REYNOLDS STRESS OVER Z-DIRECTION
!                             : AVERAGED MAGNITUDE OF STRAIN OVER Z-DIRECTION

END MODULE