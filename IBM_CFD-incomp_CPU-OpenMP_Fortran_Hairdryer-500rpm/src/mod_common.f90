MODULE MOD_COMMON

!======================= PARAMETERS ==================================
!-----GEOMETRY
      INTEGER*8, PARAMETER :: M1 = 105,M2 = 97,M3 = 97
      INTEGER*8, PARAMETER :: M1M=M1-1,M2M=M2-1,M3M=M3-1
      INTEGER*8, PARAMETER :: M3MH=M3M/2+1
!-----IMMERSED BOUNDARY
      INTEGER*8, PARAMETER :: MBODY0= 239616
      INTEGER*8, PARAMETER :: MBODY = 1246004	   
!-----HISTORY
      INTEGER*8, PARAMETER :: MHIST =  1000000    !# of time steps
!     M1,M2,M3    : Memory sizes for X,Y,Z direction grids
!                   (should be same to N1,N2,N3 in main code)
!     M1M,M2M,M3M : The number of cells for X,Y,Z directions
!     M3MH        : The number of wavenumber for Z direction for FFT to solve the Poisson equation
!     MBODY0,MBODY: Memory sizes for IBM interpolation and total (interpolation+inner) points
!     MHIST       : Maximum number of time steps

!========================= INPUTS ======================================
      REAL*8       :: RE
      INTEGER*8    :: IRESET,IREAD,IAVG,IPZERO
      INTEGER*8    :: NTST,NPRINT,NPRIAVG,NPIN
      INTEGER*8    :: IDTOPT
      INTEGER*8    :: NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR
      INTEGER*8    :: ILES,ISGSMOD,IDYNMON,IFILTER
      INTEGER*8    :: IPS
      INTEGER*8    :: IBMON,MASSON,IINTP,IMOVINGON,IFUNCorPTS
      INTEGER*8    :: NTEST,NTPRINT
      INTEGER*8    :: ITR(500),JTR(500),KTR(500)
      REAL*8       :: EPS_PTR,DTT,CFLFAC
      REAL*8       :: RESID1,WWSOR,PRA,CSGS,CSGSPS
      REAL*8       :: REI,PRAI
      CHARACTER*30 :: gridfile,fileprevel
!     Re                       : Reynolds number
!     IRESET,IREAD,IAVG,IPZERO : Reset the time, Read a previous field,
!                                Average option, Initialize pressure to zero
!     NTST,NPRINT,NPRIAVG,NPIN : Time steps for calculation,
!                                Period for writing the instantaneous field,
!                                Period for writing the average field,
!                                Period for writing the time history
!     IDTOPT                   : Time step size option

      !==Multigrid option==!
!     NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR : The number of coarsening level, W cycle option,
!                                        The number of iteration at the lowest level,
!                                        Maxumum iteration of Poisson equation,
!                                        Initial guess for Poisson equation iteration,
!                                        SOR acceleration factor

      !==Les option==!
!     ILES,ISGSMOD,IDYNMON,IFILTER     : LES option, SGS model option,
!                                        Dynamic SGS model option, Filtering direction

      !==Passive scalar option==!
!     IPS                        : Passive scalar option

      !==IBM option==!
!     IBMON,MASSON,IINTP         : IBM option, Masson option, Interpolation option
!     IMOVINGON :   MOVING IBM ON(1)/OFF(0)
!     IFUNCorPTS:   FUNCTION(1)/POINTS(2) FOR GEOMETRY

      !==Time trace option==!
!     NTEST,NTPRINT              : Time trace option, Period for time trace
!     ITR(500),JTR(500),KTR(500) : Point index for time trace

!     EPS_PTR,DTT,CFLFAC           : Initial perturbation,
!                                    Initial time step size, Maximum CFL number
!     RESID1,WWSOR,PRA,CSGS,CSGSPS : Poisson convergence criterion,
!                                    Over-relaxation parameter for SOR,
!                                    Prandtl number, SGS coefficient,
!                                    SGS coefficient for passive scalar
!     REI,PRAI                     : Reynolds number inverse, Prandtl number interse
!     gridfile,fileprevel          : grid file, field file readed when IREAD=1

!==================== GRID & MESH & DOMAIN =============================
      INTEGER*8    :: N1,N1M,N2,N2M,N3,N3M
      REAL*8       :: XL,YL,ZL
      REAL*8       :: X(M1),Y(M2),Z(M3)
      INTEGER*8    :: IPV(M1M),JPV(M2M),KPV(M3M),IMV(M1M),JMV(M2M),KMV(M3M)
      REAL*8       :: FIXIL(M1M),FIXIU(M1M),FIXJL(M2M),FIXJU(M2M),FIXKL(M3M),FIXKU(M3M)
      REAL*8       :: SDX(0:M1),SDY(0:M2),SDZ(0:M3)
      REAL*8       :: SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3)
      REAL*8       :: VDX(M1),VDY(M2),VDZ(M3)
      REAL*8       :: VVDX(M1),VVDY(M2),VVDZ(M3)
      REAL*8       :: XMP(0:M1),YMP(0:M2),ZMP(0:M3)
!     N1,N1M,N2,N2M,N3,N3M             : The number of grid points for X1,X2,X3 directions
!     XL,YL,ZL                         : Domain sizes
!     X(M1),Y(M2),Z(M3)                : Positions of grid points
!     IPV(M1M),JPV(M2M),KPV(M3M),IMV(M1M),JMV(M2M),KMV(M3M): Next and previous grid index
!     FIXIL(M1M),FIXIU(M1M),FIXJL(M2M),FIXJU(M2M),FIXKL(M3M),FIXKU(M3M) : Treat domain boundary
!     SDX(0:M1),SDY(0:M2),SDZ(0:M3)    : Grid sizes from cell face to face
!     SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3) : Inverse of grid sizes (SDX,SDY,SDZ)
!     VDX(M1),VDY(M2),VDZ(M3)          : Grid sizes from cell center to center
!     VVDX(M1),VVDY(M2),VVDZ(M3)       : Inverse of grid sizes (VDX,VDY,VDZ)
!     XMP(0:M1),YMP(0:M2),ZMP(0:M3)    : Position of cell center points

!================= BOUNDARY CONDITION VARIABLES ========================
      INTEGER*8    :: IPX,IPZ,I_BGPX,K_BGPZ
      INTEGER*8    :: JUT,JUB,JWT,JWB,KUT,KUB,KVT,KVB,JVT
      REAL*8       :: DUOUT(0:M2,0:M3),DVOUT(0:M2,0:M3),DWOUT(0:M2,0:M3)
      REAL*8       :: UOUT(M2M,M3M),VOUT(M2M,M3M),WOUT(M2M,M3M)
      REAL*8       :: VTOP(M1M,M3M),DVTOP(M1M,M3M)
!     IPX,IPZ,I_BGPX,K_BGPZ : Periodic boundary condition option for X and Z directions,
!                             Starting index for X and Z directions (depending on a boundary condition)
!     JUT,JUB,JWT,JWB,KUT,KUB,KVT,KVB,JVT ! Boundary condition option
!     DUOUT(0:M2,0:M3),DVOUT(0:M2,0:M3),DWOUT(0:M2,0:M3) : Convective boundary treatment
!                                                          at the intermediate step in fraction step method
!     UOUT(M2M,M3M),VOUT(M2M,M3M),WOUT(M2M,M3M)
!     VTOP(M1M,M3M),DVTOP(M1M,M3M)

!======================= IBM VARIABLES =================================
      INTEGER*8    :: NFC_INTP(3),NFC_INNER(3)
      INTEGER*8    :: IUVW(3),JUVW(3),KUVW(3)
!     NFC_INTP(3),NFC_INNER(3) : The number of interpolation and inner points for IBM
!     IUVW(3),JUVW(3),KUVW(3)  : Initial I,J,K index for searching forcing points

!======================= LES VARIABLES =================================
      INTEGER*8    :: NNZ
!     NNZ : Number of zero turbulent viscosity point

!======================= N-S EQUATION ==================================
      INTEGER*8    :: NTIME,IHIST
      INTEGER*8    :: M,MSUB
      INTEGER*8    :: NV,NAV
      REAL*8       :: DT,TIME,SUBDT
      REAL*8       :: DVMAX,CFLMAX,QMMAX
      REAL*8       :: ALPHA,GAMMA(3),RO(3)
      REAL*8       :: DTCONST,DTCONSTI,TEST1,ACOEF,ACOEFI
      REAL*8       :: AIU(M1),CIU(M1),AIVW(M1),CIVW(M1), &
                      AJV(M2),CJV(M2),AJUW(M2),CJUW(M2), &
                      AKW(M3),CKW(M3),AKUV(M3),CKUV(M3)
      REAL*8       :: TOTAL_TIME_B,TOTAL_TIME_E,                          &
                      TIME_BEGIN,TIME_END,SGSTIME_B(3),SGSTIME_E(3),      &
                      RHSTIME_B(3),RHSTIME_E(3),LHSTIME_B(3),LHSTIME_E(3),&
                      POISSTIME_B(3),POISSTIME_E(3)
!     NTIME,IHIST  : Time step index for current simulation, Total time step index
!     M,MSUB       : Time step index (same to NTIME), Runge-Kutta sub-time step index
!     NV,NAV       : Field adress
!     DT,TIME      : Time step size, Time
!     DVMAX,CFLMAX,QMMAX : Maximum velocity divergence, Maximum CFL number,
!                          Maximum mass source/sink for IBM
!     ALPHA,GAMMA(3),RO(3) : RK3 coefficients
!     DTCONST,DTCONSTI,TEST1,ACOEF,ACOEFI : DT for each RK3 step, inverse of DTCONST
!                                         : Poisson convergence criterion
!     AIU(M1),CIU(M1),AIVW(M1),CIVW(M1),  : Coefficient for ADI method to solve LHS
!     AJV(M2),CJV(M2),AJUW(M2),CJUW(M2),
!     AKW(M3),CKW(M3),AKUV(M3),CKUV(M3)
!     TOTAL_TIME_B,TOTAL_TIME_E,          : Measure the real time
!     TIME_BEGIN,TIME_END,SGSTIME_B(3),SGSTIME_E(3),
!     RHSTIME_B(3),RHSTIME_E(3),LHSTIME_B(3),LHSTIME_E(3),
!     POISSTIME_B(3),POISSTIME_E(3)

!===================== CHANNEL VARIABLE ================================
      INTEGER*8    :: ICH,IDPDX
      REAL*8       :: PMI,PMIAVG,PMIU,PMIL,QVOL,QVOL1,QVOL2,QFLUX,PHCAP
!     ICH,IDPDX : Channel option, dp/dx constant option
!     PMI,PMIAVG,PMIU,PMIL,QVOL,QVOL1,QVOL2,QFLUX,PHCAP : dP/dx (mean pressure gradient),
!                                                         Averaged dP/dx, dP/dx for upper wall,
!                                                         dP/dx for lower wall,
!                                                         Volume flow rate,
!                                                         volume flow rate of previous step,
!                                                         volume flow rate of present step,
!                                                         volume flux,
!                                                         correction term for volume flow rate constant simulation

!===================== AERODYNAMIC FORCE ===============================
      REAL*8       :: CDAREA
      REAL*8       :: FORCESUM(3),FORCESUMA(3)
      REAL*8       :: DUDTA,DVDTA,DWDTA
!     CDAREA                   : Coefficient to obtain drag coefficient
!     FORCESUM(3),FORCESUMA(3) : Force obtained from the momentum forcing in IBM
!     DUDTA,DVDTA,DWDTA        : Inertia contribution in IBM forcing

!======================= FIELD AVERAGE =================================
      REAL*8       :: TIMEINIT,TIMEINITZ
      INTEGER*8    :: IHISTINIT,IHISTINITZ
!     TIMEINIT,TIMEINITZ   : Start and end time for average field
!     IHISTINIT,IHISTINITZ : Start and end time step for average field

END MODULE
