MODULE MOD_MOVINGIBM

!==================== BODY POINTS ======================================
      INTEGER*8, PARAMETER :: NCP= 1000
      INTEGER*8 :: ND1,ND2
      REAL*8    :: XFU(NCP),XFL(NCP),YFU(NCP),YFL(NCP),ZFU(NCP),ZFL(NCP)
      REAL*8    :: FD,WTH,TP
!     NCP        : Maximum number of Lagrangian points of geometry/body
!     ND1,ND2    : # of Lagrangian points of upper/lower boundary
!     XFU,YFU,ZFU: Upper boundary positions of the body
!     XFL,YFL,ZFL: Lower
!     WTH        : Half of body thickness (IPZ=0)
!                  2*ZL                   (IPZ=1)

END MODULE