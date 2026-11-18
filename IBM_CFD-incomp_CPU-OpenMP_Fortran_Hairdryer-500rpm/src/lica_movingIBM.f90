!=======================================================================
      SUBROUTINE FINDFORCING
!=======================================================================
!
!     Initialize the Immersed Boundary Method
!     (Find the immersed boundary forcing points)
!
!     Option
!       IINTP: IBM interpolation ON(1)/OFF(0)
!       ILES : LES ON(1)/OFF(0)
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
!       JTYPE(): The type of interpolation (1: 1D, 2: 2D, 3: 3D)
!       NBODY(): Number of inside points of an immersed body
!       INOUT(): 1 => Outside point, 0 => Inside point
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS

      !*  Initialize the following arrays
      INOUT     = 1      ! IN or OUT of the body (0: IN, 1: OUT)
      NFC_INTP  = 0      ! The number of interpolation points in 3 directions
      NFC_INNER = 0      ! The number of inner         points in 3 directions
      NNZ       = 0      ! The number of nut zero point
      MPI       = 0      ! The direction of from in to out. for the interpolation
      JTYPE     = 0      ! Interpolation type (1: 1D, 2: 2D, 3: 3D)
      GFI       = 0.     ! Geometric factor for interpolation

      CALL FIND_INOUT           !* Find in or out of the body
      IF (IINTP.EQ.1) THEN
         CALL FINDBODY_INTP     !* Find forcing points (Interpolation O)
         CALL INIT_SURFACE_VELOCITY
         CALL IBM_INTERPOLATION !* Calculate geometric factor for interpolation
         !CALL GET_GFIMAX
      ELSEIF (IINTP.EQ.0) THEN
         CALL FINDBODY_NOINTP   !* Find forcing points (Interpolation X)
         CALL INIT_SURFACE_VELOCITY
      ENDIF

      IF (ILES.EQ.1) THEN
         CALL UPDATE_ROTOR_CENTER_CACHE
         CALL NUTZERO_INIT      !* Find NUT zero points
         CALL WALLF_DVM         !*
      ENDIF

      RETURN
      END SUBROUTINE FINDFORCING

!=======================================================================
      SUBROUTINE FIND_INOUT
!=======================================================================
!
!     Find in or out of an immersed body
!
!     Variables
!       INOUT():
!         0 => IN,  if FUNCBODY <= 1.E-10
!         1 => OUT, if FUNCBODY >  1.E-10
!         1.E-10 instead of 0. to prevent wrong surface recognizaion
!         caused by round-off error
!       FUNCBODY(): Function or point data of body
!                   It is defined in lica_cylinder.f90 or lica_sphere.f90.
!       IUVW(),JUVW(),KUVW(): Initial I,J,K index for searching forcing points
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY : INOUT
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,L
      INTEGER*8    :: NBODY(3),IS
      REAL   *8    :: FUNCBODY

      NBODY  = 0
      IUVW(1)= 2
      IUVW(2)= 1
      IUVW(3)= 1
      JUVW(1)= 1
      JUVW(2)= 2
      JUVW(3)= 1
      KUVW(1)= 1
      KUVW(2)= 1
      KUVW(3)= 2

      IF (IPX.EQ.1) IUVW(1)=1
      IF (IPZ.EQ.1) KUVW(3)=1

      IS  = 0
!$OMP PARALLEL DO PRIVATE(I,J) REDUCTION(+:NBODY)
      DO 20 K = 0,N3
      DO 20 J = 0,N2
      DO 20 I = 0,N1
      IF (FUNCBODY(X(I),YMP(J),ZMP(K),IS) .LE. 1.E-10) THEN
         NBODY(1)=NBODY(1)+1      ! for u
         INOUT(I,J,K,1) = 0
      ENDIF
      IF (FUNCBODY(XMP(I),Y(J),ZMP(K),IS) .LE. 1.E-10) THEN
         NBODY(2)=NBODY(2)+1      ! for v
         INOUT(I,J,K,2) = 0
      ENDIF
      IF (FUNCBODY(XMP(I),YMP(J),Z(K),IS) .LE. 1.E-10) THEN
         NBODY(3)=NBODY(3)+1      ! for w
         INOUT(I,J,K,3) = 0
      ENDIF
 20   CONTINUE

!      OPEN(94,FILE='#u_in_body_2d.dat')
!      OPEN(95,FILE='#v_in_body_2d.dat')
!      OPEN(96,FILE='#w_in_body_2d.dat')
!
!      DO L=1,3
!      DO I=1,N1M
!      DO J=1,N2M
!        IF(INOUT(I,J,1,L) .EQ. 0 ) THEN
!           IF (L.EQ.1) THEN
!            WRITE(94,102) X(I),YMP(J)
!           ELSEIF (L.EQ.2) THEN
!            WRITE(95,102) XMP(I),Y(J)
!           ELSE
!            WRITE(96,102) XMP(I),YMP(J)
!           ENDIF
!        ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
! 102  FORMAT(2ES13.5)
!       CLOSE(94)
!       CLOSE(95)
!       CLOSE(96)

      WRITE(*,*) ' '
      WRITE(*,*) '******** # OF FORCING INFORMATION *******'
      WRITE(*,35) NBODY(1)
      WRITE(*,36) NBODY(2)
      WRITE(*,37) NBODY(3)
 35   FORMAT(' # of Forcing pts For U  = ',I15)
 36   FORMAT(' # of Forcing pts For V  = ',I15)
 37   FORMAT(' # of Forcing pts For W  = ',I15)

      RETURN
      END

!=======================================================================
      SUBROUTINE FINDBODY_INTP
!=======================================================================
!
!     Find inner and interpolation points
!
!     Variables
!       INOUT_SUM: the number of directions in which the INOUT value
!                  changes before and after one cell at a point.
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,L
      INTEGER*8    :: IM,IP,JM,JP,KM,KP,INOUT_SUM
      REAL   *8    :: FUNCBODY

      !* Find the interpolation points among the inner points:
      !  the point where INOUT value changes in one direction or more
      !  (INOUT_SUM> 0) while it is the point inside the body (INOUT = 0)
      DO 30 L=1,3   ! for (u,v,w)
      DO 30 K=KUVW(L),N3M
      KM=K-1
      KP=K+1
      DO 30 J=JUVW(L),N2M
      JM=JMV(J)
      JP=JPV(J)
      DO 30 I=IUVW(L),N1M
      IM=I-1
      IP=I+1

      INOUT_SUM=(ABS(INOUT(IP,J,K,L)-INOUT(IM,J,K,L)))*(1.-FIXIL(I))*(1.-FIXIU(I)) &
               +(ABS(INOUT(I,JP,K,L)-INOUT(I,JM,K,L)))*(1.-FIXJL(J))*(1.-FIXJU(J)) &
               +(ABS(INOUT(I,J,KP,L)-INOUT(I,J,KM,L)))*(1.-FIXKL(K))*(1.-FIXKU(K))

      IF (INOUT(I,J,K,L) .EQ. 0 .AND. INOUT_SUM .GT. 0) THEN
         IF(NFC_INTP(L)>=MBODY0)ERROR STOP 'IBM interpolation capacity exceeded'
         NFC_INTP(L)=NFC_INTP(L)+1
         IFC(NFC_INTP(L),L)=I
         JFC(NFC_INTP(L),L)=J
         KFC(NFC_INTP(L),L)=K
         JTYPE(NFC_INTP(L),L)=INOUT_SUM
         MPI(NFC_INTP(L),L,1)=INOUT(IP,J,K,L)+INOUT(IM,J,K,L)*(-1)
         MPI(NFC_INTP(L),L,2)=INOUT(I,JP,K,L)+INOUT(I,JM,K,L)*(-1)
         MPI(NFC_INTP(L),L,3)=INOUT(I,J,KP,L)+INOUT(I,J,KM,L)*(-1)
      ENDIF
 30   CONTINUE

      !* Find the inner points:
      !  the point where the INOUT value does not change (INOUT_SUM = 0)
      !  while it is the point inside the body (INOUT = 0)
      DO 40 L=1,3      ! for (u,v,w)
      DO 40 K=KUVW(L),N3M
         KM=K-1
         KP=K+1
      DO 40 J=JUVW(L),N2M
         JM=JMV(J)
         JP=JPV(J)
      DO 40 I=IUVW(L),N1M
         IM=I-1
         IP=I+1

      INOUT_SUM=(ABS(INOUT(IP,J,K,L)-INOUT(IM,J,K,L)))*(1-FIXIL(I))*(1-FIXIU(I)) &
               +(ABS(INOUT(I,JP,K,L)-INOUT(I,JM,K,L)))*(1-FIXJL(J))*(1-FIXJU(J)) &
               +(ABS(INOUT(I,J,KP,L)-INOUT(I,J,KM,L)))*(1-FIXKL(K))*(1-FIXKU(K))

      IF (INOUT(I,J,K,L).EQ.0 .AND. INOUT_SUM.EQ.0) THEN
         IF(NFC_INTP(L)+NFC_INNER(L)>=MBODY)ERROR STOP 'IBM inner capacity exceeded'
         NFC_INNER(L)=NFC_INNER(L)+1
         IFC(NFC_INTP(L)+NFC_INNER(L),L)=I
         JFC(NFC_INTP(L)+NFC_INNER(L),L)=J
         KFC(NFC_INTP(L)+NFC_INNER(L),L)=K
      ENDIF

 40   CONTINUE

!       OPEN(95,FILE='#u_intp_2d.dat')
!       OPEN(96,FILE='#v_intp_2d.dat')
!       OPEN(97,FILE='#w_intp_2d.dat')
!        DO L= 1,3
!           DO K=1,NFC_INTP(L)
!             IF (L.EQ.1) THEN
!               WRITE(95,103) X(IFC(K,L)),YMP(JFC(K,L))
!             ELSEIF (L.EQ.2) THEN
!               WRITE(96,103) XMP(IFC(K,L)),Y(JFC(K,L))
!             ELSEIF (L.EQ.3) THEN
!               WRITE(97,103) XMP(IFC(K,L)),YMP(JFC(K,L))
!             ENDIF
!           ENDDO
!        ENDDO
!       CLOSE(95)
!       CLOSE(96)
!       CLOSE(97)
!  103  FORMAT(2ES13.5)

      WRITE(*,*) ' '
      WRITE(*,*) '******** # OF FORCING INFORMATION *******'
      WRITE(*,51) NFC_INTP(1)
      WRITE(*,52) NFC_INTP(2)
      WRITE(*,53) NFC_INTP(3)
 51   FORMAT(' # of Intp pts For U     = ',I15)
 52   FORMAT(' # of Intp pts For V     = ',I15)
 53   FORMAT(' # of Intp pts For W     = ',I15)
      WRITE(*,54) NFC_INNER(1)
      WRITE(*,55) NFC_INNER(2)
      WRITE(*,56) NFC_INNER(3)
 54   FORMAT(' # of Inner pts For U    = ',I15)
 55   FORMAT(' # of Inner pts For V    = ',I15)
 56   FORMAT(' # of Inner pts For W    = ',I15)

      IF (MBODY0 .LT. MAX0(NFC_INTP(1),NFC_INTP(2),NFC_INTP(3))) THEN
        WRITE(*,*) 'MBODY0 .LT. NFC_INTP !!!!!'
        STOP
      ENDIF
      IF (MBODY .LT. MAX0(NFC_INTP(1)+NFC_INNER(1),NFC_INTP(2)+NFC_INNER(2),NFC_INTP(3)+NFC_INNER(3))) THEN
        WRITE(*,*) 'MBODY .LT. NBODY !!!!!'
        STOP
      ENDIF

      RETURN
      END

!=======================================================================
      SUBROUTINE FINDBODY_NOINTP
!=======================================================================
!
!     Find forcing points without interpolation
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,L
      REAL   *8    :: FUNCBODY

      !* Find the forcing points if a point is inside the body (INOUT=0)
      DO 30 L=1,3      ! for (u,v,w)
      DO 30 K=KUVW(L),N3M
      DO 30 J=JUVW(L),N2M
      DO 30 I=IUVW(L),N1M
         IF (INOUT(I,J,K,L).EQ.0) THEN
         IF(NFC_INTP(L)+NFC_INNER(L)>=MBODY)ERROR STOP 'IBM inner capacity exceeded'
         NFC_INNER(L)=NFC_INNER(L)+1
         IFC(NFC_INTP(L)+NFC_INNER(L),L)=I
         JFC(NFC_INTP(L)+NFC_INNER(L),L)=J
         KFC(NFC_INTP(L)+NFC_INNER(L),L)=K
         ENDIF
 30   CONTINUE

      WRITE(*,*) ' '
      WRITE(*,*) '******** # OF FORCING INFORMATION *******'
      WRITE(*,45) NFC_INNER(1)
      WRITE(*,46) NFC_INNER(2)
      WRITE(*,47) NFC_INNER(3)
 45   FORMAT(' # of Inner pts For U    = ',I15)
 46   FORMAT(' # of Inner pts For V    = ',I15)
 47   FORMAT(' # of Inner pts For W    = ',I15)

      RETURN
      END

!=======================================================================
      SUBROUTINE IBM_INTERPOLATION
!=======================================================================
!
!     J.I.Lee 2004.2.4; 2006.7.18
!
!     Compute geometric factor associated with interpolation
!
!     The '+' sign may be altered by '-' sign
!      GFI( , ,0,0,0) : for (I  ,J  ,K  )
!      GFI( , ,0,0,1) : for (I  ,J  ,K+1)
!      GFI( , ,0,0,2) : for (I  ,J  ,K+2)
!      GFI( , ,0,1,0) : for (I  ,J+1,K  )
!      GFI( , ,0,1,1) : for (I  ,J+1,K+1)
!      GFI( , ,0,1,2) : for (I  ,J+1,K+2)
!      GFI( , ,0,2,0) : for (I  ,J+2,K  )
!      GFI( , ,0,2,1) : for (I  ,J+2,K+1)
!      GFI( , ,0,2,2) : for (I  ,J+2,K+2)
!      GFI( , ,1,0,0) : for (I+1,J  ,K  )
!      GFI( , ,1,0,1) : for (I+1,J  ,K+1)
!      GFI( , ,1,0,2) : for (I+1,J  ,K+2)
!      GFI( , ,1,1,0) : for (I+1,J+1,K  )
!      GFI( , ,1,1,1) : for (I+1,J+1,K+1)
!      GFI( , ,1,1,2) : for (I+1,J+1,K+2)
!      GFI( , ,1,2,0) : for (I+1,J+2,K  )
!      GFI( , ,1,2,1) : for (I+1,J+2,K+1)
!      GFI( , ,1,2,2) : for (I+1,J+2,K+2)
!      GFI( , ,2,0,0) : for (I+2,J  ,K  )
!      GFI( , ,2,0,1) : for (I+2,J  ,K+1)
!      GFI( , ,2,0,2) : for (I+2,J  ,K+2)
!      GFI( , ,2,1,0) : for (I+2,J+1,K  )
!      GFI( , ,2,1,1) : for (I+2,J+1,K+1)
!      GFI( , ,2,1,2) : for (I+2,J+1,K+2)
!      GFI( , ,2,2,0) : for (I+2,J+2,K  )
!      GFI( , ,2,2,1) : for (I+2,J+2,K+1)
!      GFI( , ,2,2,2) : for (I+2,J+2,K+2)
!
!      The paper 'An Immersed boundary finite...' written by Kim et al.
!      selected the no slip condition point for 2D & 3D interpolation
!      so that the line connecting the no slip condition point and
!      the forcing point to be perpendicular to the body surface.
!      But this code find the no slip condition point on the line
!      starting from the forcing point and whose direction vector is
!      (MPI( , ,3)MPI( , ,3)MPI( , ,3))
!
!      This is because of the following reason
!      1. In the two schemes of determing the no slip condition point
!         one is not better than the other in view of numerical
!         accuracy or stablilty
!      2. The MPI direction method is easier to embody in code the
!         the perpendiculr method
!      3. The perpendicular point may not exist in peculiar case
!         but the MPI direction point always exists
!
!      Note that
!      1. FUNCBODY function requrements
!      2. Secant method iteration constants by trial & error
!      (See below..)
!
!      2006.7.18
!      Generalized IBM Interpolation Scheme
!      which Guarantees GFI<=1. (in uniform grid spacing)
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,L,N,II,JJ,KK
      INTEGER*8    :: IPAIR,JPAIR,KPAIR
      REAL*8       :: XX(-1:M1+1,3),YY(0:M2,3),ZZ(-1:M3+1,3)
      REAL*8       :: P1,P2,P3
      REAL*8       :: X0,Y0,Z0,X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3
      INTEGER*8    :: MMAX,MMM,MM,MAB!,ITR
      REAL*8       :: XTEMP1,YTEMP1,ZTEMP1,XTEMP2,YTEMP2,ZTEMP2,FFS,FFE,FF1,FF2
      REAL*8       :: XX1,YY1,ZZ1,XX2,YY2,ZZ2,DDX,DDY,DDZ
      REAL*8       :: DX1,DX2,DX3,DY1,DY2,DY3,DZ1,DZ2,DZ3
      REAL*8       :: A0,B0,C0,A1,B1,C1
      REAL*8       :: FUNCBODY

!==== Define new coordinate array for shortening coding
!     (to prevent memory overhead)
      MAB = 0
      !OPEN(55,FILE='u_surf_2D.dat')
      !OPEN(66,FILE='v_surf_2D.dat')
      !OPEN(77,FILE='w_surf_2D.dat')

      DO I=1,N1
      XX(I,1)=X(I)     ! for u
      XX(I,2)=XMP(I)   ! for v
      XX(I,3)=XMP(I)   ! for w
      ENDDO
      XX(0,2)=XMP(0)
      XX(0,3)=XMP(0)

      IF(IPX .EQ. 1) THEN
        XX(0,1)=X(1)-(X(N1)-X(N1M))
        XX(-1,1)=X(1)-(X(N1)-X(N1-2))
        XX(N1+1,1)=X(N1)+(X(2)-X(1))

        XX(0,2)=X(1)-0.5*(X(N1)-X(N1M))
        XX(-1,2)=X(1)-(X(N1)-X(N1M))-0.5*(X(N1M)-X(N1M-1))
        XX(N1,2)=X(N1)+0.5*(X(2)-X(1))
        XX(N1+1,2)=X(N1)+X(2)-X(1)+0.5*(X(3)-X(2))

        XX(0,3)=X(1)-0.5*(X(N1)-X(N1M))
        XX(-1,3)=X(1)-(X(N1)-X(N1M))-0.5*(X(N1M)-X(N1M-1))
        XX(N1,3)=X(N1)+0.5*(X(2)-X(1))
        XX(N1+1,3)=X(N1)+(X(2)-X(1))+0.5*(X(3)-X(2))
      ENDIF

      DO J=1,N2
      YY(J,1)=YMP(J)
      YY(J,2)=Y(J)
      YY(J,3)=YMP(J)
      ENDDO
      YY(0,1)=YMP(0)
      YY(0,3)=YMP(0)

      DO K=1,N3
      ZZ(K,1)=ZMP(K)
      ZZ(K,2)=ZMP(K)
      ZZ(K,3)=Z(K)
      ENDDO
      ZZ(0,1)=ZMP(0)
      ZZ(0,2)=ZMP(0)

      IF(IPZ .EQ. 1) THEN
        ZZ(0,3)=Z(1)-(Z(N3)-Z(N3M))
        ZZ(-1,3)=Z(1)-(Z(N3)-Z(N3-2))
        ZZ(N3+1,3)=Z(N3)+(Z(2)-Z(1))

        ZZ(0,1)=Z(1)-0.5*(Z(N3)-Z(N3M))
        ZZ(-1,1)=Z(1)-(Z(N3)-Z(N3M))-0.5*(Z(N3M)-Z(N3M-1))
        ZZ(N3,1)=Z(N3)+0.5*(Z(2)-Z(1))
        ZZ(N3+1,1)=Z(N3)+Z(2)-Z(1)+0.5*(Z(3)-Z(2))


        ZZ(0,2)=Z(1)-0.5*(Z(N3)-Z(N3M))
        ZZ(-1,2)=Z(1)-(Z(N3)-Z(N3M))-0.5*(Z(N3M)-Z(N3M-1))
        ZZ(N3,2)=Z(N3)+0.5*(Z(2)-Z(1))
        ZZ(N3+1,2)=Z(N3)+(Z(2)-Z(1))+0.5*(Z(3)-Z(2))

      ENDIF

      DO 500 L=1,3               ! L : component (u,v,w)
      DO 500 N=1,NFC_INTP(L)

         II=IFC(N,L)
         JJ=JFC(N,L)
         KK=KFC(N,L)
         X1=XX(II,L)
         Y1=YY(JJ,L)
         Z1=ZZ(KK,L)
         P1=FLOAT(MPI(N,L,1))    ! To simplify code statements
         P2=FLOAT(MPI(N,L,2))
         P3=FLOAT(MPI(N,L,3))

         IPAIR=II+MPI(N,L,1)
         JPAIR=JJ+MPI(N,L,2)
         KPAIR=KK+MPI(N,L,3)
         X2=XX(IPAIR,L)
         Y2=YY(JPAIR,L)
         Z2=ZZ(KPAIR,L)

      MMAX = 3

      DO MMM=1,MMAX

        IF(MMM .EQ. 1) THEN
          XTEMP1=X1
          YTEMP1=Y1
          ZTEMP1=Z1
          XTEMP2=X2
          YTEMP2=Y2
          ZTEMP2=Z2
          FFS=FUNCBODY(X1,Y1,Z1,MAB)
          FFE=FUNCBODY(X2,Y2,Z2,MAB)
          IF(FFS*FFE .GT. 0.) THEN
            DO I=0,2
              DO J=0,2
                DO K=0,2
                  GFI(N,L,I,J,K) = 0.
                ENDDO
              ENDDO
            ENDDO
            GFI(N,L,0,0,0)= 1.
            GOTO 45
          ENDIF
        ELSE
          XTEMP1=XX1
          YTEMP1=YY1
          ZTEMP1=ZZ1
          XTEMP2=XX2
          YTEMP2=YY2
          ZTEMP2=ZZ2
        ENDIF

        DO MM= 0,19
          DDX=(XTEMP2-XTEMP1)
          DDY=(YTEMP2-YTEMP1)
          DDZ=(ZTEMP2-ZTEMP1)

          XX1=XTEMP1+DDX*MM/20.
          XX2=XTEMP1+DDX*(MM+1)/20.
          YY1=YTEMP1+DDY*MM/20.
          YY2=YTEMP1+DDY*(MM+1)/20.
          ZZ1=ZTEMP1+DDZ*MM/20.
          ZZ2=ZTEMP1+DDZ*(MM+1)/20.

          FF1=FUNCBODY(XX1,YY1,ZZ1,MAB)
          FF2=FUNCBODY(XX2,YY2,ZZ2,MAB)

          IF(FF1 .EQ. 0.) THEN
            X0=XX1
            Y0=YY1
            Z0=ZZ1
            GOTO 33
          ELSEIF(FF2 .EQ. 0.) THEN
            X0=XX2
            Y0=YY2
            Z0=ZZ2
            GOTO 33
          ELSEIF(FF1*FF2 .LT. 0.) THEN
            X0= 0.5*(XX1+XX2)
            Y0= 0.5*(YY1+YY2)
            Z0= 0.5*(ZZ1+ZZ2)
            IF(MMM .EQ. MMAX) GOTO 33
            GOTO 22
          ENDIF
        ENDDO

 22     CONTINUE
      ENDDO
 33   CONTINUE
      CALL SET_SURFACE_VELOCITY(N,L,X1,Y1,Z1,X0,Y0,Z0,X2,Y2,Z2)

!     Generalized Interpolation Scheme
!     This scheme guarantees thant GFI<=1.
!     2006.7.15

      X3=XX(II+MPI(N,L,1)*2,L)
      Y3=YY(JJ+MPI(N,L,2)*2,L)
      Z3=ZZ(KK+MPI(N,L,3)*2,L)

      DX1= ABS(X1-X0)
      DX2= ABS(X2-X0)
      DX3= ABS(X3-X0)
      DY1= ABS(Y1-Y0)
      DY2= ABS(Y2-Y0)
      DY3= ABS(Y3-Y0)
      DZ1= ABS(Z1-Z0)
      DZ2= ABS(Z2-Z0)
      DZ3= ABS(Z3-Z0)

      IF (MPI(N,L,1) == 0) THEN
         A0=1.
         A1=1.
      ELSE IF (DX2 >= DX1) THEN
         A0=DX2/(DX1+DX2)
         A1=1.
      ELSE IF (DX2 <  DX1) THEN
         A0=0.5
         A1=(DX3-DX1)/(DX3-DX2)
      END IF

      IF (MPI(N,L,2) == 0) THEN
         B0=1.
         B1=1.
      ELSE IF (DY2 >= DY1) THEN
         B0=DY2/(DY1+DY2)
         B1= 1.
      ELSE IF (DY2 <  DY1) THEN
         B0= 0.5
         B1=(DY3-DY1)/(DY3-DY2)
      END IF

      IF (MPI(N,L,3) == 0) THEN
         C0=1.
         C1=1.
      ELSE IF (DZ2 >= DZ1) THEN
         C0=DZ2/(DZ1+DZ2)
         C1=1.
      ELSE IF (DZ2 <  DZ1) THEN
         C0=0.5
         C1=(DZ3-DZ1)/(DZ3-DZ2)
      END IF

      GFI(N,L,0,0,0)= 1./A0/B0/C0
      GFI(N,L,0,0,1)=-1./A0/B0/C0*(   A0)*(   B0)*(1.-C0)    &
                                 *                (   C1)
      GFI(N,L,0,0,2)=-1./A0/B0/C0*(   A0)*(   B0)*(1.-C0)    &
                                 *                (1.-C1)
      GFI(N,L,0,1,0)=-1./A0/B0/C0*(   A0)*(1.-B0)*(   C0)    &
                                 *        (   B1)
      GFI(N,L,0,1,1)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                 *        (   B1)*(   C1)
      GFI(N,L,0,1,2)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                 *        (   B1)*(1.-C1)
      GFI(N,L,0,2,0)=-1./A0/B0/C0*(   A0)*(1.-B0)*(   C0)    &
                                 *        (1.-B1)
      GFI(N,L,0,2,1)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                 *        (1.-B1)*(   C1)
      GFI(N,L,0,2,2)=-1./A0/B0/C0*(   A0)*(1.-B0)*(1.-C0)    &
                                 *        (1.-B1)*(1.-C1)
      GFI(N,L,1,0,0)=-1./A0/B0/C0*(1.-A0)*(   B0)*(   C0)    &
                                 *(   A1)
      GFI(N,L,1,0,1)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                 *(   A1)        *(   C1)
      GFI(N,L,1,0,2)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                 *(   A1)        *(1.-C1)
      GFI(N,L,1,1,0)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                 *(   A1)*(   B1)
      GFI(N,L,1,1,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(   A1)*(   B1)*(   C1)
      GFI(N,L,1,1,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(   A1)*(   B1)*(1.-C1)
      GFI(N,L,1,2,0)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                 *(   A1)*(1.-B1)
      GFI(N,L,1,2,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(   A1)*(1.-B1)*(   C1)
      GFI(N,L,1,2,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(   A1)*(1.-B1)*(1.-C1)
      GFI(N,L,2,0,0)=-1./A0/B0/C0*(1.-A0)*(   B0)*(   C0)    &
                                 *(1.-A1)
      GFI(N,L,2,0,1)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                 *(1.-A1)        *(   C1)
      GFI(N,L,2,0,2)=-1./A0/B0/C0*(1.-A0)*(   B0)*(1.-C0)    &
                                 *(1.-A1)        *(1.-C1)
      GFI(N,L,2,1,0)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                 *(1.-A1)*(   B1)
      GFI(N,L,2,1,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(1.-A1)*(   B1)*(   C1)
      GFI(N,L,2,1,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(1.-A1)*(   B1)*(1.-C1)
      GFI(N,L,2,2,0)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(   C0)    &
                                 *(1.-A1)*(1.-B1)
      GFI(N,L,2,2,1)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(1.-A1)*(1.-B1)*(   C1)
      GFI(N,L,2,2,2)=-1./A0/B0/C0*(1.-A0)*(1.-B0)*(1.-C0)    &
                                 *(1.-A1)*(1.-B1)*(1.-C1)

!         IF (ABS(GFI(N,L,0,2,0)) .GT. 1.) THEN
!           WRITE(653,'(3ES13.5)') DY1,DY2,DY3
!         ENDIF

 45   CONTINUE

!      IF (KK.EQ.1) THEN
!
!        IF (L.EQ.1) THEN
!         WRITE(55,555) X0,Y0
!        ENDIF
!        IF (L.EQ.2) THEN
!         WRITE(66,555) X0,Y0
!        ENDIF
!        IF (L.EQ.3) THEN
!         WRITE(77,555) X0,Y0
!        ENDIF
!      ENDIF
! 555  FORMAT(2ES13.5)


 500  CONTINUE
      CALL REBUILD_FLUID_WALL_STENCILS


      RETURN
      END

!!=======================================================================
!        SUBROUTINE GET_GFIMAX
!!=======================================================================
!       USE MOD_IBMPRE
!       IMPLICIT NONE
!       INTEGER*8    :: I,J,K,L,N,M,LL,MM
!       REAL*8       :: TMPGFI(MBODY0,3)
!       REAL*8       :: GFIMAX,GFIMAX_JTYPE1,GFIMAX_JTYPE2,GFIMAX_JTYPE3
!       INTEGER*8    :: NGFIGT1,NJTYPE1,NJTYPE2,NJTYPE3,NFC_INTP_SUM
!
!
!        DO L=1,3
!        DO N=1,NFC_INTP(L)
!           TMPGFI(N,L)    = GFI(N,L,0,0,0)
!           GFI(N,L,0,0,0) = 0.
!        ENDDO
!        ENDDO
!
!        GFIMAX = 0.
!        GFIMAX_JTYPE1 = 0.
!        GFIMAX_JTYPE2 = 0.
!        GFIMAX_JTYPE3 = 0.
!        NGFIGT1 = 0
!        NJTYPE1 = 0
!        NJTYPE2 = 0
!        NJTYPE3 = 0
!        NFC_INTP_SUM=NFC_INTP(1)+NFC_INTP(2)+NFC_INTP(3)
!
!        DO L=1,3
!        DO M=1,NFC_INTP(L)
!        DO K=0,2
!        DO J=0,2
!        DO I=0,2
!
!           IF ((ABS(GFI(M,L,I,J,K)) .GT. GFIMAX_JTYPE1).AND. (JTYPE(M,L) .EQ. 1)) THEN
!              GFIMAX_JTYPE1= ABS(GFI(M,L,I,J,K))
!           ENDIF
!           IF ((ABS(GFI(M,L,I,J,K)) .GT. GFIMAX_JTYPE2).AND. (JTYPE(M,L) .EQ. 2)) THEN
!              GFIMAX_JTYPE2= ABS(GFI(M,L,I,J,K))
!           ENDIF
!           IF ((ABS(GFI(M,L,I,J,K)) .GT. GFIMAX_JTYPE3).AND. (JTYPE(M,L) .EQ. 3)) THEN
!              GFIMAX_JTYPE3= ABS(GFI(M,L,I,J,K))
!           ENDIF
!
!           IF (ABS(GFI(M,L,I,J,K)) .GT. GFIMAX) THEN
!              GFIMAX=ABS(GFI(M,L,I,J,K))
!              LL=L
!              MM=M
!           ENDIF
!
!           IF (ABS(GFI(M,L,I,J,K)) .GT. 1.) NGFIGT1=NGFIGT1+1
!           IF (JTYPE(M,L) .EQ. 1) NJTYPE1=NJTYPE1+1
!           IF (JTYPE(M,L) .EQ. 2) NJTYPE2=NJTYPE2+1
!           IF (JTYPE(M,L) .EQ. 3) NJTYPE3=NJTYPE3+1
!
!        ENDDO
!        ENDDO
!        ENDDO
!        ENDDO
!        ENDDO
!
!        WRITE(*,*) ' '
!        WRITE(*,*) '***** IBM INTERPOLATION INFORMATION *****'
!        WRITE(*,101) GFIMAX
!        WRITE(*,102) JTYPE(MM,LL)
!        WRITE(*,103) X(IFC(MM,LL)),IFC(MM,LL)
!        WRITE(*,104) Y(JFC(MM,LL)),JFC(MM,LL)
!        WRITE(*,105) Z(KFC(MM,LL)),KFC(MM,LL)
!        WRITE(*,106) NGFIGT1
!        WRITE(*,107) FLOAT(NGFIGT1)/FLOAT(26*NFC_INTP_SUM)*100.
!        WRITE(*,108) GFIMAX_JTYPE1
!        WRITE(*,109) GFIMAX_JTYPE2
!        WRITE(*,110) GFIMAX_JTYPE3
!        WRITE(*,111) NJTYPE1
!        WRITE(*,112) NJTYPE2
!        WRITE(*,113) NJTYPE3
!
!
!
!        DO L=1,3
!        DO N=1,NFC_INTP(L)
!           GFI(N,L,0,0,0) = TMPGFI(N,L)
!        ENDDO
!        ENDDO
!
!
!  101   FORMAT(' GFIMAX                  = ',F15.8)
!  102   FORMAT(' JTYPE OF GFIMAX         = ',I15)
!  103   FORMAT(' X(IFC),IFC              = ',F15.8,I6)
!  104   FORMAT(' Y(JFC),JFC              = ',F15.8,I6)
!  105   FORMAT(' Z(KFC),KFC              = ',F15.8,I6)
!  106   FORMAT(' # OF GFI>1              = ',I15)
!  107   FORMAT(' % OF GFI>1              = ',F15.8)
!  108   FORMAT(' GFIMAX(JTYPE=1)         = ',F15.8)
!  109   FORMAT(' GFIMAX(JTYPE=2)         = ',F15.8)
!  110   FORMAT(' GFIMAX(JTYPE=3)         = ',F15.8)
!  111   FORMAT(' # OF JTYPE=1            = ',I15)
!  112   FORMAT(' # OF JTYPE=2            = ',I15)
!  113   FORMAT(' # OF JTYPE=3            = ',I15)
!        WRITE(*,*) ' '
!
!        RETURN
!        END

!=======================================================================
      SUBROUTINE UPDATE_ROTOR_CENTER_CACHE
      USE MOD_COMMON
      USE MOD_FLOWARRAY, ONLY: ROTOR_PHI_CENTER
      USE MOD_ROTOR_MOTION, ONLY: ROTOR_SDF
      IMPLICIT NONE
      INTEGER*8 I,J,K
!$OMP PARALLEL DO PRIVATE(I,J)
      DO K=0,N3
      DO J=0,N2
      DO I=0,N1
        ROTOR_PHI_CENTER(I,J,K)=ROTOR_SDF(XMP(I),YMP(J),ZMP(K))
      ENDDO
      ENDDO
      ENDDO
      END SUBROUTINE UPDATE_ROTOR_CENTER_CACHE

      SUBROUTINE NUTZERO_INIT
!=======================================================================
!
!     Find points to make the sub-grid scale eddy viscosity zero
!     inside and around the IB body.
!
!     NNZ: Number of TNU zero point
!     INZ, JNZ, KNZ: Index of TNU zero point
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,N,L
      !INTEGER*8    :: INZ(MBODY),JNZ(MBODY),KNZ(MBODY)
      REAL*8       :: FUNCBODY

      DO 10 K=1,N3M
      DO 10 J=1,N2M
      DO 10 I=1,N1M
      IF ((ROTOR_PHI_CENTER(I  ,J  ,K  ) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I-1,J  ,K  ) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I+1,J  ,K  ) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I  ,J-1,K  ) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I  ,J+1,K  ) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I  ,J  ,K-1) .LE. 1.E-10) .OR. &
          (ROTOR_PHI_CENTER(I  ,J  ,K+1) .LE. 1.E-10)) THEN

         NNZ=NNZ+1
         INZ(NNZ)=I
         JNZ(NNZ)=J
         KNZ(NNZ)=K
      ENDIF
 10   CONTINUE

      WRITE(*,30) NNZ
 30   FORMAT(' # of Nut  =0 Pts        = ',I15)

!      OPEN(13,FILE='ibmpre_nutzero.bin',FORM='UNFORMATTED')
!       WRITE(13) NNZ
!       WRITE(13) (INZ  (N),N=1,NNZ)
!       WRITE(13) (JNZ  (N),N=1,NNZ)
!       WRITE(13) (KNZ  (N),N=1,NNZ)
!      CLOSE(13)

      RETURN
      END

!=======================================================================
      SUBROUTINE WALLF_DVM
!=======================================================================
!
!     Find points to make the sub-grid scale eddy viscosity zero
!     inside and around the IB body.
!
!     NUTZERO_INIT and WALLF_DVM subroutines find same points, but
!     they use different types of array.
!       NUTZERO_INIT => 1D array
!       WALLF_DVM    => 3D array
!
!-----------------------------------------------------------------------
      USE MOD_COMMON
      USE MOD_FLOWARRAY
      USE MOD_IBM_SURFACE_VELOCITY, ONLY: INIT_SURFACE_VELOCITY,SET_SURFACE_VELOCITY,REBUILD_FLUID_WALL_STENCILS
      IMPLICIT NONE
      INTEGER*8    :: I,J,K,II,JJ,KK
      INTEGER*8    :: NNW,NWON
      !INTEGER*8    :: NWALL_DVM(M1M,M2M,M3M)
      REAL*8       :: FUNCBODY

      NWALL_DVM = 1

      NNW = 0
      DO K=1,N3M
      DO J=2,N2M-1
      DO I=2,N1M-1

       NWON = 0
       DO II=-1,1
          IF (ROTOR_PHI_CENTER(I+II,J,K) .LE. 1.E-10) NWON = 1
       ENDDO
       DO JJ=-1,1
          IF (ROTOR_PHI_CENTER(I,J+JJ,K) .LE. 1.E-10) NWON = 1
       ENDDO

       IF (NWON .EQ. 1) THEN
          NNW=NNW+1
          NWALL_DVM(I,J,K) = 0
       ENDIF

      ENDDO
      ENDDO
      ENDDO

      WRITE(*,10) NNW
 10   FORMAT(' # of Zero-Filter Pts    = ',I15)

!      OPEN(14,FILE='ibmpre_wallfdvm.bin',FORM='UNFORMATTED')
!      WRITE(14) ((( NWALL_DVM(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)
!      CLOSE(14)

      RETURN
      END

!=======================================================================
      SUBROUTINE READCURVE
!=======================================================================
!
!     Read Lagrangian points of geometry/body
!     : Input lagrangian points should be distributed in the
!       clockwise direction.
!
!     Variables
!       ND1,ND2    : # of Lagrangian points of upper/lower boundary
!       XFU,YFU,ZFU: Upper boundary positions of the body
!       XFL,YFL,ZFL: Lower
!       WTH        : Half of body thickness (IPZ=0)
!                    2*ZL                   (IPZ=1)
!
!-----------------------------------------------------------------------
      USE MOD_MOVINGIBM

      IMPLICIT NONE
      REAL*8    PI!,TP
      INTEGER*8 I

      OPEN(891,FILE='data_u.dat')
      OPEN(892,FILE='data_l.dat')

      TP=2*PI

      WTH=1D0

      READ(891,*) ND1
      DO I=1,ND1
      READ(891,*) XFU(I),YFU(I)
      ENDDO

      READ(892,*) ND2
      DO I=1,ND2
      READ(892,*) XFL(I),YFL(I)
      ENDDO

      CLOSE(891)
      CLOSE(892)

      RETURN
      END

!=======================================================================
      FUNCTION XC(T)
!=======================================================================
!
!     Translational distance in x-direction
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,XC

      XC=0D0

      RETURN
      END

!=======================================================================
      FUNCTION YC(T)
!=======================================================================
!
!     Translational distance in y-direction
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,YC
        REAL*8    PI

        PI = ACOS(-1.)

        YC=0.2D0*DCOS(2*PI*0.193*T)

      RETURN
      END

!=======================================================================
      FUNCTION ZC(T)
!=======================================================================
!
!     Translational distance in z-direction
!
!-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 T,ZC

      ZC=0D0

      RETURN
      END

!=======================================================================

      FUNCTION UBD(X,Y,Z)
      IMPLICIT NONE
      REAL*8 X,Y,Z,UBD
      UBD=0D0
      END FUNCTION UBD

      FUNCTION VBD(X,Y,Z)
      USE MOD_ROTOR_MOTION, ONLY: OMEGA_NOW,ROTATING_OWNER
      IMPLICIT NONE
      REAL*8 X,Y,Z,VBD
      VBD=0D0
      IF(ROTATING_OWNER(X,Y,Z)) VBD=-OMEGA_NOW*Z
      END FUNCTION VBD

      FUNCTION WBD(X,Y,Z)
      USE MOD_ROTOR_MOTION, ONLY: OMEGA_NOW,ROTATING_OWNER
      IMPLICIT NONE
      REAL*8 X,Y,Z,WBD
      WBD=0D0
      IF(ROTATING_OWNER(X,Y,Z)) WBD=OMEGA_NOW*Y
      END FUNCTION WBD
