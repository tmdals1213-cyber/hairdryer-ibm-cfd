!======================================================================
      SUBROUTINE SIDEWALLBC
	  !! Boundary conditions at the "side walls"
		! JxT : upper wall (y-direction)
		! JxB : lower wall
		! KxT : anterior wall (z-direction)
		! KxB : posterior wall
      !! Options included
	    ! y-direction : Dirichlet / Neumann 
		! z-direction : Dirichlet / Neumann (3D)
        !               Periodic (2D,3D)	    
		
				
		! IPX : periodic in x-direction
		! ICH : For channel flow
		! IDPDX : RE_TAU Constant	  
!======================================================================
       USE MOD_COMMON
       IMPLICIT NONE

        JUT = 1         ! 0: DIRICHLET
        JUB = 1         ! 1: NEUMANN
        JWT = 1
        JWB = 1
        JVT = 0
        !IF IPZ=1, THEN KUT=KUB=KVT=KVB=0
        KUT = 1
        KUB = 1
        KVT = 1
        KVB = 1

        IPX = 0        ! 0 : NOT PREODIC
        IPZ = 0        ! 1 : PREODIC

        ICH = 0        ! 0 : CHANNEL X
                       ! 1 : CHANNEL O (FOR THE CONST. MASS FLOW RATE)

        IDPDX = 0      ! 1 : RETAU CONSTANT 
          
        IF (IPZ.EQ.1) THEN
        IF (KUT.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KUT'
            STOP
        ENDIF
        IF (KUB.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KUB'
            STOP
        ENDIF
        IF (KVT.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KVT'
            STOP
        ENDIF
        IF (KVB.NE.0) THEN
            PRINT*,'INCONSISTENCY BETWEEN IPZ & KVB'
            STOP
        ENDIF
        ENDIF


      WRITE(*,101) JUT,JUB,JWT,JWB
      WRITE(*,102) KUT,KUB,KVT,KVB
      WRITE(*,103) JVT
      WRITE(*,104) IPX,IPZ
      WRITE(*,105) ICH

 101  FORMAT('JUT= ',I4,'  JUB= ',I4,'  JWT= ',I4,'  JWB= ',I4)
 102  FORMAT('KUT= ',I4,'  KUB= ',I4,'  KVT= ',I4,'  KVB= ',I4)
 103  FORMAT('JVT= ',I4)
 104  FORMAT('IPX= ',I4,'  IPZ= ',I4)
 105  FORMAT('ICH= ',I4)

      RETURN
      END

!=======================================================================
        SUBROUTINE GETCDAREA
		! CDAREA = non-dimensionalizd frontal area of a sphere
		! FRONTAL AREA of a sphere : PI*D**2/4
!=======================================================================
       USE MOD_COMMON
       IMPLICIT NONE

       CDAREA= ACOS(-1D0)/4D0 ! Reference projected rotor disc area

      RETURN
      END

!=======================================================================
        SUBROUTINE CONVBC
        USE MOD_COMMON
        IMPLICIT NONE
        DUOUT=0D0;DVOUT=0D0;DWOUT=0D0;DVTOP=0D0
        RETURN
        END
!=======================================================================
        SUBROUTINE MAKEFLD
	  !! Make an initial field
	  !  If you select Dirichlet boundary conditions at the side walls,
	  !  (See 'SUBROUTINE SIDEWALLBC')
	  !  values applied at this subroutine are imposed
	  !  at the side walls until the end of code		
!=======================================================================
       USE MOD_ROTOR_MOTION, ONLY: INFLOW_RATIO
       USE MOD_COMMON
       USE MOD_FLOWARRAY, ONLY : U,V,W,P 
       IMPLICIT NONE
       INTEGER*8     :: I,J,K
       REAL*8 :: DEMO_PSI

       IHIST= 0
       TIME = 0.
       U=0D0;V=0D0;W=0D0;P=0D0 ! quiescent reservoir including ghosts

!$OMP PARALLEL DO
        DO 10 K=0,N3
        DO 10 J=0,N2
        DO 10 I=1,N1
 10     U(I,J,K) = INFLOW_RATIO

!$OMP PARALLEL DO
        DO 20 K=0,N3
        DO 20 J=1,N2
        DO 20 I=0,N1
 20     V(I,J,K) = 0.

!$OMP PARALLEL DO
        DO 30 K=1,N3
        DO 30 J=0,N2
        DO 30 I=0,N1
 30     W(I,J,K) = 0.

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

        WRITE(*,*)' '
        WRITE(*,100)
  100   FORMAT('INITIAL FIELD:      MAKING DONE')

      RETURN
      END

!=======================================================================
      SUBROUTINE REAL_TIME(SECONDS)
!=======================================================================
      IMPLICIT NONE
      INTEGER CLOCK_COUNT,CLOCK_MAX,CLOCK_RATE
      REAL SECONDS

      CALL SYSTEM_CLOCK ( CLOCK_COUNT, CLOCK_RATE, CLOCK_MAX )

      SECONDS=REAL(CLOCK_COUNT)/REAL(CLOCK_RATE)

      RETURN
      END

!=======================================================================
      subroutine timestamp ( )
!=======================================================================
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!  Example:
!    31 May 2001   9:45:54.872 AM
!  Modified:
!    06 August 2005
!  Author: John Burkardt
!  Parameters: None
        implicit none
        character ( len = 8 ) ampm
        integer d
        integer h
        integer m
        integer mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/  &
          'January  ', 'February ', 'March    ', 'April    ',          &
          'May      ', 'June     ', 'July     ', 'August   ',          &
          'September', 'October  ', 'November ', 'December ' /)
        integer n
        integer s
        integer values(8)
        integer y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )  &
      d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      write(*,*)' '

        return
      end

!=======================================================================
! Replace the existing FUNCBODY definition, do not link both definitions.
function FUNCBODY(x,y,z,is) result(f)
  use iso_fortran_env, only: real64, int64
  use mod_rotor_motion, only: rotor_sdf
  implicit none
  real(real64), intent(in) :: x,y,z
  integer(int64), intent(in) :: is
  real(real64) :: f
  f=rotor_sdf(x,y,z)
end function

      REAL*8 FUNCTION DEMO_PSI(XX,YY,ZZ)
      USE MOD_COMMON, ONLY: EPS_PTR
      IMPLICIT NONE
      REAL*8 :: XX,YY,ZZ
      DEMO_PSI=EPS_PTR*EXP(-((XX-0.85D0)/0.4D0)**2 &
               -((YY-0.17D0)/0.25D0)**2-((ZZ-0.03D0)/0.25D0)**2)
      END FUNCTION DEMO_PSI
