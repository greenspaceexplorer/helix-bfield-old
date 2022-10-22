PROGRAM HELIXTEST
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER ::     &
      VOUT = "(5X,A14,3ES12.4E2,A1)",&
      BOUT1 = "(BX BY BZ) = (",     &
      XOUT1 = "( X  Y  Z) = (",     &
      VOUT2 = ")"
  REAL*8,PARAMETER :: TOT = 1.D-4,TOM = 1.D-2
  REAL*8 :: I,A
  REAL*8 :: NTURNS
  REAL*8 :: X,Y,Z,RHO
  REAL*8 :: BRHO,BZ,MODB
  REAL*8 :: BX,BY
  

! Test individual coil subroutine
  WRITE(*,*) "Single coil test:"
  I = 90.09D0
  A = 25.5524D0
  NTURNS = 1.D0
  
  ! First test with on-axis coordinates

  X = 0.D0
  Y = 0.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)

  CALL COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X*TOM,Y*TOM,Z*TOM,VOUT2
  WRITE(*,VOUT) BOUT1,BX*TOT,BY*TOT,BZ*TOT,VOUT2
  ! Next test with off axis coordinates

  X = 5.D0
  Y = 7.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)

  CALL COILB(I,A,NTURNS,X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X*TOM,Y*TOM,Z*TOM,VOUT2
  WRITE(*,VOUT) BOUT1,BX*TOT,BY*TOT,BZ*TOT,VOUT2

! Test field for magnet with no shifts
  WRITE(*,*)
  WRITE(*,*) "Full magnet test:"

  X = 0.D0
  Y = 0.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)
  CALL HEAT_MAG(X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X*TOM,Y*TOM,Z*TOM,VOUT2
  WRITE(*,VOUT) BOUT1,BX*TOT,BY*TOT,BZ*TOT,VOUT2
  ! Next test with off axis coordinates

  X = 5.D0
  Y = 7.D0
  Z = 10.D0

  RHO = SQRT(X*X+Y*Y)

  CALL HEAT_MAG(X,Y,Z,BRHO,BZ,MODB)

  IF(RHO.EQ.0.D0)THEN
      BX = BRHO
      BY = BRHO
  ELSE
      BX = BRHO*X/RHO
      BY = BRHO*Y/RHO
  ENDIF
  WRITE(*,*)
  WRITE(*,VOUT) XOUT1,X*TOM,Y*TOM,Z*TOM,VOUT2
  WRITE(*,VOUT) BOUT1,BX*TOT,BY*TOT,BZ*TOT,VOUT2


  STOP
END PROGRAM HELIXTEST
