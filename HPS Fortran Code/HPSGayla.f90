! ****************************************************************************
!        	HORIZONTAL PLANE SOURCE (HPS) MODEL
!	     ANALYTICAL THREE DIMENSINAL SOLUTION FOR
!	         GROUND-WATER CONTAMINANT TRANSPORT
!
!      INPUT PARAMETERS:
!      NUMPT = NUMBER OF CALCULATION POINTS
!      NINT = NUMBER OF TIMESTEPS FOR INTEGRATION
!      NUMST = NUMBER OF CHANGES IN SOURCE STRENGTH
!                   (1 FOR CONTINUOUS SOURCE)
!     NUMSRS = NUMBER OF SOURCE AREAS
!     M (NUMSRS, NUMST) = SOURCE STRENGTH (GRAM/YEAR) -
!	                              (NUMSRS X NUMST) VALUES
!      TS (NUMST) = TIMES (YEARS) FOR CHANGES IN SOURCE STRENGTH -
!	                                                 NUMST VALUES
!      XS (NUMSRS), YS( NUMSRS), ZS ( NUMSRS ) = X, Y, Z  LOCATIONS (METERS)
!                                             OF SOURCE AREA CENTERS - NUMSRS VALUES
!      X(NUMPT)  ,Y(NUMPT), Z(NUMPT)= X, Y, Z LOCATIONS (METERS) FOR CALCULATIONS - NUMPT VALUES
!      L (NUMSRS)=SOURCE LENGTH (METERS) IN X DIRECTION
!      B (NUMSRS)=SOURCE   WIDTH    (METERS)     IN    Y   DIRECTION
!      U = GROUND-WATER PORE VELOCITY (METERS/YEAR)
!      ALPHX = DISPERSIVITY (METERS) IN X DIRECTION
!      ALPHY = DISPERSIVITY (METERS) IN Y DIRECTION
!      ALPHZ = DISPERSIVITY (METERS) IN Z DIRECTION
!      H=AQUIFER THICKNESS (METERS) - (=99999  FOR  INFINITE  THICKNESS)
!      T=TIME OF CALCULATION (YEARS)
!      R=RETARDATION FACTOR
!      LAMB=CONTAMINANT DECAY RATE (YEARS-1)
!      POR=POROSITY
!   *****************************************************************
       REAL M(5,10) ,MS(5) ,L(5) ,LAMB
       !INTEGER NUMSRS
       DIMENSION TS(10) ,X(20), Y(20), Z(20), B(5), XS(5), YS(5), ZS(5)
       OPEN (4,FILE='HPSIN3.DAT' ,STATUS ='UNKNOWN')
       OPEN (8,FILE='HPS0UT.DAT' ,STATUS='UNKNOWN')
!
!   READ INPUT PARAMETERS - ALL  INPUT STATEMANTS  UNFORMATED
!
      READ (4,*) NUMPT,NINT,NUMST,NUMSRS
      READ (4,*) ((M(I,J), J=1 ,NUMST), I=1  ,NUMSRS)
      READ (4,*) (TS(I), I=1 ,NUMST)
      READ (4,*) (XS(I), YS(I), ZS(I), L(I), B(I), I=1,NUMSRS)
      READ (4,*) (X(I), Y(I), Z(I), I=1, NUMPT)
      READ (4,*) U, ALPHX, ALPHY, ALPHZ
      READ (4,*) H, T, R, LAMB, POR
      WRITE (8,10)
  10    FORMAT (// 15X, 'HORIZONTAL PLANE SOURCE (HPS) MODEL' /, 5X, &
     'ANALYTICAL THREE DIMENSIONAL SOLUTION FOR GROUND-WATER CONTAMINA &
      &NT TRANSPORT'//)
      WRITE (8,20)
  20    FORMAT(1X, 'INPUT VALUES', // ,3X, 'SOURCE DATA')
      DO 29 I=1, NUMSRS
      WRITE (8,22) XS(I), YS(I), ZS(I), L(I), B(I)
  22    FORMAT (/5X, 'X  LOCATION(M)' ,5X, 'Y LOCATION(M)' ,5X, 'Z  LOCATION(M)', &
     5X, 'LENGTH(M)', 5X, 'WIDTH(M)' , &
     /6X, F7.1, 12X, F7.1, 11X, F7.1, 10X, F7.1, 8X, F7.1)
      WRITE (8,24)
  24    FORMAT ("     ", 'TIME(YR)', "         ", 'SOURCE MASS FLUX(GM/YR)')
      WRITE (8,26) (TS(J), M(I,J), J=1, NUMST)
  26    FORMAT(7X, F7.1, 13X, E11.4)
  29    CONTINUE
      WRITE (8,30)  U
  30    FORMAT (/5X, 'GROUNDWATER VELOCITY IN X DIRECTION(M/YR)= ' ,F5.2)
      WRITE (8,40)  ALPHX, ALPHY ,ALPHZ
  40    FORMAT (5X, 'DISPERSIVITY VALUES(M)', 5X, 'X = ', F4.2, 5X, 'Y = ', F4.2, &
         & 5X, 'Z= ' ,F4.2)
      WRITE (8,50) H,POR
  50    FORMAT (5X, 'AQUIFER THICKNESS(M)= ' ,F7.2, 10X, 'POROSITY= ' ,F4.2)
      WRITE (8,60) R, LAMB
  60    FORMAT (5X,'RETARDATION FACTOR= ', F6.1, 5X, 'DECAY RATE(YR-1)= ' &
        ,F10. 8)
      WRITE (8, 70)  T
  70    FORMAT(///, 1X, 'PROGRAM RESULT, TIME(YEARS)= ', F7.1/"     ", &
            &'X LOCATION(M)' ,5x, 'Y LOCATION(M)' ,5X,&
            &'Z  LOCATION(M)',5X, 'CONCENTRATION(MG/L)', 4X,&
            &'INT. TIME LIMITS(YRS)' /83X,'START', 6X,'END'/)
!
! CALCULATE COEFFICIENTS FOR CALCULATION
!
       DX = ALPHX*U
       DY = ALPHY*U
       DZ = ALPHZ*U
       C1 = 1/ (R*POR)
       CX3 = U/R
       CX4 = 4 * DX/R
       CY3 = 4*DY/R
       CZ1 = 4*DZ/R
       CZ2 = 3.1416*DZ/R
       CZ3 = 1/H
       CZ4 = (3.1416**2)*DZ / ((H**2)*R)
       CZ5 = 3.1416/H
!
!  SPECIFY INITIAL OR CONTINUOUS SOURCE RATES FOR NUMSRS SOURCES
!
      DO IS = 1, NUMSRS
      MS(IS)=M(IS,1)
      END DO
!
!  BEGIN DO LOOP FOR NUMPT CALCULATION POINTS
!
      DO 999  IP=1,NUMPT
      C=0
!
!   BEGIN DO LOOP FOR INDIVIDUAL SOURCE CALCULATIONS
!
      DO 899 IS=1, NUMSRS
      XC = X(IP)-XS(IS)
      YC = Y(IP)-YS(IS)
      ZC = Z(IP)-ZS(IS)
      CX1 = 1/(2*L(IS))
      CX2 = L(IS)/2
      CY1 = 1/(2*B(IS))
      CY2 = B(IS)/2
!  CALCULATE TIME INTEGRATIONLIMITS
!      CT1 = 2*XC*U/R + l00 * DX/R
!      CT2 = (U/R)**2 + 4 * DX * LAMB/R
!      CT3 = (U/R)**2
!      CT4= (CT1/CT2)**2 - 4 *XC**2/CT2
!      IF (CT4.LE.0) THEN
!      Tl = T - CT1/CT3 - ((CT1/CT3)**2 - 4 * XC**2/CT3)**0.5
!      T2 = T + CT1/CT3 - ((CTl/CT3)**2 - 4 * XC**2/CT3)**0.5
!      ELSE
!      Tl = T - CT1/CT2 - CT4**0.5
!      T2 = T + CTI /CT2 - CT4**0.5
!      ENDIF
!      IF(Tl.LT.0) Tl=0
!      IF(T2.GT.T) T2=T
      T1 = 0
      TAU = T1
      T2 = T
      DELT = (T2 - T1) / NINT
!
!   CALCULATE    SOURCE    STRENGTH    AT   TIMESTEP    TAU
!
!     IF (NUMST.EQ.1)  GO TO 400
!      JS=0
!  200   JS=JS + 1
!  300   IF (TAU.GT.TS(JS+1)) GO TO 200
!       MS(IS) = M(IS,JS) + (M(IS,JS+1) - M(IS,JS)) * (TAU-TS(JS))/&
!     &	        (TS(JS+1) - TS(JS))
!      DELC=0.
!      IF (MS(IS).EQ.0)   GO TO 800
!
!   CALCULATE X-DIRECTION GREEN's FUNCTION
!
  400 TTAU = T - TAU
      XG1 = (CX2 + XC - CX3 * TTAU) /  (CX4*TTAU)**0.5
      ABSXG1 = ABS(XG1)
      IF (ABSXG1.LT.0.00001) THEN
        XGl=0.
      ELSEIF (XG1.LT.-6) THEN
        XG1 = -1
      ELSEIF(XG1.GT.6) THEN
         XG1 = 1
      ELSE
        XG1 = ERF(XG1)
      ENDIF
      XG2 = (CX2 - XC + CX3 *TTAU) / (CX4 *TTAU)**0.5
      ABSXG2 = ABS(XG2)
      IF(ABSXG2.LT.0.00001) THEN
        XG2=0.
      ELSEIF(XG2.LT.-6) THEN
        XG2 = -1
      ELSEIF(XG2.GT.6) THEN
        XG2=1
      ELSE
        XG2=ERF(XG2)
      ENDIF
      XG = CX1 * (XG1 + XG2)
!
!  CALCULATE Y-DIRECTION GREEN' FUNCTION
!
      YG1 = (CY2+YC)/SQRT(CY3*TTAU)
      ABSYG1 = ABS(YG1)
      IF(ABSYG1.LT.0.00001) THEN
        YGl = 0.
      ELSEIF(YG1.LT.-6) THEN
        YGl = -1
      ELSEIF (YG1.GT.6) THEN
        YGl=1
      ELSE
        YG1 = ERF(YG1)
      ENDIF
      YG2 = (CY2-YC) / (CY3*TTAU)**0.5
      ABSYG2 = ABS(YG2)
      IF(ABSYG2.LT.0.00001) THEN
         YG2 = 0.
      ELSEIF(YG2.LT.-6) THEN
         YG2 = -1
      ELSEIF(YG2.GT.6) THEN
         YG2 = 1
      ELSE
          YG2 = ERF(YG2)
      ENDIF
      YG = CY1 * (YG1 + YG2)
!
!  CALCULATE Z-DIRECTION GREEN'S FUNCTION FOR AQUIFER OF FINITE THICKNESS
!
      IF (H.EQ.99999) GO TO 600
      CZ6=0.
      NINF = 1+AINT((H/3.1416) * ((R*25) / (DZ*TTAU))**0.5)
      DO 500 J=1, NINF
      PRINT *, "TTAU = ", TTAU
      PRINT *, "NINF = ", NINF

      CZ7 = (J**2) * CZ4 * TTAU
      IF (CZ7.GT.25)THEN
        CZ7 = 0.
      ELSE
        CZ7 = EXP(-CZ7)
      ENDIF
      CZ8 = COS(J * CZ5 * Z(IP)) *COS(J*CZ5*ZS(IS))
      CZ6 = CZ6 + CZ7 * CZ8
      IF (CZ6.EQ.0) GO TO 550
  500 CONTINUE
      ZGABS = CZ3 * (1+2*CZ6)
  550 ZG = ABS(ZGABS)
      GO TO 700
!
!   CALCULATE   Z-DIRECTION   GREEN'S   FUNCTION    FOR   AQUIFER   OF   INFINITE  C THICKNESS
!
  600   CZ6 = (ZC**2) / (CZ1  * TTAU)
      IF (CZ6.GT.25) THEN
      CZ6=0
      ELSE
      CZ6 = EXP(-CZ6)
      ENDIF
      CZ7 = (Z(IP)+ZS(IS))**2 / (CZ1 * TTAU)
      IF(CZ7.GT.25)  THEN
      CZ7=0
      ELSE
      CZ7=EXP(-CZ7)
      ENDIF
      ZG = (CZ6+CZ7) / (4 * CZ2 * TTAU)**0.5
!
!  CALCULATE DELTA AND CUMULATIVE CONCENTRATION FOR TIMESTEP TAU
!
  700 TLB = LAMB * TTAU
       PRINT *,' TLB = ', TLB
       IF (TLB.GT.25) THEN
       DELC=0.
       ELSE
       DELC = MS(IS) * XG * YG * ZG * EXP(-LAMB * TTAU)
       ENDIF
       IF (TAU.EQ. T-0.1*DELT) THEN
       DELC=0.5*DELC
       ENDIF
  800 C = C + DELC
   PRINT *,' DELC = ', DELC
   PRINT *,' C = ', C
!
!   CALCULATE NEWTIMESTEP
!
    TAU = TAU + DELT
    PRINT *,' TAU = ', TAU
    IF(TAU.EQ.T) TAU = T-0.1*DELT
    IF(TAU.GT.T2) GO TO  899
    GO TO 400
   899 CONTINUE
!
!   CALCULATE   CONCENTRATION AT CALCULATION POINT AND TIME
!
   CABS = C1 * C * DELT
   900 C = ABS(CABS)
   999 WRITE (8,99) X(IP), Y(IP), Z(IP), C, T1, T2
   99 FORMAT (9X, F7.1, 8X, F7.1, 8X, F7.1, 16X, E11.5, 11X, F7.1, 4X, F7.1)
      STOP
      END

