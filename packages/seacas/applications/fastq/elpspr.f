C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ELPSPR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2,
     &   IPNTR3, IP3, XCEN, YCEN, THETA1, THETA2, TANG, ICCW, ICW,
     &   AVALUE, BVALUE, ERR)
C***********************************************************************

C  SUBROUTINE ELPSPR = THIS ROUTINE CALCULATES THE ELIPSE PARAMETERS

C***********************************************************************

C  VARIABLES USED:
C     TANG   = TOTAL ANGLE SCRIBED BY THE ARC
C     THETA1 = FIRST CCW ANGLE OF THE ARC
C     THETA2 = SECOND CCW ANGLE OF THE ARC
C     IPNTR1 = POINTER TO FIRST COORDINATE VALUE
C     IPNTR2 = POINTER TO SECOND COORDINATE VALUE
C     IPNTR3 = POINTER TO THIRD COORDINATE VALUE
C     IP3    = THE THIRD POINT NUMBER  (CAN BE NEGATED)

C***********************************************************************

      DIMENSION COOR (2, MP), LINKP (2, MP)

      LOGICAL ERR

      PI = ATAN2(0.0, -1.0)

      TWOPI = PI + PI
      ERR = .FALSE.

C  ELIPSE GOES FROM 1ST POINT TO 2ND IN *COUNTER-CLOCKWISE* DIRECTION.

      XCEN = COOR (1, IPNTR3)
      YCEN = COOR (2, IPNTR3)

C  CHECK TO MAKE SURE THAT THE BEGINNING AND ENDING RADIUS EXIST

      IF (( (COOR (1, IPNTR1) .EQ. XCEN).AND.
     &   (COOR (2,IPNTR1) .EQ. YCEN)).OR.
     &   ((COOR (1, IPNTR2) .EQ. XCEN).AND.
     &   (COOR (2, IPNTR2) .EQ. YCEN)))THEN
         CALL PLTFLU
         WRITE (*, 10000)ABS (KNUM)
         ERR = .TRUE.
         GOTO 100
      ENDIF
      THETA1 = ATAN2 (COOR (2, IPNTR1) - YCEN, COOR (1, IPNTR1) - XCEN)
      THETA2 = ATAN2 (COOR (2, IPNTR2) - YCEN, COOR (1, IPNTR2) - XCEN)

C  NOW CALCULATE THE MAJOR AXIS (AVALUE) AND THE MINOR AXIS (BVALUE)

      X1 = COOR (1, IPNTR1) - XCEN
      Y1 = COOR (2, IPNTR1) - YCEN
      X2 = COOR (1, IPNTR2) - XCEN
      Y2 = COOR (2, IPNTR2) - YCEN

C  CHOOSE THE APPROPRIATE ELIPSE DEFINITION

      IF (Y1 * Y1 .EQ. Y2 * Y2) THEN
         AVALUE = SQRT (X1 * X1 + Y1 * Y1)
         BVALUE = AVALUE
      ELSEIF ((Y1 .EQ. 0.) .AND. (X2 .EQ. 0.)) THEN
         AVALUE = X1
         BVALUE = Y2
      ELSEIF ((Y2 .EQ. 0.) .AND. (X1 .EQ. 0.)) THEN
         AVALUE = X2
         BVALUE = Y1
      ELSE
         RATIO = SQRT (ABS ( (X1 * X1 - X2 * X2) /
     &      (Y1 * Y1 - Y2 * Y2) ))
         IF (RATIO .EQ. 0.) THEN
            AVALUE = SQRT (X1 * X1 + Y1 * Y1)
            BVALUE = AVALUE
         ELSEIF (RATIO .EQ. 1.0) THEN
            AVALUE = SQRT (X1 * X1 + Y1 * Y1)
            BVALUE = AVALUE
         ELSE
            IF ((Y2 .EQ. 0.) .OR. ((X2 .EQ. 0.) .AND. (Y1 .NE. 0.)) )
     &         THEN
               IF (RATIO .GT. 1.) THEN
                  VX = 1.
                  VY = - (1./RATIO**2) * (X1 / Y1)
               ELSE
                  VX = 1.
                  VY = - (RATIO**2) * (X1 / Y1)
               ENDIF
               D0 = SQRT (X1 * X1 + Y1 * Y1)
               A8 = ACOS ( ((X1 * VX) + (Y1 * VY)) /
     &            (D0 * SQRT (VX * VX + VY * VY)) )
               A7 = PI - A8
               A2 = ABS (ATAN2 (Y1, X1))
               THETA1 = ABS(ATAN2 (VY, 1.))
            ELSE
               IF (RATIO .GT. 1.) THEN
                  VX = 1.
                  VY = - (1./RATIO**2) * (X2 / Y2)
               ELSE
                  VX = 1.
                  VY = - (RATIO**2) * (X2 / Y2)
               ENDIF
               D0 = SQRT (X2 * X2 + Y2 * Y2)
               A8 = ACOS ( ((X2 * VX) + (Y2 * VY)) /
     &            (D0 * SQRT (VX * VX + VY * VY)) )
               A7 = PI - A8
               A2 = ABS (ATAN2 (Y2, X2))
               THETA1 = ABS(ATAN2 (VY, 1.))
            ENDIF

            RADMAX = MAX(A7,A8)
            CALL ETHETA (A7, A8, A2, THETA1, RADMAX, THETA, ERR)
            IF (ERR) THEN
               WRITE (*, 10010) ABS (KNUM)
               GOTO 100
            ENDIF

            CVALUE = D0 * SIN (A8 - THETA) / SIN (A2 - A8 + THETA)
            BVALUE = SQRT (ABS (CVALUE **2 / (RATIO **2 - 1)) )
            AVALUE = BVALUE * RATIO
         ENDIF
      ENDIF

C  NOW GET THE ANGLES GOING THE RIGHT WAY

      THETA1 = ATAN2 (COOR (2, IPNTR1) - YCEN, COOR (1, IPNTR1) - XCEN)
      THETA2 = ATAN2 (COOR (2, IPNTR2) - YCEN, COOR (1, IPNTR2) - XCEN)
      IF (IPNTR1 .EQ. IPNTR2) THEN
         THETA2 = THETA1 + TWOPI
      ELSEIF ( (IP3 .GE. 0) .AND. (THETA2 .LE. THETA1) ) THEN
         THETA2 = THETA2+TWOPI
      ELSEIF ( (IP3 .LT. 0) .AND. (THETA1 .LE. THETA2) ) THEN
         THETA1 = THETA1+TWOPI
      ENDIF
      TANG = THETA2 - THETA1
      IF (IP3 .LT. 0) THEN
         ICCW = IPNTR2
         ICW = IPNTR1
      ELSE
         ICCW = IPNTR1
         ICW = IPNTR2
      ENDIF

  100 CONTINUE

      RETURN

10000 FORMAT (' CENTER POINT FOR LINE', I5, ' LIES ON ONE OF',
     &   ' THE ENDPOINTS')
10010 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE AN ELIPSE')

      END
