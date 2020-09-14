C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INVMAP (X0, Y0, X1, Y1, X2, Y2, X3, Y3, X4, Y4, SXI,
     &   SETA, INSIDE)
C***********************************************************************

C  THIS IS A TEST OF THE INVERTED MAPPING OF AN ELEMENT

C***********************************************************************

      DOUBLE PRECISION AX, BX, CX, DX, AY, BY, CY, DY
      DOUBLE PRECISION ALPHA, BETA, GAMMA, RAD
      DOUBLE PRECISION XI, ETA, XI1, ETA1, XI2, ETA2

      LOGICAL INSIDE

      EPS = 1.E-3
      EPS2 = 1.E-10

C  GET THE A, B, C, AND D VALUES FOR X AND Y.

      AX = X1 - X0
      BX = X2 - X1
      CX = X1 - X2 + X3 -X4
      DX = X4 - X1
      AY = Y1 - Y0
      BY = Y2 - Y1
      CY = Y1 - Y2 + Y3 -Y4
      DY = Y4 - Y1

C  CALCULATE THE ALPHA, BETA, AND GAMMA VALUES.

      ALPHA = (CY * DX) - (CX * DY)
      BETA  = (AX * CY) - (AY * CX) + (BY * DX) - (BX * DY)
      GAMMA = (AX * BY) - (AY * BX)

C  CALCULATE THE XI AND ETA VALUES.

      IF (ALPHA .EQ. 0.) THEN
         ETA = -GAMMA / BETA
         IF ((ETA .EQ. 0) .AND. (BX .EQ. 0)) THEN
            XI = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA .EQ. 1.)) THEN
            XI = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA)) THEN
            XI = -1000.
         ELSE
            XI = (- AX - (DX * ETA)) / (BX + (CX * ETA))
         ENDIF
      ELSE
         RAD = BETA**2 - (4. * ALPHA * GAMMA)
         IF (RAD .LT. 0.) THEN

C**               NEGATIVE RADICAL PROBLEM AS IT APPEARS THAT
C**               THIS MAY OCCUR - IT JUST MEANS THAT THE POINT
C**               TRULY IS NOT IN THE ELEMENT.

C            CALL MESAGE ('** ERROR - NEGATIVE RADICAL IN INVMAP **')
            INSIDE = .FALSE.
            GOTO 100
         ENDIF
         RAD = DSQRT (RAD)
         ETA1 = (- BETA + RAD) / (2. * ALPHA)
         ETA2 = (- BETA - RAD) / (2. * ALPHA)

         IF ((ABS(ETA1) .LT. EPS2) .AND. (ABS(BX) .LT. EPS2)) THEN
            XI1 = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA1 .EQ. 1.)) THEN
            XI1 = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA1)) THEN
            XI1 = -1000.
         ELSE
            XI1 = (- AX - (DX * ETA1)) / (BX + (CX * ETA1))
         ENDIF

         IF ((ABS(ETA2) .LT. EPS2) .AND. (ABS(BX) .LT. EPS2)) THEN
            XI2 = (Y0 - Y1) / (Y2 - Y1)
         ELSE IF ((BX .EQ. -CX) .AND. (ETA2 .EQ. 1.)) THEN
            XI2 = (Y0 - Y3)/(Y4 - Y3)
         ELSE IF (BX .EQ. (-CX * ETA2)) THEN
            XI2 = -1000.
         ELSE
            XI2 = (- AX - (DX * ETA2)) / (BX + (CX * ETA2))
         ENDIF

         D1 = DSQRT (ETA1*ETA1 + XI1*XI1)
         D2 = DSQRT (ETA2*ETA2 + XI2*XI2)
         IF (D1 .LT. D2) THEN
            ETA = ETA1
            XI = XI1
         ELSE
            ETA = ETA2
            XI = XI2
         ENDIF
      END IF

C  CHECK TO SEE IF ETA AND XI ARE WITHIN THE ELEMENT

      IF (.NOT. ((ETA .LE. 1.0 + EPS) .AND.
     &   (ETA .GE. 0.0 - EPS)) ) THEN
         INSIDE = .FALSE.
         GOTO 100
      ELSE IF (.NOT. ((XI .LE. 1.0 + EPS) .AND.
     &   (XI .GE. 0.0 - EPS)) ) THEN
         INSIDE = .FALSE.
         GOTO 100
      ELSE
         INSIDE = .TRUE.
      ENDIF
      SXI = XI
      SETA = ETA

  100 CONTINUE
      RETURN

      END
