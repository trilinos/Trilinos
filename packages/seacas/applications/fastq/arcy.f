C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ARCY (XCEN, YCEN, THETA1, THETA2, XK, XA, X, YTRY, ERR)
C***********************************************************************

C  SUBROUTINE ARCY = ITERATIVELY SOLVES THE LOGARITHMIC SPIRAL PROBLEM
C                    TO DETERMINE A Y VALUE GIVEN AN X THAT INTERSECTS
C                    THE ARC

C***********************************************************************

      LOGICAL ERR

C  START WITH 10 INCREMENTS, EACH PASS INCREMENTS DECREASE TEN FOLD

      ANGINC = (THETA2 - THETA1) * .05
      ANG = THETA1
      F1 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG)
      IF (F1 .EQ. 0.0) THEN
         THETA = ANG
         GO TO 110
      END IF
      ANG2 = ANG + ANGINC
  100 CONTINUE
      F2 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG2)
      IF (F2 .EQ. 0.0) THEN
         THETA = ANG2
      ELSE IF (F1*F2 .LT. 0.0) THEN
         THETA = SOLVE(XA, XK, X, XCEN, YCEN, ANG, ANG2)
      ELSE
         ANG = ANG2
         ANG2 = ANG2 + ANGINC
         IF (ANG2 .LE. THETA2) GO TO 100
         ERR = .TRUE.
         GO TO 120
      END IF

  110 CONTINUE
      YTRY = (XA * EXP(XK * THETA)) * SIN(THETA) + YCEN

  120 CONTINUE

C  FIND THE SECOND ROOT IF THE FIRST ONE HAS BEEN LOCATED

      IF(.NOT.ERR)THEN
         ANG=THETA+ANGINC
         F1 = SPIRAL (XA, XK, X, XCEN, YCEN, ANG)
         F2 = SPIRAL (XA, XK, X, XCEN, YCEN, THETA2)
         IF (F1 .EQ. 0.0) THEN
            THETA = ANG
         ELSEIF (F2 .EQ. 0.0) THEN
            THETA = THETA2
         ELSE IF (F1*F2 .LT. 0.0) THEN
            THETA = SOLVE(XA, XK, X, XCEN, YCEN, ANG, THETA2)
         ELSE
            GO TO 130
         END IF
      END IF

      YTRY2 = (XA * EXP(XK * THETA)) * SIN(THETA) + YCEN
      YTRY = MAX(YTRY,YTRY2)
  130 CONTINUE
      RETURN

      END
