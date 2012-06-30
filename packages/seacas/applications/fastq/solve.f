C $Id: solve.f,v 1.1 1990/11/30 11:16:05 gdsjaar Exp $
C $Log: solve.f,v $
C Revision 1.1  1990/11/30 11:16:05  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]SOLVE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      REAL FUNCTION SOLVE (XA, XK, X, XCEN, YCEN, R1, R2)
C***********************************************************************
C
C  FUNCTION SOLVE = FINDS A SOLUTION TO THE SPIRAL EQUATION
C                   GIVEN AN INTERVAL THAT CONTAINS THE SOLUTION
C
C***********************************************************************
C
      EPS = 1.E-6
C
      F1 = SPIRAL (XA, XK, X, XCEN, YCEN, R1)
      IF (ABS(F1) .LT. EPS) THEN
         SOLVE = R1
         GO TO 110
      END IF
      F2 = SPIRAL (XA, XK, X, XCEN, YCEN, R2)
      IF (ABS(F2) .LT. EPS) THEN
         SOLVE = R2
         GO TO 110
      END IF
C
  100 CONTINUE
      IF (ABS(R1 - R2) .LT. EPS) THEN
         SOLVE = (R1 + R2)/2.0
      ELSE
         R3 = (R1 + R2)/2.0
         F3 = SPIRAL (XA, XK, X, XCEN, YCEN, R3)
C
         IF (ABS(F3) .LT. EPS) THEN
            SOLVE = R3
            GO TO 110
         ELSE IF (F1/F3 .LT. 0.0) THEN
            R2 = R3
            F2 = F3
         ELSE
            R1 = R3
            F1 = F3
         END IF
         GO TO 100
      END IF
C
  110 CONTINUE
      END
