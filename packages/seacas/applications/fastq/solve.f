C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      REAL FUNCTION SOLVE (XA, XK, X, XCEN, YCEN, R1, R2)
C***********************************************************************

C  FUNCTION SOLVE = FINDS A SOLUTION TO THE SPIRAL EQUATION
C                   GIVEN AN INTERVAL THAT CONTAINS THE SOLUTION

C***********************************************************************

      EPS = 1.E-6

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

  100 CONTINUE
      IF (ABS(R1 - R2) .LT. EPS) THEN
         SOLVE = (R1 + R2)/2.0
      ELSE
         R3 = (R1 + R2)/2.0
         F3 = SPIRAL (XA, XK, X, XCEN, YCEN, R3)

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

  110 CONTINUE
      END
