C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      REAL FUNCTION ESOLVE (A7, A8, A2, ANG1, ANG2)
C***********************************************************************

C  FUNCTION ESOLVE = FINDS A SOLUTION TO THE ELIPSE EQUATION
C                    GIVEN AN INTERVAL THAT CONTAINS THE SOLUTION

C***********************************************************************

      EPS = 1.E-6

      F1 = ELIPSE (A7, A8, A2, ANG1)
      IF (ABS(F1) .LT. EPS) THEN
         ESOLVE = ANG1
         GO TO 110
      END IF
      F2 = ELIPSE (A7, A8, A2, ANG2)
      IF (ABS(F2) .LT. EPS) THEN
         ESOLVE = ANG2
         GO TO 110
      END IF

  100 CONTINUE
      IF (ABS(ANG1 - ANG2) .LT. EPS) THEN
         ESOLVE = (ANG1 + ANG2)/2.0
      ELSE
         ANG3 = (ANG1 + ANG2)/2.0
         F3 = ELIPSE (A7, A8, A2, ANG3)

         IF (ABS(F3) .LT. EPS) THEN
            ESOLVE = ANG3
            GO TO 110
         ELSE IF (F1/F3 .LT. 0.0) THEN
            ANG2 = ANG3
            F2 = F3
         ELSE
            ANG1 = ANG3
            F1 = F3
         END IF
         GO TO 100
      END IF

  110 CONTINUE

      RETURN

      END
