C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ETHETA (A7, A8, A2, THETA1, THETA2, THETA, ERR)
C***********************************************************************

C  SUBROUTINE ETHETA = ITERATIVELY SOLVES THE ELIPTICAL PROBLEM OF
C                      FINDING AN "A" DISTANCE GIVEN TWO POINTS ON
C                      THE ELIPSE AND A CENTER POINT

C***********************************************************************

      LOGICAL ERR

C  START WITH 20 INCREMENTS, EACH PASS INCREMENTS DECREASE TEN FOLD

      ANGINC = (THETA2 - THETA1) * .05
      ANG = THETA1
      F1 = ELIPSE (A7, A8, A2, ANG)
      IF (F1 .EQ. 0.0) THEN
         THETA = ANG
         GO TO 110
      END IF
      ANG2 = ANG + ANGINC
  100 CONTINUE
      F2 = ELIPSE (A7, A8, A2, ANG2)
      IF (F2 .EQ. 0.0) THEN
         THETA = ANG2
      ELSE IF (F1*F2 .LT. 0.0) THEN
         THETA = ESOLVE (A7, A8, A2, ANG, ANG2)
      ELSE
         ANG = ANG2
         ANG2 = ANG2 + ANGINC
         IF (ANG2 .LE. THETA2) GO TO 100
         ERR = .TRUE.
         GO TO 120
      END IF

  110 CONTINUE

  120 CONTINUE

C  FIND THE SECOND ROOT IF THE FIRST ONE HAS BEEN LOCATED

      IF(.NOT.ERR)THEN
         ANG=THETA+ANGINC
         F1 = ELIPSE (A7, A8, A2, ANG)
         F2 = ELIPSE (A7, A8, A2, THETA2)
         IF (F1 .EQ. 0.0) THEN
            THETA = ANG
         ELSEIF (F2 .EQ. 0.0) THEN
            THETA = THETA2
         ELSE IF (F1*F2 .LT. 0.0) THEN
            THETA = ESOLVE (A7, A8, A2, ANG, THETA2)
         ELSE
            GO TO 130
         END IF
      END IF

  130 CONTINUE
      RETURN

      END
