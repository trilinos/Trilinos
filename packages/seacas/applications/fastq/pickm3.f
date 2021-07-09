C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PICKM3 (N, X, Y, ANGLE, M1, M2, IFIRST)
C***********************************************************************

C  SUBROUTINE PICKM3 = DETERMINES A REASONABLE SHAPE FOR A LOGICAL
C                      TRIANGLE WITH PERIMETER GIVEN IN X AND Y

C***********************************************************************

      PARAMETER  (RLARGE = 100000.)
      DIMENSION X (N), Y (N), ANGLE (N)
      DIMENSION SMANG (7), INDEX (7)
      DIMENSION ISORT (3)

C  FORM THE LIST OF SMALLEST ANGLES

      NSA = 6
      DO 100 I = 1, NSA
         SMANG (I) = 10.
         INDEX (I) = 0
  100 CONTINUE

      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      AGOLD = ATAN2 (Y (1) - Y (N), X (1) - X (N))

      DO 130 J = 1, N

C  GET THE ANGLE FORMED BY THIS SET OF POINTS

         NEXT = J + 1
         IF  (NEXT .GT. N) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J) ,  X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI)DIFF = DIFF - TWOPI
         IF (DIFF .LT.  - PI)DIFF = DIFF + TWOPI
         ANGLE (J) = PI - DIFF
         AGOLD = AGNEW

C  SORT THIS ANGLE AGAINST PREVIOUS ANGLES TO SEE IF IT IS ONE OF
C  THE SMALLEST

         SMANG (NSA + 1) = ANGLE (J)
         INDEX (NSA + 1) = J
         DO II = 1, NSA
            I = NSA + 1 - II
            IF  (SMANG (I + 1) .GE. SMANG (I)) GO TO 120
            TEMP = SMANG (I)
            ITEMP = INDEX (I)
            SMANG (I) = SMANG (I + 1)
            INDEX (I) = INDEX (I + 1)
            SMANG (I + 1) = TEMP
            INDEX (I + 1) = ITEMP
         end do
  120    CONTINUE

  130 CONTINUE

C  DETERMINE OPTIMUM ORIGIN / SHAPE COMBINATION FOR A TRIANGLE

      ATOL = PI * 150. / 180.

C  FIND SIDE DIVISION USING 5 SMALLEST ANGLES AND CHECK CONDITION

      DO 140 I = 1,  3
         ISORT (I) = INDEX (I)
  140 CONTINUE
      DO 160 I = 1,  2
         DO 150 J = I + 1,  3
            IF (ISORT (I) .GT. ISORT (J)) THEN
               ITMP = ISORT (I)
               ISORT (I) = ISORT (J)
               ISORT (J) = ITMP
            ENDIF
  150    CONTINUE
  160 CONTINUE
      I1 = ISORT (1)
      I2 = ISORT (2)
      I3 = ISORT (3)
      MM1 = I2  -  I1
      IF (MM1 .LT. 0) MM1 = N  +  MM1
      MM2 = I3  -  I2
      IF (MM2 .LT. 0) MM2 = N  +  MM2
      MM3 = N  -  MM1  -  MM2
      MAX = MAX0 (MM1,  MM2,  MM3)
      IF (MAX .LE. N - MAX - 2) THEN

C  ADD UP ASSIGNED ANGLES

         IFIRST = I1
         M1 = MM1
         M2 = MM2
         GBEST = ANGLE (I1) + ANGLE (I2) + ANGLE (I3)
      ELSE
         IFIRST = 1
         GBEST = RLARGE
      END IF

C  LIMIT THE SIZE OF ANY SIDE

      MMAX = (N - 2) / 2

C  GO AROUND THE PERIMETER USING THE 10 SMALLEST ANGLES AS POSSIBLE
C  STARTING POINTS,  AND THEN FIND THE BEST COMBINATION OF SIDE LENGTHS

      DO 190 ISA = 1, NSA
         IF (SMANG (ISA) .LE. ATOL) THEN
            I1 = INDEX (ISA)
            SUM1 = ANGLE (I1)
            IF (SUM1  .GE.  GBEST) GO TO 190

C  ASSIGN A TRIAL SECOND NODE

            DO 180 N1 = 2, MMAX
               I2 = I1 + N1
               IF (I2 .GT. N)I2 = I2 - N
               SUM2 = SUM1 + ANGLE (I2)
               IF (SUM2  .GE.  GBEST) GO TO 180

C  ASSIGN A TRIAL THIRD NODE

               DO 170 N2 = 2, N - N1 - 2
                  I3 = I2 + N2
                  IF (I3 .GT. N)I3 = I3 - N
                  GVAL = SUM2 + ANGLE (I3)
                  IF (GVAL  .GE.  GBEST) GO TO 170

C  FIND SIDE DIVISION AND CHECK CONDITION

                  MM1 = I2  -  I1
                  IF  (MM1 .LT. 0) MM1 = N + MM1
                  MM2 = I3  -  I2
                  IF  (MM2 .LT. 0) MM2 = N + MM2
C ... Guess by GDS, MM1 substituted for MM?
                  MM3 = N - MM1  - MM2
                  MAX = MAX0 (MM1, MM2, MM3)
                  IF  (MAX .LE. N - MAX - 2) THEN

C  ADD UP ASSIGNED ANGLES AND COMPARE TO PREVIOUS TRIALS

                     IFIRST = I1
                     M1 = MM1
                     M2 = MM2
                     GBEST = GVAL
                  ENDIF
  170          CONTINUE
  180       CONTINUE
         ENDIF
  190 CONTINUE

      RETURN

      END
