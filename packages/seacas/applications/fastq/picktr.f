C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PICKTR (NPER, X, Y, NID, ANGLE, HALFC, I1, I2, I3, I4,
     &   I5, I6, I7, I8)
C***********************************************************************

C  SUBROUTINE PICKTR = DETERMINES A REASONABLE SHAPE FOR A BACK-TO-BACK
C                      SET OF TRIANGLES (TRANSITION REGION)

C***********************************************************************

      PARAMETER (RLARGE = 1000000.)
      DIMENSION X(NPER), Y(NPER), NID(NPER), ANGLE(NPER)
      DIMENSION SMANG(7), INDEX(7)
      DIMENSION ISORT(4)

      LOGICAL HALFC

      PI = ATAN2(0.0, -1.0)
      PID2 = 0.5 * PI
      TWOPI = 2.0 * PI

C  FORM THE LIST OF SMALLEST ANGLES

      NSA = 6
      DO 100 I = 1, NSA
         SMANG(I) = 10.
         INDEX(I) = 0
  100 CONTINUE

      AGOLD = ATAN2 (Y (1) - Y (NPER), X (1) - X (NPER))

      DO 130 J = 1, NPER

C  GET THE ANGLE FORMED BY THIS SET OF POINTS

         NEXT = J + 1
         IF (NEXT .GT. NPER) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J), X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI) DIFF = DIFF - TWOPI
         IF (DIFF .LT. - PI) DIFF = DIFF + TWOPI
         ANGLE (J) = PI - DIFF
         AGOLD = AGNEW

C  SORT THIS ANGLE AGAINST PREVIOUS ANGLES TO SEE IF IT IS ONE OF
C  THE SMALLEST

         SMANG (NSA + 1) = ANGLE (J)
         INDEX (NSA + 1) = J
         DO 110 II = 1, NSA
            I = NSA + 1 - II
            IF  (SMANG (I + 1) .GE. SMANG (I)) GO TO 120
            TEMP = SMANG (I)
            ITEMP = INDEX (I)
            SMANG (I) = SMANG (I + 1)
            INDEX (I) = INDEX (I + 1)
            SMANG (I + 1) = TEMP
            INDEX (I + 1) = ITEMP
  110    CONTINUE
  120    CONTINUE

  130 CONTINUE

C  DETERMINE TWO/FOUR BEST CORNER POINTS FOR SEMICIRCLE/TRANSITION REGION

      ATOL = PI * 150. / 180.

C  FIND SIDE DIVISION USING 4 SMALLEST ANGLES AND CHECK CONDITION

      DO 140 I = 1, 4
         ISORT (I) = INDEX (I)
  140 CONTINUE
      DO 160 I = 1, 3
         DO 150 J = I + 1, 4
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
      I4 = ISORT (4)
      M1 = I2 - I1
      IF (M1 .LT. 0) M1 = NPER + M1
      M2 = I3 - I2
      IF (M2 .LT. 0) M2 = NPER + M2
      M3 = I4 - I3
      IF (M3 .LT. 0) M3 = NPER + M3
      M4 = NPER - M1 - M2 - M3

C  USE THE LONGEST SIDE THAT DOES NOT HAVE OPPOSITE
C  MATCHES AS THE CHOICE FOR THE BASE (OF TRANSITIONS)
C  THE BASE MUST BE AT LEAST 4 INTERVALS LONG

      IF ( (M1 .EQ. M3) .AND. (.NOT. HALFC)) THEN
         MMAX = MAX0 (M2, M4)
         IF (MMAX .GE. 4) THEN
            IF (M2 .EQ. MMAX) THEN
               IFIRST = I2
               MBASE = M2
            ELSE
               IFIRST = I4
               MBASE = M4
            ENDIF
         ENDIF
      ELSEIF ( (M2 .EQ. M4) .AND. (.NOT. HALFC)) THEN
         MMAX = MAX0 (M1, M3)
         IF (MMAX .GE. 4) THEN
            IF (M1 .EQ. MMAX) THEN
               IFIRST = I1
               MBASE = M1
            ELSE
               IFIRST = I3
               MBASE = M3
            ENDIF
         ENDIF
      ELSE
         MMAX = MAX0 (M1, M2, M3, M4)
         IF (MMAX .GE. 4) THEN
            IF (M1 .EQ. MMAX) THEN
               IFIRST = I1
               MBASE = M1
            ELSEIF (M2 .EQ. MMAX) THEN
               IFIRST = I2
               MBASE = M2
            ELSEIF (M3 .EQ. MMAX) THEN
               IFIRST = I3
               MBASE = M3
            ELSEIF (M4 .EQ. MMAX) THEN
               IFIRST = I4
               MBASE = M4
            ENDIF
         ENDIF
      ENDIF
      IF (MMAX .GE. 4) THEN
         IF (HALFC) THEN
            GBEST = ANGLE (I1) + ANGLE (I2) + ABS (PI - ANGLE (I3))
     &         + ABS (PI - ANGLE (I4))
         ELSE
            GBEST = ANGLE (I1) + ANGLE (I2) + ANGLE (I3) + ANGLE (I4)
         ENDIF
      ELSE
         IFIRST = 1
         GBEST = RLARGE
      END IF

C  GO AROUND THE PERIMETER USING THE 6 SMALLEST ANGLES AS POSSIBLE
C  STARTING POINTS, AND THEN FIND THE BEST COMBINATION OF SIDE LENGTHS

      DO 200 ISA = 1, NSA
         IF (SMANG (ISA) .LE. ATOL) THEN
            I1 = INDEX (ISA)
            SUM1 = ANGLE (I1)
            IF (HALFC) THEN
               IF (SUM1 .GT. GBEST) GO TO 200
            ELSE
               IF (SUM1 .GE. GBEST) GO TO 200
            ENDIF

C  ASSIGN A TRIAL SECOND NODE

            DO 190 N1 = 1, NPER - 4
               I2 = I1 + N1
               IF (I2 .GT. NPER) I2 = I2 - NPER
               SUM2 = SUM1 + ANGLE (I2)
               IF (HALFC) THEN
                  IF (SUM2 .GT. GBEST) GO TO 190
               ELSE
                  IF (SUM2 .GE. GBEST) GO TO 190
               ENDIF

C  ASSIGN A TRIAL THIRD NODE

               DO 180 N2 = 1, NPER - N1 - 3
                  I3 = I2 + N2
                  IF (I3 .GT. NPER) I3 = I3 - NPER
                  IF (HALFC) THEN
                     SUM3 = SUM2 + ABS (PI - ANGLE (I3))
                  ELSE
                     SUM3 = SUM2 + ANGLE (I3)
                  END IF
                  IF (HALFC) THEN
                     IF (SUM3 .GT. GBEST) GO TO 180
                  ELSE
                     IF (SUM3 .GE. GBEST) GO TO 180
                  ENDIF

C  ASSIGN A TRIAL FOURTH NODE

                  DO 170 N3 = 1, NPER - N1 - N2 - 2
                     I4 = I3 + N3
                     IF (I4 .GT. NPER) I4 = I4 - NPER
                     IF (HALFC) THEN
                        GVAL = SUM3 + ABS (PI - ANGLE (I4))
                     ELSE
                        GVAL = SUM3 + ANGLE (I4)
                     END IF
                     IF (HALFC) THEN
                        IF (GVAL .GT. GBEST) GO TO 170
                     ELSE
                        IF (GVAL .GE. GBEST) GO TO 170
                     ENDIF

C  FIND SIDE DIVISION AND CHECK CONDITION

                     M1 = I2 - I1
                     IF (M1 .LT. 0) M1 = NPER + M1
                     M2 = I3 - I2
                     IF (M2 .LT. 0) M2 = NPER + M2
                     M3 = I4 - I3
                     IF (M3 .LT. 0) M3 = NPER + M3
                     M4 = NPER - M1 - M2 - M3

C  USE THE LONGEST SIDE THAT DOES NOT HAVE OPPOSITE
C  MATCHES AS THE CHOICE FOR THE BASE (OF TRANSITIONS)
C  THE BASE MUST BE AT LEAST 4 INTERVALS LONG

                     IF ( (M1 .EQ. M3) .AND. (.NOT. HALFC)) THEN
                        MMAX = MAX0 (M2, M4)
                        IF (MMAX .GE. 4) THEN
                           IF (M2 .EQ. MMAX) THEN
                              IFIRST = I2
                              MBASE = M2
                           ELSE
                              IFIRST = I4
                              MBASE = M4
                           ENDIF
                        ENDIF
                     ELSEIF ( (M2 .EQ. M4) .AND. (.NOT. HALFC)) THEN
                        MMAX = MAX0 (M1, M3)
                        IF (MMAX .GE. 4) THEN
                           IF (M1 .EQ. MMAX) THEN
                              IFIRST = I1
                              MBASE = M1
                           ELSE
                              IFIRST = I3
                              MBASE = M3
                           ENDIF
                        ENDIF
                     ELSE
                        MMAX = MAX0 (M1, M2, M3, M4)
                        IF (MMAX .GE. 4) THEN
                           IF (M1 .EQ. MMAX) THEN
                              IFIRST = I1
                              MBASE = M1
                           ELSEIF (M2 .EQ. MMAX) THEN
                              IFIRST = I2
                              MBASE = M2
                           ELSEIF (M3 .EQ. MMAX) THEN
                              IFIRST = I3
                              MBASE = M3
                           ELSEIF (M4 .EQ. MMAX) THEN
                              IFIRST = I4
                              MBASE = M4
                           ENDIF
                        ENDIF
                        IF (MMAX .GE. 4)GBEST = GVAL
                     ENDIF
  170             CONTINUE
  180          CONTINUE
  190       CONTINUE
         ENDIF
  200 CONTINUE

C  ROTATE THE PERIMETER AND THE ANGLES SO THE BASE LEADS THE LIST

      IF (IFIRST .NE. 1) CALL FQ_ROTATE (NPER, X, Y, NID, IFIRST)
      DO 220 I = 1, IFIRST - 1
         AHOLD = ANGLE (1)
         DO 210 J = 1, NPER - 1
            ANGLE (J) = ANGLE (J + 1)
  210    CONTINUE
         ANGLE (NPER) = AHOLD
  220 CONTINUE

C  DECIDE THE TRIANGLE CORNERS

      GBEST = RLARGE

C  PICK AN ARBITRARY BASE CENTER (I3)

      DO 250 I = 3, MBASE - 1

C  FOR THIS BASE CENTER, PICK AN ARBITRARY I2 LOCATION

         DO 240 J = 2, I - 1

C  FOR THIS COMBINATION OF I3 AND I2, PICK AN ARBITRARY I4 LOCATION

            DO 230 K = I + 1, MBASE

C  CALCULATE I6 AND I8 AND ADD ANGLES TO FIND MINIMUM SUM

               KN = MBASE + 1 - K
               KK = I - J
               KL = J - 1
               KM = K - I
               MLEFT = NPER - MBASE
               KO = (MLEFT - KN + KL - KK - KM) / 2
               KP = KN + KO - KL

C  PROTECT AGAINST THE IMPOSSIBLE LENGTH PROBLEMS
C  AND THE ODD NUMBER IN THE PERIMETER INPUT ERRORS

               IF ( (KO .GT. 0) .AND. (KP .GT. 0)) THEN
                  IF (KP + KL .EQ. KN + KO) THEN

C  NOW GET THE END POINTS GIVEN THESE SIDE LENGTHS

                     J6 = MBASE + 1 + KO
                     J7 = MBASE + 1 + KO + KM
                     J8 = MBASE + 1 + KO + KM + KK

C  GET THE BASE ANGLE OF THE DIVIDER LINE

                     THETA1 = ATAN2 (Y (I + 1) - Y (I),
     &                  X (I + 1) - X (I))
                     THETA2 = ATAN2 (Y (J7) - Y (I), X (J7) - X (I))
                     THETAB = ABS (THETA2 - THETA1)
                     IF (THETAB .GT. PI) THETAB = THETAB - PI
                     IF (THETAB .LT. PID2) THETAB = PI - THETAB

C  GET THE TOP ANGLE OF THE DIVIDER LINE

                     THETA1 = ATAN2 (Y (J7 + 1) - Y (J7),
     &                  X (J7 + 1) - X (J7))
                     THETAT = ABS (THETA2 - THETA1)
                     IF (THETAT .GT. PI) THETAT = THETAT - PI
                     IF (THETAT .LT. PID2) THETAT = PI - THETAT

C  ADD THESE TO GET THE VALUE OF GVAL

                     IF (HALFC) THEN
                        GVAL = THETAB + THETAT + ABS (PI - ANGLE (J6))
     &                     + ABS (PI - ANGLE (J8)) +
     &                     (.1 * MAX0 (ABS (KK - KL),
     &                     ABS (KM - KN)) / NPER)
                     ELSE
                        GVAL = ANGLE (J6) + ANGLE (J8) + THETAB + THETAT
                     ENDIF
                     IF (GVAL .LT. GBEST) THEN
                        GBEST = GVAL
                        I1 = 1
                        I2 = J
                        I3 = I
                        I4 = K
                        I5 = MBASE + 1
                        I6 = I5 + KO
                        I7 = I6 + KM
                        I8 = I7 + KK
                     ENDIF
                  ELSE
                     CALL MESAGE ('ODD PERIMETER PROBLEMS')
                  ENDIF
               ENDIF
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
      RETURN
      END
