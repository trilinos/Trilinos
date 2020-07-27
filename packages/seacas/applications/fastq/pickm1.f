C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

c
C
C     C* FILE: [.QMESH]PICKM1.FOR
C     C* MODIFIED BY: TED BLACKER
C     C* MODIFICATION DATE: 7/6/90
C     C* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PICKM1 (N, X, Y, ANGLE, M1, IFIRST, REAL)
C***********************************************************************
C
C     SUBROUTINE PICKM1 = DETERMINES A REASONABLE SHAPE FOR A LOGICAL
C     RECTANGLE WITH PERIMETER GIVEN IN X AND Y
C
C***********************************************************************
C
      DIMENSION X (N), Y (N), ANGLE (N)
      DIMENSION SMANG (7), INDEX (7)
C
      LOGICAL REAL
C
C     FORM THE LIST OF SMALLEST ANGLES
C
      NSA = 6
      DO 100 I = 1, NSA
         SMANG (I) = 10.
         INDEX (I) = 0
 100  CONTINUE
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
      AGOLD = ATAN2 (Y (1) - Y (N), X (1) - X (N))
C
      DO 130 J = 1, N
C
C     GET THE ANGLE FORMED BY THIS SET OF POINTS
C
         NEXT = J + 1
         IF  (NEXT .GT. N) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J) ,  X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI)DIFF = DIFF - TWOPI
         IF (DIFF .LT.  - PI)DIFF = DIFF + TWOPI
         ANGLE (J) = PI - DIFF
         AGOLD = AGNEW
C
C     SORT THIS ANGLE AGAINST PREVIOUS ANGLES TO SEE IF IT IS ONE OF
C     THE SMALLEST
C
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
 120     CONTINUE
C
 130  CONTINUE
C
C     DETERMINE OPTIMUM ORIGIN / SHAPE COMBINATION
C
      ATOL = PI * 150. / 180.
      IFIRST = 1
      M1 = N / 4
      M2 = N / 2 - M1
      I2 = 1 + M1
      I3 = I2 + M2
      I4 = I3 + M1
      GBEST = ANGLE (1) + ANGLE (I2) + ANGLE (I3) + ANGLE (I4)
      BADANG = AMAX1 (ANGLE (1), ANGLE (I2), ANGLE (I3), ANGLE (I4))
C
      MMAX = N / 2 - 1
      AMAXEL = DBLE(N / 4) * DBLE( (N + 2) / 4)
      DO 150 ISA = 1, NSA
         IF (SMANG (ISA) .LE. ATOL) THEN
            I1 = INDEX (ISA)
            DO 140 M = 1, MMAX
               M2 = N / 2 - M
               I2 = I1 + M
               IF  (I2 .GT. N) I2 = I2 - N
               I3 = I2 + M2
               IF  (I3 .GT. N) I3 = I3 - N
               I4 = I3 + M
               IF  (I4 .GT. N) I4 = I4 - N
               AFAC = ANGLE (I1) + ANGLE (I2) + ANGLE (I3) + ANGLE (I4)
               ERAT = AMIN1 (AMAXEL / DBLE(M * M2) ,  5.)
               EFAC =  (ERAT + 15.) / 16.
               GVAL = AFAC * EFAC
               IF (GVAL .LT. GBEST) THEN
                  BADANG = AMAX1 (ANGLE (I1), ANGLE (I2), ANGLE (I3),
     &                 ANGLE (I4))
                  IFIRST = I1
                  M1 = M
                  GBEST = GVAL
               ENDIF
 140        CONTINUE
         ENDIF
 150  CONTINUE
      IF ( (REAL) .AND. (BADANG .GT. 2.62)) THEN
         CALL MESAGE (' **  WARNING: CORNER (S) OF THE REGION HAVE  **')
         CALL MESAGE (' **           LARGE ANGLES  (> 150 DEGREES.) **')
         CALL MESAGE (' **           POORLY FORMED MESH MAY RESULT  **')
      ENDIF
C
      RETURN
      END
