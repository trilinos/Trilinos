C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

c
C
C     C* FILE: [.QMESH]CSMANG.FOR
C     C* MODIFIED BY: TED BLACKER
C     C* MODIFICATION DATE: 7/6/90
C     C* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CSMANG (N, X, Y, ANGLE, NSANG, SMANG, INDEX)
C***********************************************************************
C
C     SUBROUTINE CSMANG = CALCULATES THE "NSANG" SMALLEST ANGLES
C     AND PLACES THEM IN THE SMANG ARRAY WITH
C     THE INDICES BEING PLACED IN THE INDEX ARRAY
C
C***********************************************************************
C
C     OBSERVATION: - IT DOES NOT MATTER THE ANGLE ORIENTATION.
C     - THE ANGLES ARE STORE IN THE ANGLE ARRAY AS THEY
C     APPEAR.
C     - THE SMALLEST ANGLES ARE IN ASCENDING ORDER.
C     - THE INDEX ARRAY RETURNS THE SMALLEST ANGLE POSITION
C     IN ASCENDING ORDER.
C
C
C     MODIFIED BY : HORACIO RECALDE             DATE:JAN 1988
C***********************************************************************
C
      DIMENSION X(N), Y(N), ANGLE(N)
      DIMENSION SMANG (NSANG + 1), INDEX (NSANG + 1)

      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI
C
C     FORM THE LIST OF SMALLEST ANGLES
C
      NSA = NSANG
      DO I = 1,NSA
         SMANG(I) = 10.
         INDEX(I) = 0
      end do
C
      AGOLD = ATAN2 (Y (1) - Y(N), X (1) - X (N))
C
      DO J = 1, N
C
C     GET THE ANGLE FORMED BY THIS SET OF POINTS
C
         NEXT = J + 1
         IF (NEXT .GT. N) NEXT = 1
         AGNEW = ATAN2 (Y (NEXT) - Y (J), X (NEXT) - X (J))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI)  DIFF = DIFF - TWOPI
         IF (DIFF .LT. -PI) DIFF = DIFF + TWOPI
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
            IF (SMANG (I + 1) .GE. SMANG (I)) GO TO 120
            TEMP = SMANG(I)
            ITEMP = INDEX(I)
            SMANG (I) = SMANG (I + 1)
            INDEX (I) = INDEX (I + 1)
            SMANG (I + 1) = TEMP
            INDEX (I + 1) = ITEMP
         end do
 120     CONTINUE
C
      end do
C
      RETURN
C
      END
