C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TRIDEL (MXND, MLN, XN, YN, ZN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, ANGLE, LNODES, BNSIZE,
     &   NLOOP, DEV1, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, GRAPH,
     &   VIDEO, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE TRIDEL = CHECKS FOR ANY TRIANGULAR SHAPED QUADS ATTACHED
C                      TO A THREE NODE ELEMENT AND DELETES THEM WHEN
C                      FOUND AND POSSIBLE

C***********************************************************************

      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
      DIMENSION NODES(4), K(3)
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), ZN(MXND), NUID(MXND)

      CHARACTER*3 DEV1
      LOGICAL ERR, DONE, GRAPH, CHECK, REDO, CCW, VIDEO, PASSED, NOROOM

      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI

      ERR = .FALSE.
      DONE = .FALSE.
      CHECK = .TRUE.
      CCW = .TRUE.
      KMAX = 30
      KOUNT = 0
      KKKADD = 0

  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. KMAX) GOTO 180
  110 CONTINUE
      REDO = .FALSE.

      DO 120 I = 1, NNN
         IF ((LXN (1, I) .GT. 0) .AND. (LXN (2, I) .GT. 0) .AND.
     &      (LXN (4, I) .EQ. 0)) THEN

C  SEE IF A 2-LINE NODE NEEDS DELETED

            IF (LXN (3, I) .LE. 0) THEN
               NODE = I
               KELEM = KXL (1, LXN (1, NODE))
               CHECK = .FALSE.
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            NODE, ERR)
               IF (ERR) GOTO 180
               CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &            NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3,
     &            DONE, CHECK, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 180
               IF (DONE) THEN

                  IF (VIDEO) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     CALL SNAPIT (3)
                  ENDIF

                  CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL,
     &               LXN, LLL, NNN, NNN, LNODES, BNSIZE, NLOOP, XMIN,
     &               XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)

                  IF ((GRAPH) .OR. (VIDEO)) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     IF (VIDEO) CALL SNAPIT (3)
                  ENDIF
                  REDO = .TRUE.
               ENDIF
            ENDIF
         ENDIF
  120 CONTINUE
      IF (REDO) GOTO 110

      DO 170 I = 1, NNN
         IF ((LXN (1, I) .GT. 0) .AND. (LXN (2, I) .GT. 0) .AND.
     &      (LXN (4, I) .EQ. 0)) THEN

C  GET THE ATTACHED LINES AND ELEMENTS
C  K(1) IS BETWEEN L1 & L2
C  K(2) IS BETWEEN L2 & L3
C  K(3) IS BETWEEN L3 & L1

            L1 = LXN (1, I)
            L2 = LXN (2, I)
            L3 = LXN (3, I)
            N1 = NXL (1, L1) + NXL (2, L1) - I
            N2 = NXL (1, L2) + NXL (2, L2) - I
            N3 = NXL (1, L3) + NXL (2, L3) - I

            IF ( (KXL (1, L1) .EQ. KXL (1, L2)) .OR.
     &         (KXL (1, L1) .EQ. KXL (2, L2)) ) THEN
               K(1) = KXL (1, L1)
            ELSEIF ( (KXL (2, L1) .EQ. KXL (1, L2)) .OR.
     &         (KXL (2, L1) .EQ. KXL (2, L2)) ) THEN
               K(1) = KXL (2, L1)
            ELSE
               ERR = .TRUE.
               CALL MESAGE ('** PROBLEMS IN TRIDEL FINDING K(1) **')
               GOTO 180
            ENDIF

            IF ( (KXL (1, L2) .EQ. KXL (1, L3)) .OR.
     &         (KXL (1, L2) .EQ. KXL (2, L3)) ) THEN
               K(2) = KXL (1, L2)
            ELSEIF ( (KXL (2, L2) .EQ. KXL (1, L3)) .OR.
     &         (KXL (2, L2) .EQ. KXL (2, L3)) ) THEN
               K(2) = KXL (2, L2)
            ELSE
               ERR = .TRUE.
               CALL MESAGE ('** PROBLEMS IN TRIDEL FINDING K(2) **')
               GOTO 180
            ENDIF

            IF ( (KXL (1, L3) .EQ. KXL (1, L1)) .OR.
     &         (KXL (1, L3) .EQ. KXL (2, L1)) ) THEN
               K(3) = KXL (1, L3)
            ELSEIF ( (KXL (2, L3) .EQ. KXL (1, L1)) .OR.
     &         (KXL (2, L3) .EQ. KXL (2, L1)) ) THEN
               K(3) = KXL (2, L3)
            ELSE
               ERR = .TRUE.
               CALL MESAGE ('** PROBLEMS IN TRIDEL FINDING K(3) **')
               GOTO 180
            ENDIF

C  NOW CHECK K(1)'S, K(2)'S, AND K(3)'S ANGLE AT THE LINE JOINT.
C  THERE ARE THREE POSSIBILITIES FOR CHANGE:
C     1) ANYTHING OVER 175 DEGREES GETS THE CORRESPONDING ELEMENT
C        DELETED
C     2) ANYTHING OVER 150 AND HOOKED TO ANOTHER 3-LINE NODE GETS
C        THE CORRESPONDING ELEMENT DELETED
C     3) AN ELONGATED ELEMENT OVER 150 DEGREES GETS A 3 ELEMENT
C        REPLACEMENT FOR THE TWO ELEMENTS THERE

            TOLER1 = 2.9670597
            TOLER2 = 2.6179939
            IF ((GRAPH) .AND. (.NOT. VIDEO)) THEN
               DIST = MAX (ABS(XN (N1) - XN (I)), ABS(XN (N2) - XN (I)),
     &            ABS(XN (N3) - XN (I)), ABS(YN (N1) - YN (I)),
     &            ABS(YN (N2) - YN (I)), ABS(YN (N3) - YN (I))) * 3.
               XMIN = XN (I) - DIST
               XMAX = XN (I) + DIST
               YMIN = YN (I) - DIST
               YMAX = YN (I) + DIST
            ENDIF

            ANG1 = ATAN2 (YN (N1) - YN (I), XN (N1) - XN (I))
            IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
            ANG2 = ATAN2 (YN (N2) - YN (I), XN (N2) - XN (I))
            IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
            ANG3 = ATAN2 (YN (N3) - YN (I), XN (N3) - XN (I))
            IF (ANG3 .LT. 0.) ANG3 = ANG3 + TWOPI

C  CHECK TO SEE IF THE NODES ARE CLOCKWISE OR COUNTERCLOCKWISE
C  (POSITIVE AREA IS CCW)

            AREA = ( (YN (N1) + YN (N3)) * .5 * (XN (N3) - XN (N1)) ) +
     &         ( (YN (N2) + YN (N1)) * .5 * (XN (N1) - XN (N2)) ) +
     &         ( (YN (N3) + YN (N2)) * .5 * (XN (N2) - XN (N3)) )

            IF (AREA .GT. 0.) THEN
               ANG12 = ANG2 - ANG1
               ANG23 = ANG3 - ANG2
               ANG31 = ANG1 - ANG3
            ELSE
               ANG12 = ANG1 - ANG2
               ANG23 = ANG2 - ANG3
               ANG31 = ANG3 - ANG1
            ENDIF
            IF (ANG12 .LT. 0.) ANG12 = ANG12 + TWOPI
            IF (ANG23 .LT. 0.) ANG23 = ANG23 + TWOPI
            IF (ANG31 .LT. 0.) ANG31 = ANG31 + TWOPI

            IF (GRAPH) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
C  5 IS PINK; 4 IS BLUE; 3 IS YELLOW; 0 IS BLACK ; 7 IS WHITE; 1 IS RED
               CALL LCOLOR ('YELOW')
               CALL D2NODE (MXND, XN, YN, I, N1)
               CALL D2NODE (MXND, XN, YN, I, N2)
               CALL D2NODE (MXND, XN, YN, I, N3)
               CALL LCOLOR ('WHITE')
               CALL SFLUSH
            ENDIF

C  NOW DO THE CHECKS FOR CHANGING THE ELEMENT

            IF (AREA .GT. 0) THEN
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(1), ANG12,
     &            TOLER1, TOLER2, N2, N1, N3, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ELSE
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(1), ANG12,
     &            TOLER1, TOLER2, N1, N2, N3, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ENDIF
            IF ((NOROOM) .OR. (ERR)) GOTO 180
            IF (DONE) GOTO 130
            IF (AREA .GT. 0) THEN
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(2), ANG23,
     &            TOLER1, TOLER2, N3, N2, N1, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ELSE
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(2), ANG23,
     &            TOLER1, TOLER2, N2, N3, N1, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ENDIF
            IF ((NOROOM) .OR. (ERR)) GOTO 180
            IF (DONE) GOTO 130
            IF (AREA .GT. 0) THEN
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(3), ANG31,
     &            TOLER1, TOLER2, N1, N3, N2, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ELSE
               CALL ADJTRI (MXND, MLN, LNODES, XN, YN, ZN, NUID, LXK,
     &            KXL, NXL, LXN, NNN, NAVAIL, IAVAIL, I, K(3), ANG31,
     &            TOLER1, TOLER2, N3, N1, N2, KREG, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, KKK, LLL, DEV1, DONE, CHECK, GRAPH,
     &            VIDEO, NOROOM, ERR, KKKADD)
            ENDIF
            IF ((NOROOM) .OR. (ERR)) GOTO 180

  130       CONTINUE
            IF (DONE) THEN

               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N1, ERR)
               IF (ERR) GOTO 180
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N2, ERR)
               IF (ERR) GOTO 180
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N3, ERR)
               IF (ERR) GOTO 180
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            I, ERR)
               IF (ERR) GOTO 180
               IF (VIDEO) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  CALL SNAPIT (3)
               ENDIF

               CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &            LLL, NNN, NNN, LNODES, BNSIZE, NLOOP, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  IF (VIDEO) CALL SNAPIT (3)
               ENDIF
               DONE = .FALSE.
               REDO = .TRUE.
               GOTO 160
            ENDIF

C  NOW CHECK THE THREE ELEMENTS TO SEE IF AN ELEMENT EXISTS WHICH:
C    (1) CONTAINS ONLY 2 OPPOSING "LARGE ANGLE" THREE-LINE NODES
C        AND AT LEAST 1 "SMALL ANGLE" FOUR- OR FIVE-LINE NODE
C    (2) ONE THAT IS REALLY SQUASHED WITH AT LEAST ONE
C        "SMALL ANGLE" FIVE-LINE NODE.
C    (3) CONTAINS A "LARGE ANGLE" THREE-LINE NODE AND TWO
C        RELATIVELY SMALL FIVE-LINE NODE AND A NORMAL FOUR
C        LINE NODE
C    (4) CONTAINS TWO "VERY SMALL" ANGLES AND TWO "VERY LARGE"
C        ANGLES
C    (5) CONTAINS TWO RELATIVELY SMALL ANGLES AND TWO RELATIVELY
C        LARGE ANGLES AND IS CONSIDERABLY SMALLER THAN IS
C        DICTATED BY THE DESIRED SIZE
C  THIS ELEMENT SHOULD BE DELETED.

            TOLER3 = 1.7453293
            TOLER4 = 1.5707963
            TOLER5 = 2.0943951
            TOLER6 = 0.9599311
            DO 150 J = 1, 3
               CALL GNXKA (MXND, XN, YN, K(J), NODES, AREA, LXK, NXL,
     &            CCW)
               IF ( (I .NE. NODES(1)) .AND. (I .NE. NODES(2)) .AND.
     &            (I .NE. NODES(3)) .AND. (I .NE. NODES(4)) ) THEN
                  CALL MESAGE ('** PROBLEMS IN TRIDEL - I IS NOT IN '//
     &               'ELEMENT K **')
                  ERR = .TRUE.
                  GOTO 180
               ENDIF

C  ARRANGE NODES SO THE COLLAPSING DIAGONAL IS FROM 1ST TO 3RD NODES
C  AND INSURE THAT THE NODE TO BE DELETED IS NOT A BOUNDARY NODE

               CALL NXKORD (NODES, I)
               N1 = NODES(1)
               N2 = NODES(2)
               N3 = NODES(3)
               N4 = NODES(4)

               X21 = XN (N2) - XN (N1)
               X32 = XN (N3) - XN (N2)
               X43 = XN (N4) - XN (N3)
               X14 = XN (N1) - XN (N4)
               Y21 = YN (N2) - YN (N1)
               Y32 = YN (N3) - YN (N2)
               Y43 = YN (N4) - YN (N3)
               Y14 = YN (N1) - YN (N4)
               D21 = SQRT (X21 **2 + Y21 **2)
               D32 = SQRT (X32 **2 + Y32 **2)
               D43 = SQRT (X43 **2 + Y43 **2)
               D14 = SQRT (X14 **2 + Y14 **2)
               DMAX = MAX (D21, D32, D43, D14) * 1.5
               IF (LXN (3, N1) .EQ. 0) THEN
                  THETA1 = PI
               ELSE
                  THETA1 = ACOS (- ( (X21 * X14) + (Y21 * Y14) ) /
     &               (D21 * D14))
               ENDIF
               IF (LXN (3, N2) .EQ. 0) THEN
                  THETA2 = PI
               ELSE
                  THETA2 = ACOS (- ( (X32 * X21) + (Y32 * Y21) ) /
     &               (D32 * D21))
               ENDIF
               IF (LXN (3, N3) .EQ. 0) THEN
                  THETA3 = PI
               ELSE
                  THETA3 = ACOS (- ( (X43 * X32) + (Y43 * Y32) ) /
     &               (D43 * D32))
               ENDIF
               IF (LXN (3, N3) .EQ. 0) THEN
                  THETA4 = PI
               ELSE
                  THETA4 = ACOS (- ( (X14 * X43) + (Y14 * Y43) ) /
     &               (D14 * D43))
               ENDIF

C  TEST CASE ONE

               IF ( (LXN(2, N1) .GT. 0) .AND.
     &            (LXN (2, N3) .GT. 0) .AND.
     &            (LXN (4, N3) .EQ. 0) .AND.
     &            (LXN (4, N2) .NE. 0) .AND.
     &            (LXN (4, N4) .NE. 0) .AND.
C     &            ( (LXN (4, N2) .LT. 0) .OR.
C     &            (LXN (4, N4) .LT. 0) ) .AND.
     &            ((THETA1 .GT. TOLER3) .OR. (THETA3 .GT. TOLER3)) .AND.
     &            ((THETA2 .LT. TOLER4) .OR. (THETA4 .LT. TOLER4)) .AND.
     &            (K (J) .NE. KKKADD))
     &            THEN
                  PASSED = .TRUE.

C  TEST CASE 2

               ELSEIF ( (LXN(2, N1) .GT. 0) .AND.
     &            (LXN (2, N3) .GT. 0) .AND.
     &            (LXN (4, N3) .GE. 0) .AND.
     &            (LXN (4, N2) .NE. 0) .AND.
     &            (LXN (4, N4) .NE. 0) .AND.
     &            ( (LXN (4, N2) .LT. 0) .OR.
     &            (LXN (4, N4) .LT. 0) ) .AND.
     &            ((THETA1 .GT. TOLER5) .OR. (THETA3 .GT. TOLER5)) .AND.
     &            ((THETA2 .LT. TOLER6) .OR. (THETA4 .LT. TOLER6)) .AND.
     &            (K (J) .NE. KKKADD) )
     &            THEN
                  PASSED = .TRUE.

C  TEST CASE 3

               ELSEIF ( (LXN(2, N1) .GT. 0) .AND.
     &            (LXN (2, N3) .GT. 0) .AND.
     &            (LXN (4, N3) .GE. 0) .AND.
     &            (LXN (4, N2) .LT. 0) .AND.
     &            (LXN (4, N4) .LT. 0) .AND.
     &            ((THETA1 .GT. TOLER3) .OR. (THETA3 .GT. TOLER3)) .AND.
     &            ((THETA2 .LT. TOLER4) .OR. (THETA4 .LT. TOLER4)) .AND.
     &            (K (J) .NE. KKKADD) )
     &            THEN
                  PASSED = .TRUE.

C  TEST CASE 4

               ELSEIF ( (LXN(2, N1) .GT. 0) .AND.
     &            (LXN (2, N3) .GT. 0) .AND.
     &            (THETA1 .GT. TOLER5) .AND.
     &            (THETA3 .GT. TOLER5) .AND.
     &            (THETA2 .LT. TOLER6) .AND.
     &            (THETA4 .LT. TOLER6) .AND.
     &            (K (J) .NE. KKKADD) )
     &            THEN
                  PASSED = .TRUE.

C  TEST CASE 5

               ELSEIF ( (LXN(2, N1) .GT. 0) .AND.
     &            (LXN (2, N3) .GT. 0) .AND.
     &            (THETA1 .GT. TOLER3) .AND.
     &            (THETA3 .GT. TOLER3) .AND.
     &            (THETA2 .LT. TOLER4) .AND.
     &            (THETA4 .LT. TOLER4) .AND.
     &            (DMAX .LT. BNSIZE (1, N1)) .AND.
     &            (K (J) .NE. KKKADD) )
     &            THEN
                  PASSED = .TRUE.

               ELSE
                  PASSED = .FALSE.
               ENDIF

               IF (PASSED) THEN
                  IF (GRAPH) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     CALL LCOLOR ('PINK ')
                     CALL D2NODE (MXND, XN, YN, N1, N2)
                     CALL D2NODE (MXND, XN, YN, N2, N3)
                     CALL D2NODE (MXND, XN, YN, N3, N4)
                     CALL D2NODE (MXND, XN, YN, N4, N1)
                     CALL LCOLOR ('WHITE')
                     CALL SFLUSH
                  ENDIF
                  NODE = N1
                  KELEM = K(J)
  140             CONTINUE
                  CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &               N1, ERR)
                  IF (ERR) GOTO 180
                  CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &               N2, ERR)
                  IF (ERR) GOTO 180
                  CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &               N3, ERR)
                  IF (ERR) GOTO 180
                  CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &               N4, ERR)
                  IF (ERR) GOTO 180
                  CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &               NNN, NAVAIL, IAVAIL, NODE, KELEM, NODE1, NODE3,
     &               DONE, CHECK, NOROOM, ERR)
                  IF ((NOROOM) .OR. (ERR)) GOTO 180
                  IF (DONE) THEN

                     IF (VIDEO) THEN
                        CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                     YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                        CALL SNAPIT (3)
                     ENDIF

                     CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL,
     &                  LXN, LLL, NNN, NNN, LNODES, BNSIZE, NLOOP, XMIN,
     &                  XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)

                     IF ((GRAPH) .OR. (VIDEO)) THEN
                        CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                     YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                        IF (VIDEO) CALL SNAPIT (3)
                     ENDIF

C  CHECK TO SEE IF WE HAVE MADE A 2-LINE NODE

                     IF (LXN (3, NODE1) .LE. 0) THEN
                        NODE = NODE1
                        KELEM = KXL (1, LXN (1, NODE))
                        CHECK = .FALSE.
                        GOTO 140
                     ELSEIF (LXN (3, NODE3) .LE. 0) THEN
                        NODE = NODE3
                        KELEM = KXL (1, LXN (1, NODE))
                        CHECK = .FALSE.
                        GOTO 140
                     ENDIF

                     CHECK = .TRUE.
                     DONE = .FALSE.
                     REDO = .TRUE.
                     GOTO 160
                  ENDIF
               ENDIF
  150       CONTINUE
  160       CONTINUE
         ENDIF
  170 CONTINUE

      CALL TRIFIX (MXND, MLN, XN, YN, ZN, NUID, LXK, KXL, NXL, LXN,
     &   NNN, LLL, KKK, NAVAIL, IAVAIL, ANGLE, LNODES, BNSIZE,
     &   NLOOP, DEV1, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, GRAPH,
     &   VIDEO, NOROOM, ERR)
      IF ((NOROOM) .OR. (ERR)) GOTO 180
      IF (REDO) GOTO 100
  180 CONTINUE

      RETURN

      END
