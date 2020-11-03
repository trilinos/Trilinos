C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TRIFIX (MXND, MLN, XN, YN, ZN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, ANGLE, LNODES, BNSIZE,
     &   NLOOP, DEV1, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAXZ, GRAPH,
     &   VIDEO, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE TRIFIX = CHECKS ALL ELEMENTS FOR ANY TRIANGULAR SHAPED
C                      LONG ELEMENT AND DELETES THEM WHEN
C                      FOUND AND POSSIBLE

C***********************************************************************

      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
      DIMENSION NODES(4)
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), ZN(MXND), NUID(MXND)

      CHARACTER*3 DEV1
      LOGICAL ERR, DONE, GRAPH, REDO, CCW
      LOGICAL VIDEO, NOROOM

      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI

      ERR = .FALSE.
      DONE = .FALSE.
      CCW = .TRUE.
      KMAX = 30
      KOUNT = 0

C  TOLERANCE IS SET AT 150 DEGREES

      TOLER = 2.6179939

  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. KMAX) GOTO 140
      REDO = .FALSE.

      DO 130 KELEM = 1, KKK
         IF (LXK (1, KELEM) .GT. 0) THEN
            CALL GNXKA (MXND, XN, YN, KELEM, NODES, AREA, LXK, NXL, CCW)
            DO 110 I = 1, 4
               I1 = NODES (I)
               IF (I .EQ. 1) THEN
                  I0 = NODES (4)
                  I2 = NODES (2)
               ELSEIF (I .EQ. 4) THEN
                  I0 = NODES (3)
                  I2 = NODES (1)
               ELSE
                  I0 = NODES (I - 1)
                  I2 = NODES (I + 1)
               ENDIF

               ANG1 = ATAN2 (YN (I0) - YN (I1), XN (I0) - XN (I1))
               IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
               ANG2 = ATAN2 (YN (I2) - YN (I1), XN (I2) - XN (I1))
               IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
               ANG = ANG1 - ANG2
               IF (ANG .LT. 0.) ANG = ANG + TWOPI

               CALL LONGEL (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL,
     &            NXL, LXN, NNN, NAVAIL, IAVAIL, I1, KELEM, ANG, TOLER,
     &            I0, I2, KREG, XMIN, XMAX, YMIN, YMAX, KKK, LLL,
     &            DONE, GRAPH, VIDEO, NOROOM, ERR, KKKADD)
               IF ((NOROOM) .OR. (ERR)) GOTO 140

               IF (DONE) THEN

                  IF ((GRAPH) .AND. (.NOT. VIDEO)) THEN
                     DIST = MAX (ABS(XN (I0) - XN (I1)),
     &                  ABS(XN (I2) - XN (I1)), ABS(YN (I0) - YN (I1)),
     &                  ABS(YN (I2) - YN (I1))) * 3.
                     XMIN = XN (I1) - DIST
                     XMAX = XN (I1) + DIST
                     YMIN = YN (I1) - DIST
                     YMAX = YN (I1) + DIST
                  ENDIF
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
                  DONE = .FALSE.
                  REDO = .TRUE.
                  GOTO 120
               ENDIF

  110       CONTINUE
  120       CONTINUE
         ENDIF
  130 CONTINUE

      IF (REDO) GOTO 100
  140 CONTINUE

      RETURN

      END
