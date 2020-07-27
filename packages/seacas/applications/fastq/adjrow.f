C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADJROW (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, ANGLE, BNSIZE, LNODES, NLOOP, IAVAIL, NAVAIL, XMIN, XMAX,
     &   YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN, LLLOLD, NNNOLD,
     &   NODE, NADJ1, NADJ2, NNN2, GRAPH, VIDEO, KREG, DEFSIZ, ADJTED,
     &   NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE ADJROW = ADJUSTS A ROW OF ELEMENTS BETWEEN TWO CORNERS

C***********************************************************************

      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES

      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)

      LOGICAL ERR, GRAPH, ADJTED, VIDEO, NOROOM

      CHARACTER*3 DEV1

      DATA TMIN1 /.80/, TMIN2 /.3/, WMIN1 /1.25/, WMIN2 /1.35/

      PI = ATAN2(0.0, -1.0)
      CALL GETIME (TIME1)
      ERR = .FALSE.
      EPS = .0523599

C  START BY SETTING UP THE LIMITS OF THE SEARCH

      IF (NADJ1 .EQ. NADJ2) THEN
         N2 = LNODES (3, NADJ1)
         KOUNT = 0
  100    CONTINUE
         KOUNT = KOUNT + 1
         IF ((ANGLE (N2) .GE. PI - EPS) .AND.
     &      (ANGLE (N2) .LE. PI + EPS)) THEN
            NADJ1 = N2
            NADJ2 = N2
            TEPS = .95 * (PI -
     &         ( (DBLE(NLOOP - 2) * PI) / DBLE(NLOOP)))
            IF (TEPS .LE. EPS) EPS = TEPS
            GOTO 110
         ELSEIF (N2 .EQ. NADJ2) THEN
            TEPS = .95 * (PI -
     &         ( (DBLE(NLOOP - 2) * PI) / DBLE(NLOOP)))
            IF (TEPS .LE. EPS) EPS = TEPS
            GOTO 110
         ELSEIF (KOUNT .GT. NLOOP) THEN
            CALL MESAGE ('** PROBLEMS IN ADJROW WITH LOOP NOT '//
     &         'CLOSING **')
            ERR = .TRUE.
            GOTO 160
         ELSE
            N2 = LNODES (3, N2)
            GOTO 100
         ENDIF
      ENDIF

  110 CONTINUE
      N1 = LNODES (3, NADJ1)
      ADJTED = .FALSE.

  120 CONTINUE
      IF (N1 .EQ. NADJ2) GOTO 150

C  CHECK A STRING OF CONCAVE (< PI) INTERIOR ANGLES FOR NEEDING A
C  TUCK INSERTED SOMEWHERE

      IF ((ANGLE (N1) .LT. PI - EPS) .AND. (LNODES (8, N1) .GT. 1) .AND.
     &   (LXN (4, N1) .EQ. 0) .AND. (LXN (3, N1) .GT. 0)) THEN

C  ADDED UP THE TURNING ANGLE AND THE AVERAGE SIZE REDUCTION

         TANG = 0.
         KANG = 0
         RATIO = 0.
         N11 = N1
  130    CONTINUE
         TANG = TANG + (PI - ANGLE (N11) )
         KANG = KANG + 1
         N0 = LNODES (2, N11)
         N2 = LNODES (3, N11)
         DIST = .5 * (SQRT ( ((XN (N0) - XN (N11)) ** 2) +
     &      ((YN (N0) - YN (N11)) ** 2) ) +
     &      SQRT ( ((XN (N2) - XN (N11)) ** 2) +
     &      ((YN (N2) - YN (N11)) ** 2) ) )
         IF (DEFSIZ .GT. 0.) THEN
            IF (DIST .LT. DEFSIZ) THEN
               RATIO = RATIO + ( DIST / BNSIZE (1, N11) )
            ELSE
               RATIO = RATIO + ( DIST / DEFSIZ)
            ENDIF
         ELSE
            RATIO = RATIO + ( DIST / BNSIZE (1, N11) )
         ENDIF
         N11 = LNODES (3, N11)
         IF ((N11 .NE. NADJ2) .AND. (ANGLE (N11) .LT. PI - EPS) .AND.
     &      (LXN (4, N11) .EQ. 0) .AND. (LXN (3, N11) .GT. 0) .AND.
     &      (LNODES (8, N11) .GT. 1)) GOTO 130
         KANG = KANG

C  NOW SEE IF THIS PORTION OF THE ROW NEEDS ADJUSTED WITH A TUCK(S)

         IF (KANG .GE. 1) THEN
            RATIO = RATIO / DBLE(KANG)

C**               THIS CRITERIA SHOULD HELP ALLEVIATE THE LONG SKINNY
C**               ELEMENT FORMATIONS WHEN TRANSITIONING.

            IF ( ((RATIO .LT. TMIN1) .AND. (TANG .GT. 1.2217)) .OR.
     &         ((RATIO .LT. TMIN2) .AND. (TANG .GT. .9)) ) THEN
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  IF (VIDEO) CALL SNAPIT (1)
               ENDIF
               N11OLD = N11
               CALL ADDTUK (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &            LXN, LNODES, ANGLE, NLOOP, IAVAIL, NAVAIL, LLL, KKK,
     &            NNN, TANG, KANG, N1, N11, NODE, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, GRAPH, VIDEO, DEV1, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 160

C  MAKE SURE THAT THE TUCK DOES NOT ELIMINATE THE END NODES FOR THE LOOP

               IF (N11 .NE. N11OLD) THEN
                  IF (NADJ2 .EQ. N11OLD) NADJ2 = N11
                  IF (NODE .EQ. N11OLD) NODE = N11
               ENDIF

               NNNOLD = NNN
               LLLOLD = LLL
               ADJTED = .TRUE.

            ENDIF
         ENDIF
         N1 = N11
         GOTO 120

C  CHECK A STRING OF CONVEX (> PI) INTERIOR ANGLES FOR NEEDING A
C  WEDGE INSERTED SOMEWHERE

      ELSEIF ((ANGLE (N1) .GE. PI + EPS) .AND. (LXN (3, N1) .GT. 0)
     &   .AND. (LXN (4, N1) .EQ. 0)) THEN

C  ADD UP THE TURNING ANGLE AND THE AVERAGE SIZE REDUCTION

         TANG = 0.
         KANG = 0
         RATIO = 0.
         IDEPTH = 0
         N11 = N1
  140    CONTINUE
         TANG = TANG + (ANGLE (N11) - PI)
         KANG = KANG + 1
         N0 = LNODES (2, N11)
         N2 = LNODES (3, N11)
         DIST = .5 * (SQRT ( ((XN (N0) - XN (N11)) ** 2) +
     &      ((YN (N0) - YN (N11)) ** 2) ) +
     &      SQRT ( ((XN (N2) - XN (N11)) ** 2) +
     &      ((YN (N2) - YN (N11)) ** 2) ) )
         IF (DEFSIZ .GT. 0.) THEN
            IF (DIST .GT. DEFSIZ) THEN
               RATIO = RATIO + ( DIST / BNSIZE (1, N11) )
            ELSE
               RATIO = RATIO + ( DIST / DEFSIZ)
            ENDIF
         ELSE
            RATIO = RATIO + ( DIST / BNSIZE (1, N11) )
         ENDIF
         N11 = LNODES (3, N11)
         IDEPTH = MAX (IDEPTH, LNODES (8, N11))
         IF ((N11 .NE. NADJ2) .AND. (ANGLE (N11) .GE. PI + EPS) .AND.
     &      (LXN (4, N11) .EQ. 0) .AND. (LXN (3, N11) .GT. 0)) GOTO 140

C  NOW SEE IF THIS PORTION OF THE ROW NEEDS ADJUSTED WITH A WEDGE(S)

         IF (KANG .GE. 1) THEN
            RATIO = RATIO / DBLE(KANG)
            IF ( ( ((RATIO .GT. WMIN1) .AND. (IDEPTH .GT. 1)) .OR.
     &         ((RATIO .GT. WMIN2) .AND. (IDEPTH .EQ. 1)) )
     &         .AND. (TANG .GT. 1.2217)) THEN
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  IF (VIDEO) CALL SNAPIT (1)
               ENDIF
               CALL ADDWDG (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &            LXN, LNODES, ANGLE, BNSIZE, NLOOP, IAVAIL, NAVAIL,
     &            LLL, KKK, NNN, LLLOLD, NNNOLD, TANG, KANG, N1, N11,
     &            XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, GRAPH, VIDEO,
     &            DEV1, KREG, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 160
               NNNOLD = NNN
               LLLOLD = LLL
               ADJTED = .TRUE.

            ENDIF
         ENDIF
         N1 = N11
         GOTO 120
      ELSE
         N1 = LNODES (3, N1)
         GOTO 120
      ENDIF

C  NOW SMOOTH, CALCULATE THE NEW ANGLES, AND PLOT IF NEEDED

  150 CONTINUE
      IF (ADJTED) THEN
         CALL GETIME (TIME2)
         TIMEAJ = TIMEAJ + TIME2 - TIME1
         CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &      LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX, YMIN,
     &      YMAX, ZMIN, ZMAX, DEV1, KREG)
         CALL GETIME (TIME1)
         CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NLOOP,
     &      ANGLE, LNODES, N1, LLL, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     &      DEV1, KREG, ERR)
         IF (ERR) GOTO 160
         IF ((GRAPH) .OR. (VIDEO)) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
      ENDIF

  160 CONTINUE

      CALL GETIME (TIME2)
      TIMEAJ = TIMEAJ + TIME2 - TIME1
      RETURN

      END
