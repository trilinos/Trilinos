C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PINCH (MXND, MXCORN, MLN, NUID, XN, YN, ZN, LXK, KXL,
     &   NXL, LXN, ANGLE, LNODES, BNSIZE, NODE, NLOOP, KKKOLD, LLLOLD,
     &   NNNOLD, IAVAIL, NAVAIL, DONE, XMIN, XMAX, YMIN, YMAX, ZMIN,
     &   ZMAX, DEV1, LLL, KKK, NNN, LCORN, NCORN, NADJ1, NADJ2, GRAPH,
     &   VIDEO, KREG, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE PINCH = COLLAPSES A CORNER WITH A SMALL ANGLE CLOSED

C***********************************************************************

      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES

      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION LCORN (MXCORN)
      DIMENSION L1LIST(20)

      LOGICAL DONE, NOROOM, ERR, FOUND, DDONE, PWEDGE, WEDGOK
      LOGICAL GRAPH, VIDEO, PGRAPH, ONLYC, BOK, PPOSBL, CHECK

      CHARACTER*3 DEV1

      PI = ATAN2(0.0, -1.0)

      CALL GETIME (TIME1)
      PGRAPH = .FALSE.
      PWEDGE = .TRUE.
      FOUND = .FALSE.
      CHECK = .FALSE.
      ONLYC = .TRUE.
      DONE = .FALSE.
      ERR = .FALSE.

C  SEE IF ONLY 2 NODES ARE LEFT ON THE LOOP

      IF (NLOOP .EQ. 2) THEN
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NODE, ERR)
         IF (ERR) GOTO 210
         CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &      LNODES (3, NODE), ERR)
         IF (ERR) GOTO 210
         CALL CLOSE2 (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &      LNODES, IAVAIL, NAVAIL, NNN, LLL, NODE, XMIN, XMAX, YMIN,
     &      YMAX, ZMIN, ZMAX, PGRAPH, VIDEO, DEV1, KREG, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 210
         NLOOP = 0
         FOUND = .TRUE.
         DONE = .TRUE.
         NNN2 = NNN
         CALL GETIME (TIME2)
         TIMEP = TIMEP + TIME2 - TIME1
         CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, LLL,
     &      NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, DEV1, KREG)
         CALL GETIME (TIME1)
         IF ((PGRAPH) .OR. (VIDEO)) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &         ZMIN, ZMAX, LLL, DEV1, KREG)
            IF (VIDEO) CALL SNAPIT (1)
         ENDIF
         GOTO 200
      ENDIF
      IF (PGRAPH) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, LLL, DEV1, KREG)
      ENDIF

C  GET THE CORNERS THAT CAN BE ADJUSTED

      N1OLD = 0
  100 CONTINUE
      IF (GRAPH) CALL LCOLOR ('YELOW')
      CALL GETCRN (MXND, MXCORN, MLN, LNODES, NCORN, LCORN,
     &   ANGLE, XN, YN, LXN, NLOOP, NODE, ONLYC, PPOSBL, GRAPH,
     &   ERR)
      IF (NCORN .EQ. 2) THEN
         IDIF = MIN0 (LNODES (7, LCORN(1)), LNODES (7, LCORN(2)) )
      ELSEIF (NCORN .EQ. 3) THEN
         ILOW = MIN0 (LNODES (7, LCORN(1)), LNODES (7, LCORN(2)),
     &      LNODES (7, LCORN(3)) )
         IHIGH = MAX0 (LNODES (7, LCORN(1)), LNODES (7, LCORN(2)),
     &      LNODES (7, LCORN(3)) )
      ENDIF
      IF (ERR) GOTO 210
      IF (GRAPH) CALL LCOLOR ('WHITE')

C  FOR NORMAL REGIONS,
C  TOLER1 IS SET AT 20 DEGREES. (A 3 DEGREE IRREGULAR NODE IS FORMED)
C  TOLER2 IS SET TO 50 DEGREES. (A 4+ DEGREE IRREGULAR NODE IS HELPED)

C  THEY ARE SET FOR AN UNEQUAL SEMICIRCLE TO 30 AND 60 RESPECTIVELY
C  THEY ARE SET FOR AN EQUAL SEMICIRCLE TO 35 AND 70 RESPECTIVELY
C  THEY ARE SET FOR A 3-2-1 TRIANGLE TO 35 AND 70 RESPECTIVELY

      IF (NCORN .EQ. 2) THEN
         IDIF = MIN0 (LNODES (7, LCORN(1)), LNODES (7, LCORN(2)) )
         IF (LNODES (7, LCORN(1)) .EQ. LNODES (7, LCORN (2)) ) THEN
            TOLER1 = .5235988
            TOLER2 = 1.2217305
         ELSE
            TOLER1 = .6108652
            TOLER2 = 1.0471976
         ENDIF
      ELSEIF ((NCORN .EQ. 3) .AND. (ILOW .EQ. 1) .AND.
     &   (IHIGH .EQ. 3)) THEN
         TOLER1 = .6108652
         TOLER2 = 1.0471976
      ELSE
         TOLER1 = .3490659
         TOLER2 = .8726646
      ENDIF

C  NOW MAKE SURE THAT A WEDGE CAN BE ALLOWED

      IF (NLOOP .LE. 4) THEN
         KNEG = 0
         DO 110 I = 1, NCORN
            IF (ANGLE (LCORN (I)) .LT. 0.) KNEG = KNEG + 1
  110    CONTINUE
         IF (KNEG .GE. 2) THEN
            WEDGOK = .FALSE.
         ELSE
            WEDGOK = .TRUE.
         ENDIF
      ELSE
         WEDGOK = .TRUE.
      ENDIF

C  NOW SORT THE CORNERS SO THE SMALLEST REMAINING ONE GOES FIRST

  120 CONTINUE
      J = 0
      DO 130 I = 1, NCORN
         IF (LCORN (I) .GT. 0) THEN
            IF (J .EQ. 0) THEN
               J = I
            ELSEIF (ANGLE (LCORN (I)) .LT. ANGLE (LCORN (J))) THEN
               J = I
            ENDIF
         ENDIF
  130 CONTINUE

      IF (J .GT. 0) THEN
         N1 = LCORN (J)
         LCORN (J) = - LCORN (J)
         N0 = LNODES (2, N1)
         N2 = LNODES (3, N1)

C  CHECK TO MAKE SURE THAT A 1-1-1-1 RECTANGLE ISN'T BEING CLOSED

         IF ((NLOOP .LE. 4) .AND. (NCORN .GE. 4)) GOTO 200

C  CHECK TO MAKE SURE THAT A 4 - 1 - 1 TRIANGLE ISN'T BEING CLOSED

C         ELSEIF ((NCORN .EQ. 3) .AND. (NLOOP .EQ. 6) .AND.
C     &      (ILOW .EQ. 1) .AND. (IHIGH .EQ. 4) ) THEN
C            GOTO 200
C         ENDIF

C  CHECK TO MAKE SURE THAT THE ANGLE IS ELIGIBLE FOR PINCHING AND
C  THAT A CLOSURE DOESN'T FORM A DEGENERATE ELEMENT ALONG THE BOUNDARY

         CALL BPINCH (MXND, MLN, LNODES, XN, YN, LXN, NXL, ANGLE,
     &      N0, N1, N2, NLOOP, TOLER1, TOLER2, BOK, ERR)
         IF (ERR) GOTO 210
         IF (BOK) THEN
            IF (NCORN .EQ. 2) IDIF = IDIF - 1

C  CHECK TO SEE IF A WEDGE NEEDS TO BE ADDED BEFORE THE THING IS PINCHED

            DIST01 = SQRT ( ((XN (N1) - XN (N0)) **2) +
     &         ((YN (N1) - YN (N0)) **2) )
            DIST21 = SQRT ( ((XN (N1) - XN (N2)) **2) +
     &         ((YN (N1) - YN (N2)) **2) )
            FACT = 2.5
            IF ((WEDGOK) .AND. (DIST01 .GT. FACT * DIST21) .AND.
     &         (KXL (1, LNODES (5, N0)) .GT. 0) ) THEN
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, N1), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, N1)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, LNODES (2, N1))), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, N1), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, N1)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, LNODES (3, N1))), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N1, ERR)
               IF (ERR) GOTO 210
               CALL WEDGE (MXND, MLN, NUID, LXK, KXL, NXL, LXN, XN, YN,
     &            LNODES, BNSIZE, IAVAIL, NAVAIL, LLL, KKK, NNN, LLLOLD,
     &            NNNOLD, N1, N6, NLOOP, PWEDGE, GRAPH, VIDEO, NOROOM,
     &            ERR)

C  WATCH FOR THE REPEATING CASE

               IF (N1 .EQ. N1OLD) THEN
                  BNSIZE (2, N1) = BNSIZE (2, N1) * 3.
                  BNSIZE (2, LNODES (3, N6)) =
     &               BNSIZE (2, LNODES (3, N6)) * 3.
                  BNSIZE (2, N6) = BNSIZE (2, N6) * 3.
               ENDIF
               N1OLD = N1

               IF ((NOROOM) .OR. (ERR)) GOTO 210
               IF (VIDEO) CALL SNAPIT (2)
               IF (NODE .EQ. N1) NODE = LNODES (2, N2)
               IF (NADJ1 .EQ. N1) NADJ1 = LNODES (2, N2)
               IF (NADJ2 .EQ. N1) NADJ2 = LNODES (2, N2)
               ANGLE (LNODES (2, N2)) = ANGLE (N1)
               N1 = LNODES (2, N2)
               N0 = LNODES (2, N1)
               ANGLE (N1) = PI
               ANGLE (N0) = PI

            ELSEIF ((WEDGOK) .AND. (DIST21 .GT. FACT * DIST01) .AND.
     &         (LXN (3, N2) .NE. 0) .AND.
     &         (KXL (1, LNODES (5, N1)) .GT. 0) ) THEN
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, N2), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, N2)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, LNODES (2, N2))), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, N2), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, N2)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, LNODES (3, N2))), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N2, ERR)
               IF (ERR) GOTO 210
               CALL WEDGE (MXND, MLN, NUID, LXK, KXL, NXL, LXN, XN, YN,
     &            LNODES, BNSIZE, IAVAIL, NAVAIL, LLL, KKK, NNN, LLLOLD,
     &            NNNOLD, N2, N6, NLOOP, PWEDGE, GRAPH, VIDEO, NOROOM,
     &            ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 210
               IF (VIDEO) CALL SNAPIT (2)
               IF (NODE .EQ. N1) NODE = LNODES (3, N6)
               IF (NADJ2 .EQ. N1) NADJ2 = LNODES (3, N6)
               AHOLD = ANGLE (N2)
               N2 = LNODES (3, N1)
               ANGLE (N2) = PI
               ANGLE (LNODES (3, N2)) = PI
               ANGLE (LNODES (3, LNODES (3, N2))) = AHOLD
            ENDIF

C  PROCEED WITH THE PINCH

            LINE1 = LNODES (5, N0)
            LINE2 = LNODES (5, N1)

C  CHECK TO MAKE SURE THAT AT LEAST ONE OF THE LINES
C  IS NOT A BOUNDARY LINE AND GET THE NODE TO BE DELETED

            IF ((LXN (2, N0) .GT. 0) .OR.
     &         (LXN (2, N2) .GT. 0)) THEN

               FOUND = .TRUE.

               IF (LXN (2, N0) .GT. 0) THEN
                  NGONE = N0
                  NTHERE = N2
                  LNEW = LINE2
                  LOLD = LINE1
                  LNODES (2, NTHERE) = LNODES (2, N0)
                  LNODES (3, LNODES (2, N0)) = NTHERE
               ELSE
                  NGONE = N2
                  NTHERE = N0
                  LNEW = LINE1
                  LOLD = LINE2
                  LNODES (3, NTHERE) = LNODES (3, N2)
                  LNODES (2, LNODES (3, N2)) = NTHERE
                  LNODES (5, NTHERE) = LNODES (5, N2)
               ENDIF
               LNODES (4, N1) = - 2

C  SEE IF THE NODES BEING USED ARE IN THE CORNER LIST
C  IF THEY ARE THEN THOSE CORNERS ARE NEGATED

               DO 140 NC = 1, NCORN
                  IF ( (LCORN (NC) .EQ. NTHERE) .OR.
     &               (LCORN (NC) .EQ. NGONE) )
     &               LCORN (NC) = - IABS (LCORN (NC))
  140          CONTINUE

C  DELETE THE OLD LINE AND REDO LINK ARRAYS

               KOLD = KXL (1, LOLD)
               KNEW = KXL (1, LNEW)
               KXL (1, LNEW) = KNEW
               KXL (2, LNEW) = KOLD
               KXL (1, LOLD) = 0
               KXL (2, LOLD) = 0
               IF ((VIDEO) .OR. (PGRAPH)) THEN
                  CALL LCOLOR ('BLACK')
                  CALL D2NODE (MXND, XN, YN, NXL (1, LOLD),
     &               NXL (2, LOLD))
                  IF (GRAPH) THEN
                     CALL LCOLOR ('WHITE')
                  ELSE
                     CALL LCOLOR ('YELOW')
                  ENDIF
                  CALL SFLUSH
               ENDIF
               NXL (1, LOLD) = 0
               NXL (2, LOLD) = 0

C  FIX THE LINES PER ELEMENT ARRAY FOR THE ONE ELEMENT CHANGING

               IF (KOLD .GT. 0) THEN
                  DO 150 II = 1, 4
                     IF (LXK (II, KOLD) .EQ. LOLD) THEN
                        LXK (II, KOLD) = LNEW
                        GOTO 160
                     ENDIF
  150             CONTINUE
                  CALL MESAGE ('** PROBLEMS IN PINCH FIXING'//
     &               ' ELEMENT **')
                  ERR = .TRUE.
                  GOTO 210
  160             CONTINUE
               ENDIF

C  RECONNECT ALL LINES CONNECTING TO NGONE TO NTHERE

               CALL GETLXN (MXND, LXN, NGONE, L1LIST, NL, ERR)
               IF (ERR) THEN
                  CALL MESAGE ('** PROBLEMS IN PINCH GETTING NGONE'//
     &               'LINES **')
                  GOTO 210
               ENDIF
               DO 170 II = 1, NL
                  LL = L1LIST (II)
                  IF (NXL (1, LL) .EQ. NGONE) THEN
                     IF ((VIDEO) .OR. (PGRAPH)) THEN
                        CALL LCOLOR ('BLACK')
                        CALL D2NODE (MXND, XN, YN, NXL (1, LL),
     &                     NXL (2, LL))
                        IF (GRAPH) THEN
                           CALL LCOLOR ('WHITE')
                        ELSE
                           CALL LCOLOR ('YELOW')
                        ENDIF
                        CALL D2NODE (MXND, XN, YN, NTHERE,
     &                     NXL (2, LL))
                        CALL SFLUSH
                     ENDIF
                     NXL (1, LL) = NTHERE
                  ELSEIF (NXL (2, LL) .EQ. NGONE) THEN
                     IF ((VIDEO) .OR. (PGRAPH)) THEN
                        CALL LCOLOR ('BLACK')
                        CALL D2NODE (MXND, XN, YN, NXL (1, LL),
     &                     NXL (2, LL))
                        IF (GRAPH) THEN
                           CALL LCOLOR ('WHITE')
                        ELSE
                           CALL LCOLOR ('YELOW')
                        ENDIF
                        CALL D2NODE (MXND, XN, YN, NXL (1, LL),
     &                     NTHERE)
                        CALL SFLUSH
                     ENDIF
                     NXL (2, LL) = NTHERE
                  ENDIF
  170          CONTINUE

C  FIX LXN ARRAY
C  UNHOOK LOLD FROM NGONE AND FROM N1

               CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NGONE,
     &            LOLD, NNN, ERR, NOROOM)
               IF ((NOROOM) .OR. (ERR)) THEN
                  CALL MESAGE ('** PROBLEMS IN PINCH DELETING '//
     &               'LOLD FROM NGONE **')
                  GOTO 210
               ENDIF
               CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &            LOLD, NNN, ERR, NOROOM)
               IF ((NOROOM) .OR. (ERR)) THEN
                  CALL MESAGE ('** PROBLEMS IN PINCH DELETING '//
     &               'LOLD FROM N1 **')
                  GOTO 210
               ENDIF

C  ADD ALL LINES STILL HOOKED TO NGONE TO THE LIST OF LINES FOR NTHERE

               DO 180 II = 1, NL
                  LL = L1LIST (II)
                  IF (LL .NE. LOLD) THEN
                     CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &                  NTHERE, LL, NNN, ERR, NOROOM)
                     IF ((NOROOM) .OR. (ERR)) THEN
                        CALL MESAGE ('** PROBLEMS IN PINCH ADDING'//
     &                     'LL TO NTHERE **')
                        GOTO 210
                     ENDIF
                  ENDIF
  180          CONTINUE

C  DELETE NGONE (UNHOOK EVERYTHING FROM IT)

               DO 190 II = 1, 3
                  LXN (II, NGONE) = 0
  190          CONTINUE
               LXN (4, NGONE) = IAVAIL
               IAVAIL = NGONE
               NAVAIL = NAVAIL+1
               NUID (NGONE) = 0
               NLOOP = NLOOP - 2

C  PLOT THE CLOSURE BEFORE SMOOTHING

               IF (VIDEO) THEN
                  CALL SFLUSH
                  CALL SNAPIT (2)
               ENDIF

C  NOW SEE IF THE CLOSURE HAS PRODUCED A 2-LINE NODE AND
C  THUS REQUIRES THAT ONE OF THE ELEMENTS MUST BE SQUASHED

               IF ((LXN (3, N1) .EQ. 0) .AND. (LXN (2, N1) .GT. 0)) THEN
                  CALL DELEM (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &               NNN, NAVAIL, IAVAIL, N1, KXL (1, LXN (1, N1)),
     &               IDUM1, IDUM2, DDONE, CHECK, NOROOM, ERR)
                  IF ((NOROOM) .OR. (ERR)) GOTO 210
                  IF (VIDEO) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     CALL SNAPIT (2)
                  ENDIF
               ENDIF

C  SEE IF ONLY 2 NODES ARE LEFT ON THE LOOP

               IF (NLOOP .EQ. 2) THEN
                  CALL CLOSE2 (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL,
     &               NXL, LXN, LNODES, IAVAIL, NAVAIL, NNN, LLL, NTHERE,
     &               XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, PGRAPH, VIDEO,
     &               DEV1, KREG, NOROOM, ERR)
                  IF ((NOROOM) .OR. (ERR)) GOTO 210
                  NLOOP = 0
                  DONE = .TRUE.
                  NNN2 = NNN
               ELSE
                  NNN2 = 1
               ENDIF

C  PERFORM THE SMOOTH ON THE MESH

               CALL GETIME (TIME2)
               TIMEP = TIMEP + TIME2 - TIME1
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            N1, ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            NTHERE, ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, NTHERE), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, NTHERE)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, LNODES (3, NTHERE))), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, NTHERE), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, NTHERE)), ERR)
               IF (ERR) GOTO 210
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, LNODES (2, NTHERE))), ERR)
               IF (ERR) GOTO 210
               CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &            LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
               CALL GETIME (TIME1)
               IF ((PGRAPH) .OR. (VIDEO)) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  IF (VIDEO) CALL SNAPIT (2)
               ENDIF

C  CALCULATE NEW ANGLES

               IF (NLOOP .GT. 0)
     &            CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL,
     &            LXN, NLOOP, ANGLE, LNODES, NTHERE, LLL, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
               IF (ERR) GOTO 210

               IF ((NODE .EQ. NGONE) .OR. (NODE .EQ. N1))
     &            NODE = NTHERE
               IF ((NADJ1 .EQ. NGONE) .OR. (NADJ1 .EQ. N1))
     &            NADJ1 = NTHERE
               IF ((NADJ2 .EQ. NGONE) .OR. (NADJ2 .EQ. N1))
     &            NADJ2 = NTHERE
               IF (DONE) GOTO 200

            ENDIF
         ENDIF
         GOTO 120
      ENDIF

  200 CONTINUE

C  NOW GO BACK AND GET THE NEW CORNERS AND TRY AGAIN IF THE FIRST
C  TIME WAS SUCCESSFUL

      IF ((FOUND) .AND. (.NOT. DONE)) THEN
         FOUND = .FALSE.
         GOTO 100
      ENDIF

C  NOW PLOT THE NEW BOUNDARY IF A PINCH HAS OCCURRED

      IF ((FOUND) .AND. (GRAPH) .AND. (.NOT. PGRAPH)) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, LLL, DEV1, KREG)
      ENDIF

C  BE SURE THAT THE LXN ARRAY WILL GET FIXED (FIXLXN) LATER UP
C  TO THE CURRENT NNN

      NNNOLD = NNN

  210 CONTINUE

      CALL GETIME (TIME2)
      TIMEP = TIMEP + TIME2 - TIME1
      RETURN

      END
