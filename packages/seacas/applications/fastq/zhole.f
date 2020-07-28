C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ZHOLE (MP, ML, MS, MR, NS, MAXNL, MAXNP, MAXPRM, NPRM,
     &   MAXNBC, MAXSBC, KNBC, KSBC, KNUM, IPOINT, COOR, IPBOUN, ILINE,
     &   LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, ISLIST, INDXH, NPPF, IFPB, LISTPB, NLPF, IFLB,
     &   LISTLB, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKPB,
     &   LINKLB, LINKSB, X, Y, NID, LISTL, MARKED, NL, LSTNBC, MXND,
     &   XN, YN, NUID, LXK, KXL, NXL, LXN, NXH, NPERIM, NNN, NNNOLD,
     &   KKK, LLL, IAVAIL, NAVAIL, JHOLE, INSIDE, EPS, NOROOM, ERR,
     &   AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR,
     &   MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &   REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE ZHOLE  =  REMESHES AROUND HOLE IN REGION

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), NINT(ML), LTYPE(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION ISLIST(4*MR), LINKP(2, MP), LINKL(2, ML)
      DIMENSION LINKS(2, MS), LISTL(MAXNL)
      DIMENSION X(MAXNP), Y(MAXNP), NID(MAXNP, MAXPRM), NPERIM(MAXPRM)
      DIMENSION LINKPB(2, MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION LINKLB(2, ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION LINKSB(2, ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION LSTNBC(MAXNBC)
      DIMENSION XN(MXND), YN(MXND), NUID(MXND), LXK(4, MXND)
      DIMENSION KXL(2, 3*MXND), NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION NXH(MXND)
      DIMENSION KLIST1(20), LINES(20), NODES(4)

      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL ADDLNK, CCW, DELETE, ERR, EVEN, LREAL, NOROOM, COUNT
      LOGICAL LPNTIN, REMESH, LCIRCL, LDEL

C  CHECK FOR INPUT ERRORS

      ERR = .FALSE.
      IF (NNN - NNNOLD .LE. 0) THEN
         CALL MESAGE ('NO NODES DEFINED IN REGION')
         ERR = .TRUE.

C  GOOD INPUT

      ELSE
         LNUM = ABS(ISLIST(INDXH))
         ADDLNK = .FALSE.
         CALL LTSORT (ML, LINKL, LNUM, LIN, ADDLNK)
         LCIRCL = NS .EQ. 1 .AND. LTYPE(LIN) .EQ. 3

C  CIRCULAR HOLE

         IF (LCIRCL) THEN
            CALL LTSORT (MP, LINKP, LCON(1, LIN), I1, ADDLNK)
            CALL LTSORT (MP, LINKP, LCON(2, LIN), I2, ADDLNK)
            IF (I1 .NE. I2) THEN
               CALL MESAGE ('CIRCULAR HOLE DOES NOT CLOSE')
               ERR = .TRUE.
               GO TO 380
            END IF
            CALL LTSORT (MP, LINKP, LCON(3, LIN), I2, ADDLNK)
            XCEN = COOR(1, I2)
            YCEN = COOR(2, I2)
            RADIUS = (COOR(1, I1) - XCEN)**2 + (COOR(2, I1) - YCEN)**2
            IF (RADIUS .LE. 0.0) THEN
               CALL MESAGE ('RADIUS HAS ZERO LENGTH')
               ERR = .TRUE.
               GO TO 380
            END IF
            XMIN = XCEN - SQRT(RADIUS)
            XMAX = XCEN + SQRT(RADIUS)
            YMIN = YCEN - SQRT(RADIUS)
            YMAX = YCEN + SQRT(RADIUS)
            NPERV = NINT(LIN)

C  NON-CIRCULAR HOLE
         ELSE
            NLP1 = NL + 1
            CCW = .TRUE.
            COUNT = .FALSE.
            EVEN = .FALSE.
            LREAL = .FALSE.

            CALL PERIM (MP, ML, MS, NS, MAXNL, MAXNP, MAXNBC, MAXSBC,
     &         KNBC, KSBC, KNUM, IPOINT, COOR, IPBOUN, ILINE, LTYPE,
     &         NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &         ILLIST, ISLIST(INDXH), NPPF, IFPB, LISTPB, NLPF, IFLB,
     &         LISTLB, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKPB,
     &         LINKLB, LINKSB, X, Y, NID(1, NPRM), NPERIM(NPRM),
     &         LISTL(NLP1), NL1, LSTNBC, MARKED, EVEN, LREAL, ERR, CCW,
     &         COUNT, NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD,
     &         LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD,
     &         NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS,
     &         SIZMIN, EMAX, EMIN)
            IF (ERR .OR. NOROOM) GO TO 380
            NPERV = NPERIM(NPRM)
            NPNT = NPERV
            XMIN = X(1)
            XMAX = XMIN
            YMIN = Y(1)
            YMAX = YMIN
            DO 100 I = 1, NPNT
               XCEN = XCEN + X(I)
               YCEN = YCEN + Y(I)
               XMIN = MIN(XMIN, X(I))
               XMAX = MAX(XMAX, X(I))
               YMIN = MIN(YMIN, Y(I))
               YMAX = MAX(YMAX, Y(I))
  100       CONTINUE
            XCEN = XCEN/DBLE(NPNT)
            YCEN = YCEN/DBLE(NPNT)
            RADIUS = SQRT((XCEN - X(1))**2 + (YCEN - Y(1))**2)
            DO 110 I = 2, NPNT
               R = SQRT((XCEN - X(I))**2 + (YCEN - Y(I))**2)
               RADIUS = MIN(RADIUS, R)
  110       CONTINUE
         END IF

C  INITIALIZE NODES PER (ON) HOLE

         DO 120 I = 1, NNN
            NXH(I) = 0
  120    CONTINUE

C  DELETE EVERYTHING ATTACHED TO NODES WITHIN HOLE

         NEAR = 0
         SMALL = 0.0
         DO 130 I = NNNOLD + 1, NNN
            IF (XN(I) .GT. XMIN .AND. XN(I) .LT. XMAX .AND.
     &         YN(I) .GT. YMIN .AND. YN(I) .LT. YMAX) THEN
               DIST = (XN(I) - XCEN)**2 + (YN(I) - YCEN)**2
               IF (LCIRCL) THEN
                  LDEL = DIST .LT. RADIUS
               ELSE
                  LDEL = LPNTIN (MAXNP, X, Y, NPNT, XN(I), YN(I))
               END IF
               IF (DIST .LT. SMALL .OR. NEAR .EQ. 0) THEN
                  NEAR = I
                  SMALL = DIST
               END IF
               IF (LDEL) THEN
                  IF (NUID(I) .EQ. 0) THEN
                     CALL DELHOL (I, MXND, LXK, KXL, NXL, LXN, NXH,
     &                  NUID, NNN, IAVAIL, NAVAIL, NOROOM, ERR)
                     IF (NOROOM .OR. ERR) GO TO 380

C  CANNOT DELETE BOUNDARY NODES

                  ELSE
                     CALL MESAGE ('HOLE CROSSES FIXED BOUNDARY')
                     ERR = .TRUE.
                     GO TO 380
                  END IF
               END IF
            END IF
  130    CONTINUE

C  PROCESS SMALL CIRCLES (I.E. SMALLER THAN AN ELEMENT)

         IF (SMALL .GT. RADIUS) THEN
            CCW = .TRUE.
            CALL GKXN (MXND, KXL, LXN, NEAR, KS1, KLIST1, ERR)
            DO 150 I = 1, KS1
               CALL GNXKA (MXND, XN, YN, KLIST1(I), NODES, AREA, LXK,
     &            NXL, CCW)
               SUM = 0.0
               DO 140 J = 1, 4
                  J1 = J + 1
                  IF (J1 .GT. 4) J1 = 1
                  SUM = SUM + ABS((XN(NODES(J1)) - XN(NODES(J)))*(YCEN
     &               - YN(NODES(J))) - (XCEN - XN(NODES(J)))
     &               *(YN(NODES(J1)) - YN(NODES(J))))
  140          CONTINUE
               SUM = SUM/2.0
               IF (ABS((AREA - SUM)/AREA) .LT. 1.0E-4) GO TO 160
  150       CONTINUE
            CALL MESAGE ('FAILED TO FIND ELEMENT SURROUNDING '//
     &         'SMALL HOLE')
            ERR = .TRUE.
            GO TO 380
  160       CONTINUE

            DO 170 I = 1, 4
               IF (NUID(NODES(I)) .EQ. 0) THEN
                  CALL DELHOL (NODES(I), MXND, LXK, KXL, NXL, LXN, NXH,
     &               NUID, NNN, IAVAIL, NAVAIL, NOROOM, ERR)
                  IF (NOROOM .OR. ERR) GO TO 380
               END IF
  170       CONTINUE
         END IF

C  SQUARE UP BOUNDARY (DELETE INTERIOR NODES WITH ONLY TWO LINES)

  180    CONTINUE
         DELETE = .FALSE.
         DO 190 I = NNNOLD + 1, NNN
            IF (NXH(I) .EQ. 1 .AND. NUID(I) .EQ. 0) THEN
               CALL GETLXN (MXND, LXN, I, LINES, NUML, ERR)
               IF (NUML .EQ. 2 .AND. NUID(I) .EQ. 0) THEN
                  CALL DELHOL (I, MXND, LXK, KXL, NXL, LXN, NXH,
     &               NUID, NNN, IAVAIL, NAVAIL, NOROOM, ERR)
                  IF (NOROOM .OR. ERR) GO TO 380
                  DELETE = .TRUE.
               END IF
            END IF
  190    CONTINUE
         IF (DELETE) GO TO 180

C  GENERATE DELETED ELEMENT BOUNDARY NODE LIST

         NH = 0
         DO 200 I = NNNOLD + 1, NNN
            IF (NXH(I) .GT. 0) THEN
               NH = NH + 1
               NXH(NH) = I
            END IF
  200    CONTINUE

C  ENSURE THAT THERE ARE A MINIMUM OF MIN(12, NPERV) INTERVALS
C     AROUND HOLE

         IF (NH .LT. MAX(12, NPERV)) THEN
            DO 210 I = NH + 1, MXND
               NXH(I) = 0
  210       CONTINUE

            DO 220 I = 1, NH
               IF (NUID(NXH(I)) .EQ. 0) THEN
                  CALL DELHOL (NXH(I), MXND - NH, LXK, KXL, NXL, LXN,
     &               NXH(NH + 1), NUID, NNN, IAVAIL, NAVAIL, NOROOM,
     &               ERR)
                  IF (NOROOM .OR. ERR) GO TO 380
                  DELETE = .TRUE.
               ELSE
                  I1 = NXH(I) + NH
                  NXH(I1) = 1
               END IF
  220       CONTINUE

            IF (DELETE) THEN
               I1 = 0
               DO 230 I = NH + 1, MXND
                  I1 = I1 + 1
                  NXH(I1) = NXH(I)
  230          CONTINUE
               GO TO 180
            ELSE
               CALL MESAGE
     &            ('INTERVAL MISMATCH BETWEEN HOLE AND BOUNDARY.')
               CALL MESAGE ('ALL ELEMENTS DELETED.')
               ERR = .TRUE.
               GO TO 380
            END IF
         END IF

C  ORDER THE INTERIOR NODE LIST

         DO 260 I = 1, NH - 1
            CALL GETLXN (MXND, LXN, NXH(I), LINES, NUML, ERR)
            DO 250 J = 1, NUML
               J1 = NXL(2, LINES(J)) + NXL(1, LINES(J)) - NXH(I)
               DO 240 K = I + 1, NH
                  IF (NXH(K) .EQ. J1) THEN
                     NXH(K) = NXH(I + 1)
                     NXH(I + 1) = J1
                     GO TO 260
                  END IF
  240          CONTINUE
  250       CONTINUE
  260    CONTINUE

C  MAKE SURE LOOP CLOSES

         CALL GETLXN (MXND, LXN, NXH(NH), LINES, NUML, ERR)
         DO 270 J = 1, NUML
            J1 = NXL(2, LINES(J)) + NXL(1, LINES(J)) - NXH(NH)
            IF (NXH(1) .EQ. J1) THEN
               GO TO 280
            END IF
  270    CONTINUE
         CALL MESAGE ('HOLE PERIMETER DOES NOT CLOSE')
         ERR = .TRUE.
         GO TO 380
  280    CONTINUE

C  MAKE SURE HOLE PERIMETER IS DEFINED COUNTER-CLOCKWISE

         PI = ACOS(-1.0)
         TWOPI = PI + PI
         SPIRO = 0.0
         AGOLD = ATAN2(YN(NXH(1)) - YN(NXH(NH)),
     &      XN(NXH(1)) - XN(NXH(NH)))
         DO 290 I = 1, NH
            IF (I .EQ. NH) THEN
               NEXT = 1
            ELSE
               NEXT = I + 1
            END IF
            AGNEW = ATAN2(YN(NXH(NEXT)) - YN(NXH(I)),
     &         XN(NXH(NEXT)) - XN(NXH(I)))
            DIFF = AGNEW - AGOLD
            IF (DIFF .GT. PI) THEN
               DIFF = DIFF - TWOPI
            ELSE IF (DIFF .LT .-PI) THEN
               DIFF = DIFF + TWOPI
            END IF
            IF (ABS(ABS(DIFF) - PI) .LT .1.0E-3) THEN
               CALL MESAGE ('PERIMETER CONTAINS SWITCHBACKS')
               ERR = .TRUE.
               GO TO 380
            ENDIF
            SPIRO = SPIRO + DIFF
            AGOLD = AGNEW
  290    CONTINUE

         IF (SPIRO .LT .0.0) THEN
            DO 300 I = 1, NH/2
               ITEMP = NXH(I)
               NXH(I) = NXH(NH - I + 1)
               NXH(NH - I + 1) = ITEMP
  300       CONTINUE
         ELSE IF ((ABS(SPIRO) .LT .PI) .OR.
     &      (ABS(SPIRO) .GT. (3.*PI))) THEN
            CALL MESAGE
     &         ('UNABLE TO DETERMINE CW OR CCW SENSE OF HOLE')
            ERR = .TRUE.
            GO TO 380
         ENDIF

C  FIND THE BEST STARTING POINT ON THE CIRCULAR HOLE

         IF (NNN + NH .GT. MXND) THEN
            NOROOM = .TRUE.
            GO TO 380
         END IF

C  GENERATE THE PERIMETER OF THE HOLE

         EVEN = .TRUE.
         CCW = .TRUE.
         LREAL = .TRUE.
         COUNT = .TRUE.
         IF (NPERV .NE. NH) THEN
            IF (LCIRCL) THEN
               NINT(LIN) = NH
            ELSE
               TDIST = 0.0
               DO 310 I = 1, NL1
                  LNUM = LISTL(NL + I)
                  CALL LTSORT (ML, LINKL, LNUM, LIN, ADDLNK)
                  I1 = LCON (1, LIN)
                  I2 = LCON (2, LIN)
                  I3 = LCON (3, LIN)
                  CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
                  CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
                  IF (I3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
                     IF (I3 .LT. 0) J3 = -J3
                  ELSE
                     J3 = 0
                  END IF
                  CALL LINLEN (MP, COOR, LINKP, JHOLE, ILINE(LIN),
     &               LTYPE(LIN), I3, J1, J2, J3, DIST, ERR)
                  IF (ERR) GO TO 380
                  TDIST = TDIST + DIST
  310          CONTINUE
               NSUM = 0
               DO 320 I = 1, NL1 - 1
                  LNUM = LISTL(NL + I)
                  CALL LTSORT (ML, LINKL, LNUM, LIN, ADDLNK)
                  I1 = LCON (1, LIN)
                  I2 = LCON (2, LIN)
                  I3 = LCON (3, LIN)
                  CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
                  CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
                  IF (I3 .NE. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
                     IF (I3 .LT. 0) J3 = -J3
                  ELSE
                     J3 = 0
                  END IF
                  CALL LINLEN (MP, COOR, LINKP, JHOLE, ILINE(LIN),
     &               LTYPE(LIN), I3, J1, J2, J3, DIST, ERR)
                  IF (ERR) GO TO 380
                  NINT(LIN) = INT(NH * DIST/TDIST + 0.5)
                  NSUM = NSUM + NINT(LIN)
  320          CONTINUE
               LNUM = LISTL(NL + NL1)
               CALL LTSORT (ML, LINKL, LNUM, LIN, ADDLNK)
               NINT(LIN) = NH - NSUM
            END IF
            CALL MESAGE ('INTERVALS MODIFIED FOR HOLE')
         END IF
         NLP1 = NL + 1
         CALL PERIM (MP, ML, MS, NS, MAXNL, MAXNP, MAXNBC, MAXSBC, KNBC,
     &      KSBC, KNUM, IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &      FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST,
     &      ISLIST(INDXH), NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF,
     &      IFSB, LISTSB, LINKP, LINKL, LINKS, LINKPB, LINKLB, LINKSB,
     &      X, Y, NID(1, NPRM), NPERIM(NPRM), LISTL(NLP1), NL1, LSTNBC,
     &      MARKED, EVEN, LREAL, ERR, CCW, COUNT, NOROOM, AMESUR, XNOLD,
     &      YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK,
     &      NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
         IF (NOROOM .OR. ERR) GO TO 380
         IF (NL1 .GE. 0 .AND. NL + NL1 .LE. MAXNL) THEN
            NL = NL + NL1
         ELSE
            CALL MESAGE ('UNABLE TO ADD HOLE LINES TO REGION LINE LIST')
            ERR = .TRUE.
            GO TO 380
         END IF

C  TACK THE HOLE LINE LIST ONTO THE BOUNDARY LINE LIST

         IF (NPERIM(NPRM) .NE. NH) THEN
            CALL MESAGE ('INTERVAL MISMATCH ON HOLE PERIMETER')
            ERR = .TRUE.
            GO TO 380
         END IF

         ISTART = 0
         DIST = 0.0
         DO 340 I = 1, NH
            SUM = 0.0
            I1 = I - 1
            DO 330 J = 1, NH
               I1 = I1 + 1
               IF (I1 .GT. NH) I1 = 1
               I2 = NXH(I1)
               SUM = SUM + (XN(I2) - X(J))**2 + (YN(I2) - Y(J))**2
  330       CONTINUE

            IF (SUM .LT. DIST .OR. ISTART .EQ. 0) THEN
               DIST = SUM
               ISTART = I
            END IF
  340    CONTINUE

         NNNX = NNN
         DO 350 J = 1, NH
            NNN = NNN + 1
            XN(NNN) = X(J)
            YN(NNN) = Y(J)
            NUID(NNN) = NID(J, NPRM)
  350    CONTINUE

C  FIRST ROW OF ELEMENTS

         CALL INNERH (MXND, NXH, NUID, LXK, KXL, NXL, LXN, KKK, LLL,
     &      NNN, NNNX, NH, ISTART, IAVAIL, NAVAIL, NOROOM, ERR)
         IF (NOROOM .OR. ERR) GO TO 380

C  INSERT INNER NECKLACE OF ELEMENTS

         ISTART = 1
         DO 370 J = 1, INSIDE
            NNNX = NNN
            DO 360 I = 1, NH
               NNN = NNN + 1
               NXH(I) = NNNX - NH + I
               XN(NNN) = XN(NXH(I))
               YN(NNN) = YN(NXH(I))
               NUID(NNN) = NUID(NXH(I))
               NUID(NXH(I)) = 0
               LXN(2, NXH(I)) = ABS(LXN(2, NXH(I)))
  360       CONTINUE
            CALL INNERH (MXND, NXH, NUID, LXK, KXL, NXL, LXN, KKK, LLL,
     &         NNN, NNNX, NH, ISTART, IAVAIL, NAVAIL, NOROOM, ERR)
            IF (NOROOM .OR. ERR) GO TO 380
  370    CONTINUE
      END IF

  380 CONTINUE
      RETURN
      END
