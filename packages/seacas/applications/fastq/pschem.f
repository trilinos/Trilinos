C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C**               SIDE BOUNDARY FLAG WOULD REDIMENSION FOREVER.
C**               THE KKSBC VARIABLE WAS CHANGED TO BE SET AT THE
C**               BEGINNING OF THE ROUTINE INSTEAD OF RIGHT BEFORE THE
C**               ZHOLE PROCESSING WAS STARTED.

      SUBROUTINE PSCHEM (MP, ML, MS, MR, N, IPOINT, COOR, IPBOUN, ILINE,
     &   LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &   ILLIST, IREGN, NSPR, IFSIDE, ISLIST, NPPF, IFPB, LISTPB, NLPF,
     &   IFLB, LISTLB, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS, LINKR,
     &   LINKSC, LINKPB, LINKLB, LINKSB, IFHOLE, NHPR, IHLIST, MAXNBC,
     &   KNBC, MAXSBC, KSBC, MXND, NNN, NNNOLD, KKK, KKKOLD, LLL, X, Y,
     &   NID, LISTL, XN, YN, NUID, LXK, KXL, NXL, LXN, LSTNBC, NXH,
     &   NPERIM, MARKED, IAVAIL, NAVAIL, MXNL, MXNPER, NPER, MAXPRM,
     &   NPRM, MSC, ISCHM, SCHEME, SCHSTR, RECT, M1, INSIDE, JJHOLE,
     &   KKSBC, DEV1, EIGHT, NINE, STEP, L, NL, MCOM, CIN, IIN, RIN,
     &   KIN, ICOM, JCOM, XMIN, XMAX, YMIN, YMAX, ICODE, NOROOM, ERR,
     &   AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR,
     &   MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &   REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  PSCHEM - PROCESS A COMPLETE SCHEME

C***********************************************************************

      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IREGN(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION ISCHM(MSC), SCHEME(MSC)
      DIMENSION NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR), LINKPB(2, MP)
      DIMENSION LINKLB(2, ML), LINKSB(2, ML)
      DIMENSION N(29), CIN(MCOM), IIN(MCOM), RIN(MCOM), KIN(MCOM)
      DIMENSION X(MXNPER), Y(MXNPER), NID(MXNPER*MAXPRM)
      DIMENSION LISTL(MXNL), MARKED(3, MXNL)
      DIMENSION XN(MXND), YN(MXND), NUID(MXND), LXK(4, MXND)
      DIMENSION KXL(2, 3*MXND), NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION LSTNBC(MAXNBC), NXH(MXND), NPERIM(MAXPRM)

      DIMENSION ILPC(10)

      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(2 * NPEOLD), BMESUR(NPNOLD)

      CHARACTER*72 ADDSTR, CIN, DEFSCH, SCHEME, SCHOLE, SCHSTR
      CHARACTER DEV1*3

      LOGICAL DOLINK, EIGHT, ERR, NINE, NOROOM, STEP
      LOGICAL ADDLNK, ACTIVE, IANS, DONE, DOSMOO, DOTILT, LACT(10)
      LOGICAL RECT, REMESH

      DATA IEXIT, IOVER, IQUIT /1, 2, 3/

      ALPHA = 0.7
      ASMALL = 45.0
      RO = 1.
      TOL = .03
      WF = 0.7
      ICODE = 0
      ISTYPE = 1
      J = 1
      NACALL = 0
      NEWSGN = 0
      DOTILT = .TRUE.
      DOSMOO = .TRUE.
      ERR = .FALSE.
      NOROOM = .FALSE.
      KKSBC = KSBC
      NLEFTP = 0

      NIT = 5 * NPER/2
      CALL MNORM (MXND, XN, YN, NXL, LLL, STDLEN)
      EPS = 0.03 * STDLEN

      CALL STRIPB (SCHSTR, I1, LENSCH)
  100 CONTINUE

      IISIGN = NEWSGN
      NEWSGN = 0

C  ACT ON NEXT COMMAND

C  A - ALPHA CONTROL FOR APALSM

      IF ((SCHSTR(J:J) .EQ. 'A') .OR. (SCHSTR(J:J) .EQ. 'a')) THEN
         IF (IISIGN .GE. 0) THEN
            ALPHA = MIN(ALPHA + 0.1, 1.0)
         ELSE
            ALPHA = MAX(ALPHA - 0.1, 0.0)
         END IF
         IF (ISTYPE .LE. 2) DOSMOO = .TRUE.
         WRITE(*, 10010) ALPHA

C  B - BLACKER TRANSITION REGION TEST

      ELSE IF ((SCHSTR(J:J) .EQ. 'B') .OR.
     &   (SCHSTR(J:J) .EQ. 'b')) THEN
         CONTINUE

C  C - SEMI-CIRCLE REGION TEST

      ELSE IF ((SCHSTR(J:J) .EQ. 'C') .OR.
     &   (SCHSTR(J:J) .EQ. 'c')) THEN
         CONTINUE

C  D - DELETE WORST RHOMBUS

      ELSE IF ((SCHSTR(J:J) .EQ. 'D') .OR.
     &   (SCHSTR(J:J) .EQ. 'd')) THEN
         LIMIT = 0
         IF (DOTILT) THEN
            CALL RESTA (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &         KKKOLD, NAVAIL, IAVAIL, NNN, LIMIT, IREST, TILT, ERR,
     &         NOROOM)
            IF (NOROOM) THEN
               GO TO 140
            ELSE IF (ERR) THEN
               CALL MESAGE ('ERROR DURING SHAPE SORTING OF ELEMENTS')
               GO TO 140
            END IF
            DOTILT = .FALSE.
         END IF
         ATILT = MIN(45.0, TILT*0.667)
         ATILT = (ASMALL/45.0)*ATILT
         CALL SQUASH (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &      KKKOLD, NNN, NAVAIL, IAVAIL, ATILT, DONE, NOROOM, ERR)
         IF (NOROOM) THEN
            GO TO 140
         ELSE IF (ERR) THEN
            CALL MESAGE ('ERROR DURING DELETION OF ELEMENT')
            GO TO 140
         END IF
         IF (DONE) THEN
            DOSMOO = .TRUE.
            RECT = .FALSE.
            ACTIVE = .TRUE.
            CALL MESAGE ('ELEMENT DELETED')
         ELSE
            CALL MESAGE ('NO ELEMENT(S) DELETED')
         END IF

C  E - EXIT SAVING SCHEME AND REGION

      ELSE IF ((SCHSTR(J:J) .EQ. 'E') .OR.
     &   (SCHSTR(J:J) .EQ. 'e')) THEN

C  SAVE THE SCHEME USED IF STEPPING THROUGH

         IF (STEP) THEN
            DOLINK = .TRUE.
            NOLD10 = N(10)
            NOLD24 = N(24)
            JJ = ABS(IREGN(L))
            CALL STRLNG (SCHSTR, LTRY)
            SCHSTR(LTRY:LTRY) = ' '
            CALL INSCHM (MR, MSC, N(10), N(24), JJ, SCHSTR, ISCHM,
     &         SCHEME, LINKSC, DEFSCH, NOROOM, DOLINK)
            IF (NOROOM) THEN
               N(10) = NOLD10
               N(24) = NOLD24
               WRITE(*, 10020) SCHSTR(1:LENSCH), ABS(IREGN(L))
               GO TO 140
            END IF
         END IF
         ICODE = IEXIT
         GO TO 140

C  F - CONTROL UNDER- OR OVER-RELAXATION

      ELSE IF ((SCHSTR(J:J) .EQ. 'F') .OR.
     &   (SCHSTR(J:J) .EQ. 'f')) THEN
         RODEL = .25
         IF (IISIGN .GE. 0) THEN
            RO = RO + RODEL
         ELSE
            RO = MAX(RO - RODEL, RODEL)
         END IF
         DOSMOO = .TRUE.
         WRITE(*, 10030) RO

C  H - INDICATES A HELP MESSAGE RESPONSE IF STEPPING

      ELSE IF ((SCHSTR(J:J) .EQ. 'H') .OR.
     &   (SCHSTR(J:J) .EQ. 'h')) THEN
         IF (STEP) CALL HELP_FQ (3)

C  I - CHANGE MAX SMOOTHING ITERATIONS

      ELSE IF ((SCHSTR(J:J) .EQ. 'I') .OR.
     &   (SCHSTR(J:J) .EQ. 'i')) THEN
         IF (IISIGN .GE. 0) THEN
            NIT = INT(DBLE(NIT)*1.500 + 0.51)
         ELSE
            NIT = INT(DBLE(NIT)*.6667 + 0.51)
         END IF
         DOSMOO = .TRUE.
         WRITE(*, 10050) NIT

C  J - CONTROL SMOOTHING TOLERANCE

      ELSE IF ((SCHSTR(J:J) .EQ. 'J') .OR.
     &   (SCHSTR(J:J) .EQ. 'j')) THEN
         IF (IISIGN .GE. 0) THEN
            TOL = TOL*1.259921
            EPS = EPS*1.259921
         ELSE
            TOL = TOL/1.259921
            EPS = EPS/1.259921
         END IF
         DOSMOO = .TRUE.
         WRITE(*, 10060) TOL, EPS

C  L - INSERT ROW OF ELEMENTS AROUND A HOLE (TOO LATE NOW)

      ELSE IF ((SCHSTR(J:J) .EQ. 'L') .OR.
     &   (SCHSTR(J:J) .EQ. 'l')) THEN
         CONTINUE

C  M - LOGICAL MESH SIDES CHOSEN BY QMESH (TOO LATE NOW)

      ELSE IF ((SCHSTR(J:J) .EQ. 'M') .OR.
     &   (SCHSTR(J:J) .EQ. 'm')) THEN
         CONTINUE

C  N - NECKLACE

      ELSE IF ((SCHSTR(J:J) .EQ. 'N') .OR.
     &   (SCHSTR(J:J) .EQ. 'n')) THEN
         CALL NCKLCE (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK, NNN,
     &      NNNOLD, LLL, NAVAIL, IAVAIL, EPS, NOROOM, ERR)
         IF (NOROOM) THEN
            GO TO 140
         ELSE IF (ERR) THEN
            CALL MESAGE ('ERROR DURING NECKLACING OF REGION')
            GO TO 140
         END IF
         ACTIVE = .TRUE.
         RECT = .FALSE.
         DOSMOO = .TRUE.
         CALL MESAGE ('NECKLACE INSTALLED')

C  O - ORIGINATE THE MESH AGAIN.

      ELSE IF ((SCHSTR(J:J) .EQ. 'O') .OR.
     &   (SCHSTR(J:J) .EQ. 'o')) THEN
         CALL MESAGE ('PROCESSING RETURNED TO ORIGINAL')
         SCHSTR = ' '
         ICODE = IOVER
         GO TO 140

C  P - PLOT

      ELSE IF ((SCHSTR(J:J) .EQ. 'P') .OR.
     &   (SCHSTR(J:J) .EQ. 'p')) THEN
         IF (STEP) THEN
            CALL PLOTL (MXND, XN, YN, NXL, ABS(IREGN(L)), XMIN, XMAX,
     &         YMIN, YMAX, LLL, DEV1)
         ELSE
            CALL MESAGE ('PLOTTING AVAILABLE ONLY IN INTERACTIVE '//
     &         'STEP PROCESSING')
         END IF

C  Q - QUIT STEP PROCESSING WITHOUT SAVING MESH

      ELSE IF ((SCHSTR(J:J) .EQ. 'Q') .OR.
     &   (SCHSTR(J:J) .EQ. 'q')) THEN
         CALL MESAGE  ('REGION PROCESSING ABORTED WITH "QUIT"')
         ICODE = IQUIT
         GO TO 140

C  R - RESTRUCTURE

      ELSE IF ((SCHSTR(J:J) .EQ. 'R') .OR.
     &   (SCHSTR(J:J) .EQ. 'r')) THEN
         LIMIT = MXND
         CALL RESTA (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &      KKKOLD, NAVAIL, IAVAIL, NNN, LIMIT, IREST, TILT, ERR,
     &      NOROOM)
         IF (NOROOM) THEN
            GO TO 140
         ELSE IF (ERR) THEN
            CALL MESAGE ('ERROR DURING RESTRUCTURE OF REGION')
            GO TO 140
         END IF
         DOTILT = .FALSE.
         IF (IREST .GE. 1) THEN
            DOSMOO = .TRUE.
            RECT = .FALSE.
            ACTIVE = .TRUE.
            CALL MESAGE ('RESTRUCTURE COMPLETED')
         END IF

C  S - SMOOTH

      ELSE IF ((SCHSTR(J:J) .EQ. 'S') .OR.
     &   (SCHSTR(J:J) .EQ. 's')) THEN
         IF (DOSMOO) THEN
            IF (ISTYPE .EQ. 1) THEN
               IF (RECT) THEN
                  CALL REPSMO (MXND, XN, YN, LXN, NNN, NNNOLD, NIT, EPS,
     &               RO, M1)
                  CALL MESAGE ('EQUIPOTENTIAL SMOOTHING COMPLETED')
               ELSE
                  IF ((NACALL/5)*5 .EQ. NACALL) THEN
                     CALL ARELAX (MXND, XN, YN, LXK, KXL, NXL, LLL,
     &                  ARFACT)
                  END IF
                  CALL APALSM (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN,
     &               NNNOLD, NIT, TOL, RO*ARFACT, ALPHA, ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('ERROR DURING AREA PULL & '//
     &                  'LAPLACIAN SMOOTHING')
                     GO TO 140
                  END IF
                  CALL MESAGE ('AREA PULL AND LAPLACIAN SMOOTHING '//
     &               'COMPLETED')
                  NACALL = NACALL + 1
               END IF
            ELSE IF (ISTYPE .EQ. 2) THEN
               IF ((NACALL/5)*5 .EQ. NACALL) THEN
                  CALL ARELAX (MXND, XN, YN, LXK, KXL, NXL, LLL, ARFACT)
               END IF
               CALL APALSM (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN,
     &            NNNOLD, NIT, TOL, RO*ARFACT, ALPHA, ERR)
               IF (ERR) THEN
                  CALL MESAGE
     &               ('ERROR DURING AREA PULL & LAPLACIAN SMOOTHING')
                  GO TO 140
               END IF
               CALL MESAGE
     &            ('AREA PULL AND LAPLACIAN SMOOTHING COMPLETED')
               NACALL = NACALL + 1
            ELSE IF (ISTYPE .EQ. 3) THEN
               CALL CIAPAL (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN,
     &            NNNOLD, NIT, EPS, RO, 0.5)
               CALL MESAGE ('CENTROID INVERSE PUSH AND LAPLACIAN '//
     &            'SMOOTHING COMPLETED')
            ELSE IF (ISTYPE .EQ. 4) THEN
               CALL CASMO (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN,
     &            NNNOLD, NIT, EPS, RO)
               CALL MESAGE
     &            ('CENTROID AREA PULL SMOOTHING COMPLETED')
            ELSE IF (ISTYPE .EQ. 5) THEN
               IF (RECT) THEN
                  CALL REPSMO (MXND, XN, YN, LXN, NNN, NNNOLD, NIT, EPS,
     &               RO, M1)
                  CALL MESAGE ('EQUIPOTENTIAL SMOOTHING COMPLETED')
               ELSE
                  CALL SMOGS (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT,
     &               EPS, RO)
                  CALL MESAGE ('LAPLACIAN SMOOTHING COMPLETED')
               END IF
            ELSE IF (ISTYPE .EQ. 6) THEN
               CALL L2SMO (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT,
     &            EPS, RO)
               CALL MESAGE ('LENGTH WEIGHTED LAPLACIAN SMOOTHING '//
     &            'COMPLETED')
            ELSE IF (ISTYPE .EQ. 7) THEN
               CALL ISOLAP (MXND, XN, YN, LXK, KXL, NXL, LXN, NNN,
     &            NNNOLD, WF, NIT, EPS, RO)
               CALL MESAGE
     &            ('LAPLACIAN-ISOPARAMETRIC SMOOTHING COMPLETED')
            END IF
            DOSMOO = .FALSE.
            ACTIVE = .TRUE.
         ELSE
            CALL MESAGE ('MESH AND/OR SMOOTHING PARAMETERS HAVE')
            CALL MESAGE ('NOT CHANGED - NO SMOOTHING ATTEMPTED')
         END IF

C  T - TRIANGULAR REGION - QUAD MESH GENERATION

      ELSE IF ((SCHSTR(J:J) .EQ. 'T') .OR.
     &   (SCHSTR(J:J) .EQ. 't')) THEN
         CONTINUE

C  U - PENTAGON REGION - QUAD MESH GENERATION

      ELSE IF ((SCHSTR(J:J) .EQ. 'U') .OR.
     &   (SCHSTR(J:J) .EQ. 'u')) THEN
         CONTINUE

C  V - CHANGE ASMALL FOR SQUASH (D)

      ELSE IF ((SCHSTR(J:J) .EQ. 'V') .OR.
     &   (SCHSTR(J:J) .EQ. 'v')) THEN
         IF (IISIGN .GE. 0) THEN
            ASMALL = MIN(ASMALL + 2.5, 80.0)
         ELSE
            ASMALL = MAX(ASMALL - 2.5, 10.0)
         END IF
         WRITE(*, 10070) ASMALL

C  W - RESTRUCTURE WORST POSSIBLE ELEMENT ONLY

      ELSE IF ((SCHSTR(J:J) .EQ. 'W') .OR.
     &   (SCHSTR(J:J) .EQ. 'w')) THEN
         LIMIT = 1
         CALL RESTA (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &      KKKOLD, NAVAIL, IAVAIL, NNN, LIMIT, IREST, TILT, ERR,
     &      NOROOM)
         IF (NOROOM) THEN
            GO TO 140
         ELSE IF (ERR) THEN
            CALL MESAGE ('ERROR DURING WORST ELEMENT RESTRUCTURE')
            GO TO 140
         END IF
         DOTILT = .FALSE.
         IF (IREST .GE. 1) THEN
            DOSMOO = .TRUE.
            RECT = .FALSE.
            ACTIVE = .TRUE.
            CALL MESAGE ('WORST ELEMENT RESTRUCTURED')
         END IF

C  X - PAVING REGION - QUAD MESH GENERATION

      ELSE IF ((SCHSTR(J:J) .EQ. 'X') .OR.
     &   (SCHSTR(J:J) .EQ. 'x')) THEN
         CONTINUE

C  Y - CONTROL UNDER- OR OVER-RELAXATION

      ELSE IF ((SCHSTR(J:J) .EQ. 'Y') .OR.
     &   (SCHSTR(J:J) .EQ. 'y')) THEN
         WFDEL = .1
         IF (IISIGN .GE. 0) THEN
            WF = WF + WFDEL
         ELSE
            WF = MAX(WF - WFDEL, WFDEL)
         END IF
         DOSMOO = .TRUE.
         WRITE(*, 10040) WF

C  Z - PROCESS REGION WITH HOLES

      ELSE IF ((SCHSTR(J:J) .EQ. 'Z') .OR.
     &   (SCHSTR(J:J) .EQ. 'z')) THEN
         IF (NHPR(L) .EQ. 0) THEN
            CALL MESAGE ('NO HOLES DEFINED IN THIS REGION')
            GO TO 130
         ELSE IF (JJHOLE .EQ. 0) THEN
            JJHOLE = IFHOLE(L)
            MXHOLE = JJHOLE + NHPR(L) - 1
         ELSE
            JJHOLE = JJHOLE + 1
         END IF
         IF (JJHOLE .GT. MXHOLE) THEN
            CALL MESAGE ('ALL HOLES PROCESSED FOR REGION')
            GO TO 130
         END IF
         ADDLNK = .FALSE.
         CALL LTSORT (MR, LINKR, IHLIST(JJHOLE), JHOLE, ADDLNK)

C  JHOLE IS NEGATIVE FOR REGIONS ON BODY CARD WITH LESS THAN THREE INTERVALS

         JHOLE = ABS(JHOLE)
         CALL LTSORT (MR, LINKSC, ABS(IREGN(JHOLE)), IPNTR,
     &      ADDLNK)
         IF ((ABS(IREGN(JHOLE)) .LE. N(24)) .AND.
     &      (IPNTR .GT. 0)) THEN
            SCHOLE = SCHEME(IPNTR)
         ELSE
            SCHOLE = ' '
         END IF
         IF (STEP) THEN
            CALL STRCUT (SCHOLE)
            CALL STRLNG (SCHOLE, LENHOL)
            WRITE(*, 10000) SCHOLE(1:LENHOL)
            CALL INTRUP ('USE CURRENT HOLE SCHEME', IANS, MCOM, ICOM,
     &         JCOM, CIN, IIN, RIN, KIN)
            IF (.NOT.IANS) THEN
  110          CONTINUE
               IF (ICOM .LE. JCOM) THEN
                  SCHOLE = CIN(ICOM)
                  ICOM = ICOM + 1
                  IANS = .TRUE.
               ELSE
                  CALL INQSTR ('ENTER HOLE PROCESSING SCHEME: ', SCHOLE)
               END IF
               IF ((SCHOLE(1:1) .EQ. 'H') .OR.
     &            (SCHOLE(1:1) .EQ. 'h')) THEN
                  CALL MESAGE (' ')
                  CALL HELP_FQ (13)
                  CALL MESAGE (' ')
                  GO TO 110
               END IF
            END IF
         END IF

         CALL STRCUT (SCHOLE)
         CALL STRLNG (SCHOLE, LENHOL)
         INSIDE = 0
         DO 120 M = 1, LENHOL
            IF ((SCHOLE(M:M) .EQ. 'L') .OR.
     &         (SCHOLE(M:M) .EQ. 'l')) INSIDE = INSIDE + 1
  120    CONTINUE

         NPRM = NPRM + 1

         CALL ZHOLE (MP, ML, MS, MR, NSPR(JHOLE), MXNL, MXNPER, MAXPRM,
     &      NPRM, MAXNBC, MAXSBC, KNBC, KSBC, IREGN(JHOLE), IPOINT,
     &      COOR, IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN,
     &      ISBOUN, ISIDE, NLPS, IFLINE, ILLIST, ISLIST,
     &      IFSIDE(JHOLE), NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF,
     &      IFSB, LISTSB, LINKP, LINKL, LINKS, LINKPB, LINKLB, LINKSB,
     &      X, Y, NID, LISTL, MARKED, NL, LSTNBC, MXND, XN, YN, NUID,
     &      LXK, KXL, NXL, LXN, NXH, NPERIM, NNN, NNNOLD, KKK, LLL,
     &      IAVAIL, NAVAIL, JHOLE, INSIDE, EPS, NOROOM, ERR, AMESUR,
     &      XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR,
     &      MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
         IF (NOROOM) THEN
            GO TO 140
         ELSE IF (ERR) THEN
            CALL MESAGE ('HOLE PROCESSING FAILED')
            GO TO 140
         ELSE
            DOSMOO = .TRUE.
            RECT = .FALSE.
            ACTIVE = .TRUE.
            IF (INSIDE .GT. 0) THEN
               NIT3 = 3
               RONE = 1.0
               CALL SMOGS (MXND, XN, YN, NXL, LXN, NNN, NNNOLD, NIT3,
     &            EPS, RONE)
            END IF
            IF (IPNTR .LE. 0) THEN
               ADDLNK = .TRUE.
               CALL LTSORT (MR, LINKSC, ABS(IREGN(JHOLE)), IPNTR,
     &            ADDLNK)
               ADDLNK = .FALSE.
            END IF
            if (ipntr .gt. 0) SCHEME(IPNTR) = SCHOLE
            CALL MESAGE ('HOLE PROCESSING COMPLETED')
         END IF
  130    CONTINUE

C  ( - START LOOP

      ELSE IF (SCHSTR(J:J) .EQ. '(') THEN
         IF (NLEFTP .GE. 10) THEN
            CALL MESAGE ('TOO MANY NESTED LOOPS IN THE SCHEME')
            GO TO 140
         END IF
         NLEFTP = NLEFTP + 1
         LACT(NLEFTP) = ACTIVE
         ILPC(NLEFTP) = J
         ACTIVE = .FALSE.

C  ) - END OF LOOP - CHECK FOR ACTIVITY

      ELSE IF (SCHSTR(J:J) .EQ. ')') THEN
         IF (NLEFTP .LE. 0) THEN
            CALL MESAGE ('THERE IS NO LEFT PARENTHESIS TO')
            CALL MESAGE ('MATCH THE RIGHT PARENTHESIS')
            CALL MESAGE ('")" IS THUS IGNORED')
         ELSE

C  LOOP BACK

            IF (ACTIVE) THEN
               ACTIVE = .FALSE.
               J = ILPC(NLEFTP)
               LACT(NLEFTP) = .TRUE.

C  LOOP IS COMPLETED

            ELSE
               ACTIVE = LACT(NLEFTP)
               NLEFTP = NLEFTP - 1
            END IF
         END IF

C  + SIGN

      ELSE IF (SCHSTR(J:J) .EQ. '+') THEN
         NEWSGN = +1

C  - SIGN

      ELSE IF (SCHSTR(J:J) .EQ. '-') THEN
         NEWSGN = -1

C  1, 2, ..., 6  SMOOTHING TYPE DECLARATION

      ELSE IF (SCHSTR(J:J) .EQ. '1') THEN
         IF (ISTYPE .NE. 1) THEN
            ISTYPE = 1
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "EQUIPOTENTIAL"')
      ELSE IF (SCHSTR(J:J) .EQ. '2') THEN
         IF (ISTYPE .NE. 2) THEN
            ISTYPE = 2
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "AREA PULL & LAPLACIAN"')
      ELSE IF (SCHSTR(J:J) .EQ. '3') THEN
         IF (ISTYPE .NE. 3) THEN
            ISTYPE = 3
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "CENTROID INVERSE AREA '//
     &      'PUSH & LAPLACIAN"')
      ELSE IF (SCHSTR(J:J) .EQ. '4') THEN
         IF (ISTYPE .NE. 4) THEN
            ISTYPE = 4
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "CENTROID AREA PULL"')
      ELSE IF (SCHSTR(J:J) .EQ. '5') THEN
         IF (ISTYPE .NE. 5) THEN
            ISTYPE = 5
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "LAPLACIAN"')
      ELSE IF (SCHSTR(J:J) .EQ. '6') THEN
         IF (ISTYPE .NE. 6) THEN
            ISTYPE = 6
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "LENGTH WEIGHTED '//
     &      'LAPLACIAN"')
      ELSE IF (SCHSTR(J:J) .EQ. '7') THEN
         IF (ISTYPE .NE. 7) THEN
            ISTYPE = 7
            DOSMOO = .TRUE.
         END IF
         CALL MESAGE ('SMOOTHING TYPE SET TO "LAPLACIAN-ISOPARAMETRIC"')

C  BLANK SCHEME FLAG

      ELSE IF (SCHSTR(J:J) .EQ. ' ') THEN
         IF (J .NE. 1) CALL MESAGE ('BLANK SCHEME COMMAND IGNORED')

C  ILLEGAL SCHEME FLAG

      ELSE
         WRITE(*, 10080) SCHSTR(J:J)
      END IF

C  GET NEXT SCHEME COMMAND

      J = J + 1
      IF (J .LE. LENSCH) THEN
         GO TO 100
      ELSE IF (STEP) THEN
         CALL MESAGE ('------------------------')
         IF (ICOM .LE. JCOM) THEN
            ADDSTR = CIN(ICOM)
            ICOM = ICOM + 1
         ELSE
            CALL INQSTR ('FURTHER PROCESSING STEPS: ', ADDSTR)
         END IF
         CALL STRCUT (ADDSTR)
         CALL STRLNG (ADDSTR, LENA)
         IF ((LENSCH .EQ. 1) .AND. (SCHSTR(J - 1:J - 1) .EQ. ' '))
     &      LENSCH = LENSCH - 1
         IF (LENSCH + LENA .GT. 72) THEN
            CALL MESAGE ('ERROR - SCHEME TOO LONG')
            GO TO 140
         END IF
         SCHSTR(LENSCH + 1:LENSCH + LENA) = ADDSTR(1:LENA)
         CALL STRCUT (SCHSTR)
         J = LENSCH + 1
         CALL STRLNG (SCHSTR, LENSCH)
         GO TO 100
      END IF
      ICODE = IEXIT

C  END OF THIS REGION

  140 CONTINUE
      RETURN
10000 FORMAT (' ', /, ' INITIAL MESH DEFINED USING THIS HOLE SCHEME:' /,
     &   '    ', A)
10010 FORMAT (' ALPHA SMOOTHING PARAMETER FOR EQUAL AREAS SET TO:',
     &   F6.3)
10020 FORMAT (' SCHEME: ', A, /,
     &   ' FOR REGION:', I5, /,
     &   ' CANNOT BE SAVED HERE DUE TO DIMENSIONING CONSTRAINTS')
10030 FORMAT (' RO SMOOTHING PARAMETER FOR RELAXATION SET TO:', F6.3)
10040 FORMAT (' WF SMOOTHING PARAMETER FOR ISOPARAMETIC SET TO:', F6.3)
10050 FORMAT (' NO OF SMOOTHING ITERATIONS SET TO:', I5)
10060 FORMAT (' SMOOTHING TOLERANCE SET TO:', G14.7, /,
     &   ' SMOOTHING EPSILON SET TO:', G14.7)
10070 FORMAT (' SMALLEST ANGLE OF ELEMENT TO BE DELETED SET TO:', F6.3)
10080 FORMAT (' ILLEGAL SCHEME COMMAND: ', A1)
      END
