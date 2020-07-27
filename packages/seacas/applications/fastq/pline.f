C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PLINE (MP, ML, MAXNP, MAXNBC, MAXSBC, IPOINT,
     &   COOR, LINKP, KNUM, KT, NINT, FAC, IP1, IP2, IP3, X, Y, NID,
     &   IPBC1, IPBC2, ILBC, ISBC, LINKPB, NPPF, IFPB, LISTPB, LINKLB,
     &   NLPF, IFLB, LISTLB, LINKSB, NSPF, IFSB, LISTSB, LSTNBC, KNBC,
     &   KSBC, ERR, TEST, REAL, COUNT, NOROOM, AMESUR, XNOLD, YNOLD,
     &   NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD,
     &   NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS,
     &   SIZMIN, EMAX, EMIN, GRAPH, DXMAX)
C***********************************************************************

C  SUBROUTINE PLINE = PRODUCES A NODE STRING FOR THE K'TH LINE IN THE
C                     LINE TABLE

C***********************************************************************

C  VARIABLES USED:
C     NID   = AN ARRAY OF UNIQUE NODE IDENTIFIERS.
C     REAL  = .TRUE. FOR AN ACTUAL GENERATION
C           = .FALSE. FOR A TRIAL GENERATION
C     ERR   = .TRUE. IF AN ERROR WAS ENCOUNTERED
C     IP1   = POINTER FOR THE FIRST POINT
C     IP2   = POINTER FOR THE SECOND POINT
C     IP3   = POINTER FOR THE THIRD POINT
C     MAXNP = MAXIMUM NUMBER OF NODES ON THE PERIMETER
C     NOTE: MAXNP MUST BE ADJUSTED FOR THE CURRENT LOCATION
C           IN X, Y, &NID
C     KT    = THE LINE TYPE:
C           = 1 FOR STRAIGHT LINES
C           = 2 FOR CORNER LINES
C           = 3 FOR ARC WITH CENTER GIVEN
C           = 4 FOR ARC WITH THIRD POINT ON THE ARC
C           = 5 FOR PARABOLA
C           = 6 FOR ARC WITH RADIUS GIVEN

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)
      DIMENSION X (MAXNP), Y (MAXNP), NID (MAXNP)
      DIMENSION LINKPB (2, MP), NPPF (MP), IFPB (MP), LISTPB (2, MP)
      DIMENSION LINKLB (2, ML), NLPF (ML), IFLB (ML), LISTLB (2, ML)
      DIMENSION LINKSB (2, ML), NSPF (ML), IFSB (ML), LISTSB (2, ML)
      DIMENSION LSTNBC (MAXNBC)

      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL ERR, REAL, TEST, ADDLNK, COUNT, NOROOM, REMESH
      LOGICAL GRAPH

      PI = ATAN2(0.0, -1.0)

      EPS = .01
      ERR = .FALSE.
      ADDLNK = .FALSE.
      NOROOM = .FALSE.
      TWOPI  =  PI + PI

C     COMPUTE FRACTION OF TOTAL LENGTH FOR FIRST INTERVAL

      N = IABS (NINT) + 1
      IF (N .LE. 1)THEN
         WRITE (*, 10000)KNUM
         ERR = .TRUE.
         GOTO 340
      ENDIF
      IF (N .GT. MAXNP) THEN
        WRITE (*,*) 'INTERNAL ERROR: Intervals exceed space'
        STOP
      END IF

      DFF = 1.0/DBLE(N - 1)
      IF (ABS (1.0 - FAC) .GT. 1.0E-6)
     &   DFF =  (FAC - 1.0)/ (FAC ** (N - 1) - 1.0)

C  DEFINE FIRST POINT EXACTLY AND BRANCH

      X (1) = COOR (1, IP1)
      Y (1) = COOR (2, IP1)
      IF (N .GT. 2)THEN

C  STRAIGHT LINE GENERATION

         IF (KT .EQ. 1)THEN
            YDIFF = COOR (2, IP2) - COOR (2, IP1)
            XDIFF = COOR (1, IP2) - COOR (1, IP1)
            D = SQRT (YDIFF **2 + XDIFF **2)
            IF (D .EQ. 0.)THEN
               WRITE (*, 10010)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF
            IF (REMESH) THEN

               CALL STRSIZ (MAXNP, X, Y, NINT, N, COOR(1,IP2),
     &            COOR(2,IP2), XDIFF, YDIFF, D, ERR, TEST, XNOLD, YNOLD,
     &            NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD,
     &            NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS,
     &            SIZMIN, EMAX, EMIN, GRAPH, DXMAX)
               IF (ERR) GOTO 340
            ELSE
               DEL = D * DFF
               DO 100 I = 2, N - 1
                  PART = DEL/D
                  X (I) = X (I - 1) + PART * XDIFF
                  Y (I) = Y (I - 1) + PART * YDIFF
                  DEL = DEL * FAC
  100          CONTINUE
            ENDIF

C  CORNER GENERATION

         ELSEIF (KT .EQ. 2)THEN
            XDA = COOR (1, IP3) - COOR (1, IP1)
            YDA = COOR (2, IP3) - COOR (2, IP1)
            XDB = COOR (1, IP2) - COOR (1, IP3)
            YDB = COOR (2, IP2) - COOR (2, IP3)
            DA = SQRT (XDA **2 + YDA **2)
            DB = SQRT (XDB **2 + YDB **2)
            IF ((DA .EQ. 0.).OR. (DB .EQ. 0.))THEN
               WRITE (*, 10010)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF
            D = DA + DB
            DEL = D * DFF
            SUM = 0.0

C  BREAK N INTO TWO PARTS APPROPRIATELY

            DO 110 I = 2, N - 1
               KI = I
               SUM = SUM + DEL
               IF (SUM + 0.5 * DEL .GT. DA)GOTO 120
               DEL = DEL * FAC
  110       CONTINUE

C  GENERATE FIRST SIDE OF CORNER

  120       CONTINUE
            NA = KI
            DFF = 1.0/DBLE(NA - 1)
            IF (ABS (1.0 - FAC) .GT. 1.0E-6)
     &         DFF =  (FAC - 1.0)/ (FAC ** (NA - 1) - 1.0)
            DEL = DA * DFF
            DO 130 I = 2, NA
               PART = DEL/DA
               X (I) = X (I - 1) + PART * XDA
               Y (I) = Y (I - 1) + PART * YDA
               DEL = DEL * FAC
  130       CONTINUE

C  GENERATE SECOND SIDE OF CORNER

            NB = N - KI + 1
            DFF = 1.0/DBLE(NB - 1)
            IF (ABS (1.0 - FAC) .GT. 1.0E-6)
     &         DFF =  (FAC - 1.0)/ (FAC ** (NB - 1) - 1.0)
            DEL = DB * DFF
            NAP = NA + 1
            DO 140 I = NAP, N - 1
               PART = DEL/DB
               X (I) = X (I - 1) + PART * XDB
               Y (I) = Y (I - 1) + PART * YDB
               DEL = DEL * FAC
  140       CONTINUE

C  CIRCULAR ARC GENERATION

         ELSEIF ((KT .EQ. 3).OR. (KT .EQ. 4).OR. (KT .EQ. 6))THEN

C  ARC WITH CENTER GIVEN
C  ARC GOES FROM 1ST POINT TO 2ND IN *COUNTER-CLOCKWISE* DIRECTION.

            IF (KT .EQ. 3)THEN
               XCEN = COOR (1, IABS (IP3))
               YCEN = COOR (2, IABS (IP3))

C  CIRCLE WITH THIRD POINT ON ARC.

            ELSEIF (KT .EQ. 4)THEN
               THETA1 = ATAN2 (COOR (2, IP3) - COOR (2, IP1),
     &            COOR (1, IP3) - COOR (1, IP1)) + PI/2.0
               THETA2 = ATAN2 (COOR (2, IP3) - COOR (2, IP2),
     &            COOR (1, IP3) - COOR (1, IP2)) + PI/2.0
               DET =  - COS (THETA1) * SIN (THETA2) + COS (THETA2) *
     &            SIN (THETA1)
               X1 = 0.5 * (COOR (1, IP1) + COOR (1, IP3))
               Y1 = 0.5 * (COOR (2, IP1) + COOR (2, IP3))
               X2 = 0.5 * (COOR (1, IP2) + COOR (1, IP3))
               Y2 = 0.5 * (COOR (2, IP2) + COOR (2, IP3))
               R =  ( - SIN (THETA2) * (X2 - X1) + COS (THETA2) *
     &            (Y2 - Y1))/DET
               XCEN = X1 + R * COS (THETA1)
               YCEN = Y1 + R * SIN (THETA1)

C     CIRCLE WITH RADIUS GIVEN

            ELSEIF (KT .EQ. 6)THEN
               DX = 0.5 * (COOR (1, IP2) - COOR (1, IP1))
               DY = 0.5 * (COOR (2, IP2) - COOR (2, IP1))
               CHORD = SQRT (DX * DX + DY * DY)
               R = ABS (COOR (1, IABS (IP3)))
               IF (R .LE. CHORD)THEN
                  XCEN = 0.5 * (COOR (1, IP1) + COOR (1, IP2))
                  YCEN = 0.5 * (COOR (2, IP1) + COOR (2, IP2))
               ELSE
                  ARM = SQRT (R * R - CHORD * CHORD)
                  IF (IP3 .LT. 0)THEN
                     XCEN = COOR (1, IP1) + DX + ARM * DY/CHORD
                     YCEN = COOR (2, IP1) + DY - ARM * DX/CHORD
                  ELSE
                     XCEN = COOR (1, IP1) + DX - ARM * DY/CHORD
                     YCEN = COOR (2, IP1) + DY + ARM * DX/CHORD
                  ENDIF
               ENDIF
            ENDIF
            R1 = SQRT ((COOR (1, IP1) - XCEN) **2 + (COOR (2, IP1) -
     &         YCEN) **2)
            R2 = SQRT ((COOR (1, IP2) - XCEN) **2 + (COOR (2, IP2) -
     &         YCEN) **2)
            IF ((R1 .EQ. 0.).OR. (R2 .EQ. 0.))THEN
               WRITE (*, 10020)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF
            THETA1 = ATAN2 (COOR (2, IP1) - YCEN, COOR (1, IP1) - XCEN)
            THETA2 = ATAN2 (COOR (2, IP2) - YCEN, COOR (1, IP2) - XCEN)

C  ARC WITH THE CENTER GIVEN

            IF (KT .EQ. 3)THEN
               IF ((IP3 .GE. 0).AND. (THETA2 .LE. THETA1))
     &            THETA2 = THETA2 + TWOPI
               IF ((IP3 .LT. 0).AND. (THETA1 .LE. THETA2))
     &            THETA1 = THETA1 + TWOPI
               TANG = THETA2 - THETA1

C  CIRCULAR ARC WITH 3RD POINT ON ARC - CLOCKWISE OR COUNTER-CLOCKWISE

            ELSEIF (KT .EQ. 4)THEN
               THETA3 = ATAN2 (COOR (2, IP3) - YCEN, COOR (1, IP3) -
     &            XCEN)
               IF (THETA2 .LE. THETA1)THETA2 = THETA2 + TWOPI
               IF (THETA3 .LE. THETA1)THETA3 = THETA3 + TWOPI
               TANG = THETA2 - THETA1
               IF (THETA3 .GT. THETA2)TANG = - (TWOPI - TANG)

C     CIRRCULAR ARC WITH RADIUS GIVEN - CLOCKWISE OR COUNTER-CLOCKWISE

            ELSEIF (KT .EQ. 6)THEN
               IF ((IP3 .GE. 0).AND. (THETA2 .LE. THETA1))
     &            THETA2 = THETA2 + TWOPI
               IF ((IP3 .LT. 0).AND. (THETA1 .LE. THETA2))
     &            THETA1 = THETA1 + TWOPI
               TANG = THETA2 - THETA1
            ENDIF

C  GENERATE THE CIRCLE

            ANG = THETA1
            DEL = TANG * DFF
            AA =  (LOG (R2/R1))/ (THETA2 - THETA1)
            BB = R2/EXP (AA * THETA2)
            IF (REMESH) THEN

               CALL CRCSIZ (MAXNP, X, Y, NINT, N, COOR(1,IP2),
     &            COOR(2,IP2), XCEN, YCEN, THETA1, THETA2, TANG, AA, BB,
     &            ERR, TEST, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG,
     &            BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN,
     &            REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &            GRAPH, DXMAX)
               IF (ERR) GOTO 340
            ELSE
               DO 150 I = 2, N - 1
                  ANG = ANG + DEL
                  RADIUS = BB * EXP (AA * ANG)
                  X (I) = XCEN + COS (ANG) * RADIUS
                  Y (I) = YCEN + SIN (ANG) * RADIUS
                  DEL = DEL * FAC
  150          CONTINUE
            ENDIF

C     ELIPSE

         ELSEIF (KT .EQ. 7)THEN

C  GET THE ELIPSE PARAMETERS

            CALL ELPSPR (MP, KT, KNUM, COOR, LINKP, IP1, IP2,
     &         IABS(IP3), IP3, XCEN, YCEN, THETA1, THETA2, TANG, ICCW,
     &         ICW, AVALUE, BVALUE, ERR)
            IF (ERR) GOTO 340

C  GENERATE THE ELIPSE

            ANG = THETA1
            DEL = TANG * DFF
            DO 160 I = 2, N - 1
               ANG = ANG + DEL
               RADIUS = SQRT ( (AVALUE **2 * BVALUE **2) /
     &            ( (BVALUE **2 * COS (ANG) **2) +
     &            (AVALUE **2 * SIN (ANG) **2) ) )
               X (I) = XCEN + COS (ANG) * RADIUS
               Y (I) = YCEN + SIN (ANG) * RADIUS
               DEL = DEL * FAC
  160       CONTINUE

C     PARABOLA

         ELSEIF (KT .EQ. 5)THEN
            IF (N .GT. 250)THEN
               WRITE (*, 10030)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF

C  CHECK LEGITIMACY OF DATA

            XMID =  (COOR (1, IP1) + COOR (1, IP2)) * 0.5
            YMID =  (COOR (2, IP1) + COOR (2, IP2)) * 0.5
            DOT =  (COOR (1, IP2) - COOR (1, IP1)) * (COOR (1, IP3) -
     &         XMID) + (COOR (2, IP2) - COOR (2, IP1)) *
     &         (COOR (2, IP3) - YMID)
            PERP = SQRT ((COOR (1, IP2) - COOR (1, IP1)) **2 +
     &         (COOR (2, IP2) - COOR (2, IP1)) **2) *
     &         SQRT ((COOR (1, IP3) - XMID) **2 + (COOR (2, IP3)
     &         - YMID) **2)
            IF (DOT .GE. 0.05 * PERP)THEN
               WRITE (*, 10040)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF

C  GET ARC LENGTH

            HALFW = SQRT ((COOR (1, IP2) - COOR (1, IP1)) **2 +
     &         (COOR (2, IP2) - COOR (2, IP1)) **2) * 0.5
            IF (HALFW .EQ. 0.)THEN
               WRITE (*, 10010)KNUM
               ERR = .TRUE.
               GOTO 340
            ENDIF
            HEIGHT = SQRT ((XMID - COOR (1, IP3)) **2 +
     &         (YMID - COOR (2, IP3)) **2)
            COEF = HEIGHT/HALFW **2
            TCOEF = 2.0 * COEF

C  PARC IS A STATEMENT FUNCTION

            PLEFT = PARC ( - TCOEF * HALFW, TCOEF)
            ARCTOT = 2.0 * PARC (TCOEF * HALFW, TCOEF)
            ARCDEL = DFF * ARCTOT
            ARCNXT = ARCDEL
            ARCNOW = 0.0
            THETA = ATAN2 (COOR (2, IP2) - COOR (2, IP1),
     &         COOR (1, IP2) - COOR (1, IP1))

C  CORRECT FOR ORIENTATION

            CROSS =  (COOR (1, IP3) - XMID) * (COOR (2, IP2) -
     &         COOR (2, IP1)) - (COOR (2, IP3) - YMID) *
     &         (COOR (1, IP2) - COOR (1, IP1))
            IF (CROSS .LT. 0.0)THETA = THETA + PI
            SINT = SIN (THETA)
            COST = COS (THETA)

C  FIND POINTS APPROXIMATELY BY INTEGRATION

            XL =  - HALFW
            FL = SQRT (1.0 + (TCOEF * XL) **2)
            KOUNT = 1
            DELX = 2.0 * HALFW/200.0
            DO 170 I = 1, 100
               FM = SQRT (1.0 + (TCOEF * (XL + DELX)) **2)
               XR =  - HALFW + DBLE(I) * 2.0 * DELX
               FR = SQRT (1.0 + (TCOEF * XR) **2)
               ARCOLD = ARCNOW
               ARCNOW = ARCNOW + DELX * (FL + 4.0 * FM + FR)/3.0
               IF (ARCNOW .GE. ARCNXT)THEN

C  COMPUTE POSITION IN LOCAL COORDINATE SYSTEM

                  FRAC =  (ARCNXT - ARCOLD)/ (ARCNOW - ARCOLD)
                  XK = XL + FRAC * 2.0 * DELX
                  YK = COEF * XK **2

C  CORRECT FOR ORIENTATION PROBLEM

                  IF (CROSS .LT. 0.0)XK =  - XK

C  ROTATE IN LINE WITH GLOBAL COORDINATE SYSTEM

                  ROTX = XK * COST - YK * SINT
                  ROTY = YK * COST + XK * SINT

C  RESTORE XK

                  IF (CROSS .LT. 0.0)XK =  - XK

C  TRANSLATE

                  KOUNT = KOUNT + 1
                  X (KOUNT) = ROTX + COOR (1, IP3)
                  Y (KOUNT) = ROTY + COOR (2, IP3)

C  PREPARE FOR NEXT POINT

                  IF (KOUNT .GE. N - 1)GOTO 180
                  ARCDEL = ARCDEL * FAC
                  ARCNXT = ARCNXT + ARCDEL

C  RESTART INTEGRATION

                  XR = XK
                  FR = SQRT (1.0 + (TCOEF * XR) **2)

C  CORRECT FOR INTEGRATION ERROR

                  ARCNOW = PARC (TCOEF * XR, TCOEF) - PLEFT
               ENDIF
               XL = XR
               FL = FR
  170       CONTINUE
  180       CONTINUE
         ENDIF

C  SEE IF THE 1 INTERVAL LINE HAS SOME LENGTH TO IT

      ELSE
         IF ((COOR (1, IP1) .EQ. COOR (1, IP2)) .AND. (COOR (2, IP1)
     &      .EQ. COOR (2, IP2)))THEN
            WRITE (*, 10010)KNUM
            ERR = .TRUE.
            GOTO 340
         ENDIF
      ENDIF

C     NORMAL EXIT
C     DEFINE LAST POINT EXACTLY

      X (N) = COOR (1, IP2)
      Y (N) = COOR (2, IP2)
      IF (TEST) GOTO 340

C  DEFINE UNIQUE  (IN THE WHOLE BODY) NODE NUMBERS

      LPART = 1000000000 + KNUM * 100000
      DO 190 I = 1, N
         NID (I) = LPART + I
  190 CONTINUE
      NID (1) = IPOINT (IP1)
      NID (N) = IPOINT (IP2)

C  FLAG PREVIOUSLY USED POINTS WITH NEGATIVES

      IF (NINT .LT. 0)THEN
         DO 200 I = 1, N
            NID (I) =  - NID (I)
  200    CONTINUE
      ENDIF
      IF (IPOINT (IP1) .LT. 0)NID (1) =  - IABS (NID (1))
      IF (IPOINT (IP2) .LT. 0)NID (N) =  - IABS (NID (N))
      IF (IP1 .EQ. IP2)NID (N) =  - IABS (NID (N))

C  WRITE OUT NODAL BOUNDARY CONDITIONS FOR THE FIRST POINT
C  IF THE POINT HAS NOT BEEN USED BEFORE

      IF (IPOINT (IP1) .GT. 0)THEN
         IF (IPBC1 .GT. 0)THEN
            CALL LTSORT (MP, LINKPB, IPBC1, K, ADDLNK)
            IFLAG = IPBC1
  210       CONTINUE
            DO 220 I = IFPB (K), NPPF (K) + IFPB (K) - 1
               CALL LTSORT (MP, LINKP, LISTPB (1, I), IPNTR1, ADDLNK)
               IF (IPNTR1 .EQ. IP1)THEN
                  KNBC = KNBC + 2
                  IF (KNBC .LE. MAXNBC) THEN
                     IF (REAL) THEN
                        LSTNBC (KNBC - 1) = - IFLAG
                        LSTNBC (KNBC) = IABS (NID (1))
                     ENDIF
                  ELSE
                     NOROOM = .TRUE.
                  ENDIF
                  IF (LISTPB (2, I) .GT. 0)THEN
                     CALL LTSORT (MP, LINKPB, LISTPB (2, I), K, ADDLNK)
                     IFLAG = LISTPB (2, I)
                     GOTO 210
                  ELSE
                     GOTO 230
                  ENDIF
               ENDIF
  220       CONTINUE
            WRITE (*, 10060)IFLAG
            ERR = .TRUE.
            GOTO 340
         ENDIF
  230    CONTINUE
      ENDIF

C  WRITE OUT NODAL BOUNDARY CONDITIONS FOR THE SECOND POINT
C  IF THE POINT HAS NOT BEEN USED BEFORE

      IF (IPOINT (IP2) .GT. 0)THEN
         IF (IPBC2 .GT. 0)THEN
            CALL LTSORT (MP, LINKPB, IPBC2, K, ADDLNK)
            IFLAG = IPBC2
  240       CONTINUE
            DO 250 I = IFPB (K), NPPF (K) + IFPB (K) - 1
               CALL LTSORT (MP, LINKP, LISTPB (1, I), IPNTR1, ADDLNK)
               IF (IPNTR1 .EQ. IP2)THEN
                  KNBC = KNBC + 2
                  IF (KNBC .LE. MAXNBC) THEN
                     IF (REAL) THEN
                        LSTNBC (KNBC - 1) = - IFLAG
                        LSTNBC (KNBC) = IABS (NID (N))
                     ENDIF
                  ELSE
                     NOROOM = .TRUE.
                  ENDIF
                  IF (LISTPB (2, I) .GT. 0)THEN
                     CALL LTSORT (MP, LINKPB, LISTPB (2, I), K, ADDLNK)
                     IFLAG = LISTPB (2, I)
                     GOTO 240
                  ELSE
                     GOTO 260
                  ENDIF
               ENDIF
  250       CONTINUE
            WRITE (*, 10060)IFLAG
            ERR = .TRUE.
            GOTO 340
         ENDIF
  260    CONTINUE
      ENDIF

C  WRITE OUT NODAL BOUNDARY CONDITIONS FOR THE LINE
C  IF THE LINE HAS NOT BEEN USED BEFORE

      IF (NINT .GT. 0)THEN
         IF (ILBC .GT. 0)THEN
            CALL LTSORT (ML, LINKLB, ILBC, K, ADDLNK)
            IFLAG = ILBC
  270       CONTINUE
            DO 290 I = IFLB (K), NLPF (K) + IFLB (K) - 1
               IF (LISTLB (1, I) .LE. 0)THEN
                  CALL MESAGE ('PROBLEMS WITH SIDES IN FLAG LIST'//
     &               ' IN PLINE')
               ELSE
                  IF (LISTLB (1, I) .EQ. KNUM)THEN
                     KNBC = KNBC + 1
                     IF (KNBC .LE. MAXNBC) THEN
                        IF (REAL) LSTNBC (KNBC) = - IFLAG
                     ELSE
                        NOROOM = .TRUE.
                     ENDIF
                     DO 280 J = 1, N
                        KNBC = KNBC + 1
                        IF (KNBC .LE. MAXNBC) THEN
                           IF (REAL) LSTNBC (KNBC) = IABS (NID (J))
                        ELSE
                           NOROOM = .TRUE.
                        ENDIF
  280                CONTINUE
                     IF (LISTLB (2, I) .GT. 0)THEN
                        CALL LTSORT (ML, LINKLB, LISTLB (2, I), K,
     &                     ADDLNK)
                        IFLAG = LISTLB (2, I)
                        GOTO 270
                     ELSE
                        GOTO 300
                     ENDIF
                  ENDIF
               ENDIF
  290       CONTINUE
            WRITE (*, 10050)IFLAG
            ERR = .TRUE.
            GOTO 340
         ENDIF
  300    CONTINUE
      ENDIF

C  IF COUNT, THEN COUNT THE SIDE BOUNDARY CONDITIONS FOR THE LINE
C  NOTE: IT DOES NOT MATTER IF THE LINE HAS BEEN USED BEFORE

      IF ((COUNT) .AND. (ISBC .GT. 0))THEN
         CALL LTSORT (ML, LINKSB, ISBC, K, ADDLNK)
         IFLAG = ISBC
  310    CONTINUE
         DO 320 I = IFSB (K), NSPF (K) + IFSB (K) - 1
            IF (LISTSB (1, I) .LT. 0)THEN
               CALL MESAGE ('PROBLEMS WITH SIDES IN FLAG LIST IN PLINE')
            ELSE
               IF (LISTSB (1, I) .EQ. KNUM)THEN
                  KSBC = KSBC + ((N - 1) * 3)
                  IF (KSBC .GT. MAXSBC) NOROOM = .TRUE.
                  IF (LISTSB (2, I) .GT. 0)THEN
                     CALL LTSORT (ML, LINKSB, LISTSB (2, I), K, ADDLNK)
                     IFLAG = LISTSB (2, I)
                     GOTO 310
                  ELSE
                     GOTO 330
                  ENDIF
               ENDIF
            ENDIF
  320    CONTINUE
         WRITE (*, 10070)IFLAG
         ERR = .TRUE.
         GOTO 340
      ENDIF
  330 CONTINUE

C  NORMAL COMPLETION

  340 CONTINUE

      RETURN

10000 FORMAT (' ZERO NUMBER OF INTERVALS FOR LINE', I5)
10010 FORMAT (' ZERO LINE LENGTH ENCOUNTERED FOR LINE', I5)
10020 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE AN ARC')
10030 FORMAT (' MORE THAN 250 POINTS FOR PARABOLA NOT ALLOWED - LINE',
     &   I5)
10040 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE A PARABOLA')
10050 FORMAT (' LINE BOUNDARY FLAG', I5, ' IS NOT PROPERLY LINKED')
10060 FORMAT (' POINT BOUNDARY FLAG', I5, ' IS NOT PROPERLY LINKED')
10070 FORMAT (' SIDE BOUNDARY FLAG', I5, ' IS NOT PROPERLY LINKED')

      END
