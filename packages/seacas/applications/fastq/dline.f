C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DLINE (MP, ML, COOR, LINKP, KNUM, KT, IP1, IP2, IP3,
     &   NUMPLT, X1, Y1, TEST, GETMAX, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************

C  SUBROUTINE DLINE = DRAWS A LINE ACCORDING TO THE CURRENT DEFINITION
C                     OR SIMPLY GETS THE MAX/MIN FOR LINES GETMAX=.TRUE.

C***********************************************************************

C  VARIABLES USED:
C     IP1   = POINTER FOR THE FIRST POINT
C     IP2   = POINTER FOR THE SECOND POINT
C     IP3   = POINTER FOR THE THIRD POINT

C***********************************************************************

      DIMENSION COOR (2, MP), LINKP (2, MP)

      CHARACTER*72 DUMMY

      LOGICAL NUMPLT, ADDLNK, TEST, GETMAX, ERR

      PI = ATAN2(0.0, -1.0)

      IF (TEST)WRITE (12, 10000)'SP', KNUM, ';'
      ADDLNK = .FALSE.

C  DEFINE FIRST POINT EXACTLY AND BRANCH

      CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
      CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
      IF ((IPNTR1 .LE. 0).OR. (IPNTR2 .LE. 0))GOTO 140
      IF (IP3 .NE. 0) THEN
         CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
         IF (IPNTR3 .LE. 0)GOTO 140
      ELSE
         IPNTR3 = 0
      ENDIF
      X1 = COOR (1, IPNTR1)
      Y1 = COOR (2, IPNTR1)

C  STRAIGHT LINE GENERATION

      IF (KT .EQ. 1) THEN
         X2 = COOR (1, IPNTR2)
         Y2 = COOR (2, IPNTR2)
         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, X2, XMAX)
            YMAX = AMAX1 (Y1, Y2, YMAX)
            XMIN = AMIN1 (X1, X2, XMIN)
            YMIN = AMIN1 (Y1, Y2, YMIN)
            GOTO 140
         ENDIF
         XMID =  (X1 + X2) * .5
         YMID =  (Y1 + Y2) * .5

C  CORNER GENERATION

      ELSEIF (KT .EQ. 2) THEN
         X2 = COOR (1, IPNTR3)
         Y2 = COOR (2, IPNTR3)
         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, X2, COOR (1, IPNTR2), XMAX)
            YMAX = AMAX1 (Y1, Y2, COOR (2, IPNTR2), YMAX)
            XMIN = AMIN1 (X1, X2, COOR (1, IPNTR2), XMIN)
            YMIN = AMIN1 (Y1, Y2, COOR (2, IPNTR2), YMIN)
            GOTO 140
         ENDIF
         IF (TEST)WRITE (12, 10010)
     &      'PU;PA', INT (X1 * 1000.), ', ',
     &      INT (Y1 * 1000.), ';PD;PA',
     &      INT (X2 * 1000.), ', ', INT (Y2 * 1000.), ';'
         CALL MPD2VC (1, X1, Y1, X2, Y2)
         XMID = X1 + ((X2 - X1) * .25)
         YMID = Y1 + ((Y2 - Y1) * .25)
         X1 = X2
         Y1 = Y2
         X2 = COOR (1, IPNTR2)
         Y2 = COOR (2, IPNTR2)

C  CIRCULAR ARC GENERATION

      ELSEIF ((KT .EQ. 3).OR. (KT .EQ. 4).OR. (KT .EQ. 6)) THEN
         CALL ARCPAR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2,
     &      IPNTR3, IP3, XCEN, YCEN, THETA1, THETA2, TANG, R1, R2, ERR,
     &      IDUM1, IDUM2, XK, XA)
         IF (ERR) GOTO 140

C  GENERATE THE CIRCLE

         ANG = THETA1
         DARC = .10
         INC = INT (ABS (TANG) / DARC) + 1
         IF (INC .LE. 6)INC = 6
         DEL = TANG * (1.0 / DBLE(INC))
         IEND = INC - 1
         XK =  (LOG (R2 / R1)) / (THETA2 - THETA1)
         XA = R2 / EXP (XK * THETA2)
         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, XMAX)
            YMAX = AMAX1 (Y1, YMAX)
            XMIN = AMIN1 (X1, XMIN)
            YMIN = AMIN1 (Y1, YMIN)
         ENDIF
         DO 100 I = 1, IEND
            ANG = ANG + DEL
            RADIUS = XA * EXP (XK * ANG)
            X2 = XCEN + COS (ANG) * RADIUS
            Y2 = YCEN + SIN (ANG) * RADIUS
            IF (GETMAX) THEN
               XMAX = AMAX1 (X2, XMAX)
               YMAX = AMAX1 (Y2, YMAX)
               XMIN = AMIN1 (X2, XMIN)
               YMIN = AMIN1 (Y2, YMIN)
            ELSE
               CALL MPD2VC (1, X1, Y1, X2, Y2)
            ENDIF
            IF (TEST)WRITE (12, 10010)
     &         'PU;PA', INT (X1 * 1000.), ', ',
     &         INT (Y1 * 1000.), ';PD;PA',
     &         INT (X2 * 1000.), ', ', INT (Y2 * 1000.), ';'
            X1 = X2
            Y1 = Y2
            IF (I .EQ. INC / 2) THEN
               XMID = X1
               YMID = Y1
            ENDIF
  100    CONTINUE

C  ELIPSE GENERATION

      ELSEIF  (KT .EQ. 7) THEN
         CALL ELPSPR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2, IPNTR3,
     &      IP3, XCEN, YCEN, THETA1, THETA2, TANG, IDUM1, IDUM2, AVALUE,
     &      BVALUE, ERR)
         IF (ERR) GOTO 140

C  GENERATE THE ELIPSE

         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, XMAX)
            YMAX = AMAX1 (Y1, YMAX)
            XMIN = AMIN1 (X1, XMIN)
            YMIN = AMIN1 (Y1, YMIN)
         ENDIF
         DARC = .10
         INC = MAX0 (INT (ABS (TANG) / DARC) + 1, 15)
         DEL = TANG * (1.0 / DBLE(INC))
         IEND = INC - 1
         ANG  =  THETA1
         DO 110 I  =  1, IEND
            ANG  =  ANG + DEL
            RADIUS  =  SQRT  (  (AVALUE **2 * BVALUE **2) /
     &         (  (BVALUE **2 * COS  (ANG) **2) +
     &         (AVALUE **2 * SIN  (ANG) **2) ) )
            X2 = XCEN + COS (ANG) * RADIUS
            Y2 = YCEN + SIN (ANG) * RADIUS
            IF (GETMAX) THEN
               XMAX = AMAX1 (X2, XMAX)
               YMAX = AMAX1 (Y2, YMAX)
               XMIN = AMIN1 (X2, XMIN)
               YMIN = AMIN1 (Y2, YMIN)
            ELSE
               CALL MPD2VC (1, X1, Y1, X2, Y2)
            ENDIF
            IF (TEST)WRITE (12, 10010)
     &         'PU;PA', INT (X1 * 1000.), ', ',
     &         INT (Y1 * 1000.), ';PD;PA',
     &         INT (X2 * 1000.), ', ', INT (Y2 * 1000.), ';'
            X1 = X2
            Y1 = Y2
            IF (I .EQ. INC / 2) THEN
               XMID = X1
               YMID = Y1
            ENDIF
  110    CONTINUE

C  PARABOLA

      ELSEIF (KT .EQ. 5) THEN
         N = 50
         FAC = 1.
         DFF = .02

C  CHECK LEGITIMACY OF DATA

         XMID =  (COOR (1, IPNTR1) + COOR (1, IPNTR2)) * 0.5
         YMID =  (COOR (2, IPNTR1) + COOR (2, IPNTR2)) * 0.5
         DOT =  (COOR (1, IPNTR2) - COOR (1, IPNTR1)) *
     &      (COOR (1, IPNTR3) - XMID) + (COOR (2, IPNTR2) -
     &      COOR (2, IPNTR1)) * (COOR (2, IPNTR3) - YMID)
         PERP = SQRT ((COOR (1, IPNTR2) - COOR (1, IPNTR1)) **2 +
     &      (COOR (2, IPNTR2) - COOR (2, IPNTR1)) **2) *
     &      SQRT ((COOR (1, IPNTR3) - XMID) **2 +
     &      (COOR (2, IPNTR3) - YMID) **2)
         IF (DOT .GE. 0.05 * PERP) THEN
            CALL PLTFLU
            WRITE (*, 10040)KNUM
            GOTO 140
         ENDIF
         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, XMAX)
            YMAX = AMAX1 (Y1, YMAX)
            XMIN = AMIN1 (X1, XMIN)
            YMIN = AMIN1 (Y1, YMIN)
         ENDIF

C  GET ARC LENGTH

         HALFW = SQRT ((COOR (1, IPNTR2) - COOR (1, IPNTR1)) **2 +
     &      (COOR (2, IPNTR2) - COOR (2, IPNTR1)) **2) * 0.5
         IF (HALFW .EQ. 0.) THEN
            CALL PLTFLU
            WRITE (*, 10020)KNUM
            GOTO 140
         ENDIF
         HEIGHT = SQRT ((XMID - COOR (1, IPNTR3)) **2 + (YMID -
     &      COOR (2, IPNTR3)) **2)
         COEF = HEIGHT / HALFW **2
         TCOEF = 2.0 * COEF

C  PARC IS A STATEMENT FUNCTION

         PLEFT = PARC ( - TCOEF * HALFW, TCOEF)
         ARCTOT = 2.0 * PARC (TCOEF * HALFW, TCOEF)
         ARCDEL = DFF * ARCTOT
         ARCNXT = ARCDEL
         ARCNOW = 0.0
         THETA = ATAN2 (COOR (2, IPNTR2) - COOR (2, IPNTR1),
     &      COOR (1, IPNTR2) - COOR (1, IPNTR1))

C  CORRECT FOR ORIENTATION

         CROSS =  (COOR (1, IPNTR3) - XMID) * (COOR (2, IPNTR2) -
     &      COOR (2, IPNTR1)) - (COOR (2, IPNTR3) - YMID) *
     &      (COOR (1, IPNTR2) - COOR (1, IPNTR1))
         IF (CROSS .LT. 0.0)THETA = THETA + PI
         SINT = SIN (THETA)
         COST = COS (THETA)

C  FIND POINTS APPROXIMATELY BY INTEGRATION

         XL =  - HALFW
         FL = SQRT (1.0 + (TCOEF * XL) **2)
         KOUNT = 1
         DELX = 2.0 * HALFW / 200.0
         DO 120 I = 1, 100
            FM = SQRT (1.0 + (TCOEF * (XL + DELX)) **2)
            XR =  - HALFW + DBLE(I) * 2.0 * DELX
            FR = SQRT (1.0 + (TCOEF * XR) **2)
            ARCOLD = ARCNOW
            ARCNOW = ARCNOW + DELX * (FL + 4.0 * FM + FR) / 3.0
            IF (ARCNOW .GE. ARCNXT) THEN

C  COMPUTE POSITION IN LOCAL COORDINATE SYSTEM

               FRAC =  (ARCNXT - ARCOLD) / (ARCNOW - ARCOLD)
               XK = XL + FRAC * 2.0 * DELX
               YK = COEF * XK **2

C  CORRECT FOR ORIENTATION PROBLEM

               IF (CROSS .LT. 0.0)XK = - XK

C  ROTATE IN LINE WITH GLOBAL COORDINATE SYSTEM

               ROTX = XK * COST - YK * SINT
               ROTY = YK * COST + XK * SINT

C  RESTORE XK

               IF (CROSS .LT. 0.0)XK = - XK

C  TRANSLATE

               KOUNT = KOUNT + 1
               X2 = ROTX + COOR (1, IPNTR3)
               Y2 = ROTY + COOR (2, IPNTR3)
               IF (TEST)WRITE (12, 10010)
     &            'PU;PA', INT (X1 * 1000.), ', ',
     &            INT (Y1 * 1000.), ';PD;PA',
     &            INT (X2 * 1000.), ', ', INT (Y2 * 1000.), ';'
               IF (GETMAX) THEN
                  XMAX = AMAX1 (X2, XMAX)
                  YMAX = AMAX1 (Y2, YMAX)
                  XMIN = AMIN1 (X2, XMIN)
                  YMIN = AMIN1 (Y2, YMIN)
               ELSE
                  CALL MPD2VC (1, X1, Y1, X2, Y2)
               ENDIF
               X1 = X2
               Y1 = Y2

C  PREPARE FOR NEXT POINT

               IF (KOUNT .GE. N - 1)GOTO 130
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
  120    CONTINUE
  130    CONTINUE
         XMID = COOR (1, IPNTR3)
         YMID = COOR (2, IPNTR3)
      ENDIF

C     NORMAL EXIT
C     DEFINE LAST POINT EXACTLY

      X2 = COOR (1, IPNTR2)
      Y2 = COOR (2, IPNTR2)
      IF (GETMAX) THEN
         XMAX = AMAX1 (X2, XMAX)
         YMAX = AMAX1 (Y2, YMAX)
         XMIN = AMIN1 (X2, XMIN)
         YMIN = AMIN1 (Y2, YMIN)
         GOTO 140
      ENDIF
      IF (TEST)WRITE (12, 10010)
     &   'PU;PA', INT (X1 * 1000.), ', ',
     &   INT (Y1 * 1000.), ';PD;PA',
     &   INT (X2 * 1000.), ', ', INT (Y2 * 1000.), ';'
      CALL MPD2VC (1, X1, Y1, X2, Y2)
      CALL PLTFLU

C  PLOT THE LINE NUMBER IF DESIRED

      IF (KNUM .GT. 0) THEN
         CALL MP2PT (1, XMID, YMID, X1, Y1, MASK)
         IF ((MOD (MASK, 2) .NE. 0).AND. (NUMPLT)) THEN
            CALL GETDUM (KNUM, DUMMY, LEN)
            CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
         ENDIF
      ENDIF

  140 CONTINUE

      RETURN

10000 FORMAT (A2, I6, A1)
10010 FORMAT (A5, I10, A1, I10, A6, I10, A1, I10, A1)
10020 FORMAT (' ZERO LINE LENGTH ENCOUNTERED FOR LINE', I5)
10040 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE A PARABOLA')

      END
