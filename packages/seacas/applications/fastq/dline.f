C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C $Id: dline.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: dline.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1991/03/21 15:44:35  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:06:11  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:06:10  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]DLINE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DLINE (MP, ML, COOR, LINKP, KNUM, KT, IP1, IP2, IP3,
     &   NUMPLT, X1, Y1, TEST, GETMAX, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************
C
C  SUBROUTINE DLINE = DRAWS A LINE ACCORDING TO THE CURRENT DEFINITION
C                     OR SIMPLY GETS THE MAX/MIN FOR LINES GETMAX=.TRUE.
C
C***********************************************************************
C
C  VARIABLES USED:
C     IP1   = POINTER FOR THE FIRST POINT
C     IP2   = POINTER FOR THE SECOND POINT
C     IP3   = POINTER FOR THE THIRD POINT
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LINKP (2, MP)
C
      CHARACTER*72 DUMMY
C
      LOGICAL NUMPLT, ADDLNK, TEST, GETMAX, ERR
C
      PI = ATAN2(0.0, -1.0)
C
      IF (TEST)WRITE (12, 10000)'SP', KNUM, ';'
      ADDLNK = .FALSE.
C
C  DEFINE FIRST POINT EXACTLY AND BRANCH
C
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
C
C  STRAIGHT LINE GENERATION
C
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
C
C  CORNER GENERATION
C
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
     &      'PU;PA', IFIX (X1 * 1000.), ', ',
     &      IFIX (Y1 * 1000.), ';PD;PA',
     &      IFIX (X2 * 1000.), ', ', IFIX (Y2 * 1000.), ';'
         CALL MPD2VC (1, X1, Y1, X2, Y2)
         XMID = X1 + ((X2 - X1) * .25)
         YMID = Y1 + ((Y2 - Y1) * .25)
         X1 = X2
         Y1 = Y2
         X2 = COOR (1, IPNTR2)
         Y2 = COOR (2, IPNTR2)
C
C  CIRCULAR ARC GENERATION
C
      ELSEIF ((KT .EQ. 3).OR. (KT .EQ. 4).OR. (KT .EQ. 6)) THEN
         CALL ARCPAR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2,
     &      IPNTR3, IP3, XCEN, YCEN, THETA1, THETA2, TANG, R1, R2, ERR,
     &      IDUM1, IDUM2, XK, XA)
         IF (ERR) GOTO 140
C
C  GENERATE THE CIRCLE
C
         ANG = THETA1
         DARC = .10
         INC = IFIX (ABS (TANG) / DARC) + 1
         IF (INC .LE. 6)INC = 6
         DEL = TANG * (1.0 / FLOAT (INC))
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
     &         'PU;PA', IFIX (X1 * 1000.), ', ',
     &         IFIX (Y1 * 1000.), ';PD;PA',
     &         IFIX (X2 * 1000.), ', ', IFIX (Y2 * 1000.), ';'
            X1 = X2
            Y1 = Y2
            IF (I .EQ. INC / 2) THEN
               XMID = X1
               YMID = Y1
            ENDIF
  100    CONTINUE
C
C  ELIPSE GENERATION
C
      ELSEIF  (KT .EQ. 7) THEN
         CALL ELPSPR (MP, KT, KNUM, COOR, LINKP, IPNTR1, IPNTR2, IPNTR3,
     &      IP3, XCEN, YCEN, THETA1, THETA2, TANG, IDUM1, IDUM2, AVALUE,
     &      BVALUE, ERR)
         IF (ERR) GOTO 140
C
C  GENERATE THE ELIPSE
C
         IF (GETMAX) THEN
            XMAX = AMAX1 (X1, XMAX)
            YMAX = AMAX1 (Y1, YMAX)
            XMIN = AMIN1 (X1, XMIN)
            YMIN = AMIN1 (Y1, YMIN)
         ENDIF
         DARC = .10
         INC = MAX0 (IFIX (ABS (TANG) / DARC) + 1, 15)
         DEL = TANG * (1.0 / FLOAT (INC))
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
     &         'PU;PA', IFIX (X1 * 1000.), ', ',
     &         IFIX (Y1 * 1000.), ';PD;PA',
     &         IFIX (X2 * 1000.), ', ', IFIX (Y2 * 1000.), ';'
            X1 = X2
            Y1 = Y2
            IF (I .EQ. INC / 2) THEN
               XMID = X1
               YMID = Y1
            ENDIF
  110    CONTINUE
C
C  PARABOLA
C
      ELSEIF (KT .EQ. 5) THEN
         N = 50
         FAC = 1.
         DFF = .02
C
C  CHECK LEGITIMACY OF DATA
C
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
C
C  GET ARC LENGTH
C
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
C
C  PARC IS A STATEMENT FUNCTION
C
         PLEFT = PARC ( - TCOEF * HALFW, TCOEF)
         ARCTOT = 2.0 * PARC (TCOEF * HALFW, TCOEF)
         ARCDEL = DFF * ARCTOT
         ARCNXT = ARCDEL
         ARCNOW = 0.0
         THETA = ATAN2 (COOR (2, IPNTR2) - COOR (2, IPNTR1),
     &      COOR (1, IPNTR2) - COOR (1, IPNTR1))
C
C  CORRECT FOR ORIENTATION
C
         CROSS =  (COOR (1, IPNTR3) - XMID) * (COOR (2, IPNTR2) -
     &      COOR (2, IPNTR1)) - (COOR (2, IPNTR3) - YMID) *
     &      (COOR (1, IPNTR2) - COOR (1, IPNTR1))
         IF (CROSS .LT. 0.0)THETA = THETA + PI
         SINT = SIN (THETA)
         COST = COS (THETA)
C
C  FIND POINTS APPROXIMATELY BY INTEGRATION
C
         XL =  - HALFW
         FL = SQRT (1.0 + (TCOEF * XL) **2)
         KOUNT = 1
         DELX = 2.0 * HALFW / 200.0
         DO 120 I = 1, 100
            FM = SQRT (1.0 + (TCOEF * (XL + DELX)) **2)
            XR =  - HALFW + FLOAT (I) * 2.0 * DELX
            FR = SQRT (1.0 + (TCOEF * XR) **2)
            ARCOLD = ARCNOW
            ARCNOW = ARCNOW + DELX * (FL + 4.0 * FM + FR) / 3.0
            IF (ARCNOW .GE. ARCNXT) THEN
C
C  COMPUTE POSITION IN LOCAL COORDINATE SYSTEM
C
               FRAC =  (ARCNXT - ARCOLD) / (ARCNOW - ARCOLD)
               XK = XL + FRAC * 2.0 * DELX
               YK = COEF * XK **2
C
C  CORRECT FOR ORIENTATION PROBLEM
C
               IF (CROSS .LT. 0.0)XK = - XK
C
C  ROTATE IN LINE WITH GLOBAL COORDINATE SYSTEM
C
               ROTX = XK * COST - YK * SINT
               ROTY = YK * COST + XK * SINT
C
C  RESTORE XK
C
               IF (CROSS .LT. 0.0)XK = - XK
C
C  TRANSLATE
C
               KOUNT = KOUNT + 1
               X2 = ROTX + COOR (1, IPNTR3)
               Y2 = ROTY + COOR (2, IPNTR3)
               IF (TEST)WRITE (12, 10010)
     &            'PU;PA', IFIX (X1 * 1000.), ', ',
     &            IFIX (Y1 * 1000.), ';PD;PA',
     &            IFIX (X2 * 1000.), ', ', IFIX (Y2 * 1000.), ';'
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
C
C  PREPARE FOR NEXT POINT
C
               IF (KOUNT .GE. N - 1)GOTO 130
               ARCDEL = ARCDEL * FAC
               ARCNXT = ARCNXT + ARCDEL
C
C  RESTART INTEGRATION
C
               XR = XK
               FR = SQRT (1.0 + (TCOEF * XR) **2)
C
C  CORRECT FOR INTEGRATION ERROR
C
               ARCNOW = PARC (TCOEF * XR, TCOEF) - PLEFT
            ENDIF
            XL = XR
            FL = FR
  120    CONTINUE
  130    CONTINUE
         XMID = COOR (1, IPNTR3)
         YMID = COOR (2, IPNTR3)
      ENDIF
C
C     NORMAL EXIT
C     DEFINE LAST POINT EXACTLY
C
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
     &   'PU;PA', IFIX (X1 * 1000.), ', ',
     &   IFIX (Y1 * 1000.), ';PD;PA',
     &   IFIX (X2 * 1000.), ', ', IFIX (Y2 * 1000.), ';'
      CALL MPD2VC (1, X1, Y1, X2, Y2)
      CALL PLTFLU
C
C  PLOT THE LINE NUMBER IF DESIRED
C
      IF (KNUM .GT. 0) THEN
         CALL MP2PT (1, XMID, YMID, X1, Y1, MASK)
         IF ((MOD (MASK, 2) .NE. 0).AND. (NUMPLT)) THEN
            CALL GETDUM (KNUM, DUMMY, LEN)
            CALL PLTXTH (X1, Y1, DUMMY (1:LEN))
         ENDIF
      ENDIF
C
  140 CONTINUE
C
      RETURN
C
10000 FORMAT (A2, I6, A1)
10010 FORMAT (A5, I10, A1, I10, A6, I10, A1, I10, A1)
10020 FORMAT (' ZERO LINE LENGTH ENCOUNTERED FOR LINE', I5)
10040 FORMAT (' POINTS GIVEN FOR LINE', I5, ' DO NOT DEFINE A PARABOLA')
C
      END
