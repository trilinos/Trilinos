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

C $Id: adjrow.f,v 1.4 1998/07/14 18:18:16 gdsjaar Exp $
C $Log: adjrow.f,v $
C Revision 1.4  1998/07/14 18:18:16  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.3  1998/03/23 05:17:50  gdsjaar
C Fixed data statement ordering
C
C Revision 1.2  1991/03/21 15:44:16  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:03:25  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:03:23  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]ADJROW.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADJROW (MXND, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, ANGLE, BNSIZE, LNODES, NLOOP, IAVAIL, NAVAIL, XMIN, XMAX,
     &   YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN, LLLOLD, NNNOLD,
     &   NODE, NADJ1, NADJ2, NNN2, GRAPH, VIDEO, KREG, DEFSIZ, ADJTED,
     &   NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE ADJROW = ADJUSTS A ROW OF ELEMENTS BETWEEN TWO CORNERS
C
C***********************************************************************
C
      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES
C
      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
C
      LOGICAL ERR, GRAPH, ADJTED, VIDEO, NOROOM
C
      CHARACTER*3 DEV1
C
      DATA TMIN1 /.80/, TMIN2 /.3/, WMIN1 /1.25/, WMIN2 /1.35/
C
      PI = ATAN2(0.0, -1.0)
      CALL GETIME (TIME1)
      ERR = .FALSE.
      EPS = .0523599
C
C  START BY SETTING UP THE LIMITS OF THE SEARCH
C
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
     &         ( (FLOAT(NLOOP - 2) * PI) / FLOAT (NLOOP)))
            IF (TEPS .LE. EPS) EPS = TEPS
            GOTO 110
         ELSEIF (N2 .EQ. NADJ2) THEN
            TEPS = .95 * (PI -
     &         ( (FLOAT(NLOOP - 2) * PI) / FLOAT (NLOOP)))
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
C
  110 CONTINUE
      N1 = LNODES (3, NADJ1)
      ADJTED = .FALSE.
C
  120 CONTINUE
      IF (N1 .EQ. NADJ2) GOTO 150
C
C  CHECK A STRING OF CONCAVE (< PI) INTERIOR ANGLES FOR NEEDING A
C  TUCK INSERTED SOMEWHERE
C
      IF ((ANGLE (N1) .LT. PI - EPS) .AND. (LNODES (8, N1) .GT. 1) .AND.
     &   (LXN (4, N1) .EQ. 0) .AND. (LXN (3, N1) .GT. 0)) THEN
C
C  ADDED UP THE TURNING ANGLE AND THE AVERAGE SIZE REDUCTION
C
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
C
C  NOW SEE IF THIS PORTION OF THE ROW NEEDS ADJUSTED WITH A TUCK(S)
C
         IF (KANG .GE. 1) THEN
            RATIO = RATIO / FLOAT (KANG)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/19/90
CC* MODIFICATION: ADDED THE TMIN2 CRITERIA FOR INSERTION OF A TUCK.
C**               THIS CRITERIA SHOULD HELP ALLEVIATE THE LONG SKINNY
C**               ELEMENT FORMATIONS WHEN TRANSITIONING.
C
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
C
C  MAKE SURE THAT THE TUCK DOES NOT ELIMINATE THE END NODES FOR THE LOOP
C
               IF (N11 .NE. N11OLD) THEN
                  IF (NADJ2 .EQ. N11OLD) NADJ2 = N11
                  IF (NODE .EQ. N11OLD) NODE = N11
               ENDIF
C
               NNNOLD = NNN
               LLLOLD = LLL
               ADJTED = .TRUE.
C
            ENDIF
         ENDIF
         N1 = N11
         GOTO 120
C
C  CHECK A STRING OF CONVEX (> PI) INTERIOR ANGLES FOR NEEDING A
C  WEDGE INSERTED SOMEWHERE
C
      ELSEIF ((ANGLE (N1) .GE. PI + EPS) .AND. (LXN (3, N1) .GT. 0)
     &   .AND. (LXN (4, N1) .EQ. 0)) THEN
C
C  ADD UP THE TURNING ANGLE AND THE AVERAGE SIZE REDUCTION
C
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
C
C  NOW SEE IF THIS PORTION OF THE ROW NEEDS ADJUSTED WITH A WEDGE(S)
C
         IF (KANG .GE. 1) THEN
            RATIO = RATIO / FLOAT (KANG)
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
C
            ENDIF
         ENDIF
         N1 = N11
         GOTO 120
      ELSE
         N1 = LNODES (3, N1)
         GOTO 120
      ENDIF
C
C  NOW SMOOTH, CALCULATE THE NEW ANGLES, AND PLOT IF NEEDED
C
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
C
  160 CONTINUE
C
      CALL GETIME (TIME2)
      TIMEAJ = TIMEAJ + TIME2 - TIME1
      RETURN
C
      END
