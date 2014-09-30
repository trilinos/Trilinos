C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: addrow.f,v 1.2 2004/01/21 05:18:39 gdsjaar Exp $
C $Log: addrow.f,v $
C Revision 1.2  2004/01/21 05:18:39  gdsjaar
C Initialized several variables identified by valgrind.
C
C Revision 1.1.1.1  1990/11/30 11:03:04  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:03:03  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]ADDROW.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADDROW TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE ADDROW (MXND, MXCORN, MXLOOP, MLN, MAXPRM, NUID, XN,
     &   YN, ZN, LXK, KXL, NXL, LXN, ANGLE, BNSIZE, LNODES, NBEGIN,
     &   NEND, NLOOP, NEXTN1, LINKPR, KPERIM, KKKOLD, LLLOLD, NNNOLD,
     &   IAVAIL, NAVAIL, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL,
     &   KKK, NNN, NNN2, NADJ1, NADJ2, LCORN, KLOOP, GRAPH, VIDEO, KREG,
     &   DONE, SIZEIT, NOROOM, ERR, XNOLD, YNOLD, NXKOLD, LINKEG,
     &   LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN,
     &   REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************
C
C  SUBROUTINE ADDROW = ADDS A ROW OF ELEMENTS BETWEEN TWO CORNERS
C
C***********************************************************************
C
      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES
C
      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
      DIMENSION LCORN (MXCORN)
      DIMENSION NLOOP (MXLOOP), NEXTN1 (MXLOOP)
      DIMENSION LINKPR (3, MAXPRM)
      DIMENSION X(5), Y(5)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL ERR, GRAPH, AMBIG, DONE, VIDEO, SIZEIT, NOROOM
C
      CHARACTER*3 DEV1
C
      nlold = 0
      idum1 = 0
      idum2 = 0
      CALL GETIME (TIME1)
      ERR = .FALSE.
      AMBIG = .FALSE.
      NNNOLD = MIN (NNN, NNNOLD)
      NNN2 = NNNOLD
C
C  IN THE LNODES ARRAY,
C  THE CORNER STATUS IS STORED IN LNODES (1, N1):
C      0 = NOT DECIDED
C      1 = ROW END
C      3 = ROW SIDE
C      5 = ROW CORNER
C      7 = ROW REVERSAL
C  THE PRECEDING NODE IN LNODES (2, N1),
C  THE NEXT NODE IN LNODES (3, N1),
C  THE INTERIOR/EXTERIOR STATUS OF NODE IS IN LNODES (4, N1).
C      1 = EXTERIOR OR ON THE BOUNDARY OF THE MESH
C          (NEGATED FOR SMOOTHING)
C      2 = INTERIOR TO THE MESH (NEGATED FOR SMOOTHING)
C  THE NEXT COUNTERCLOCKWISE LINE IS STORED IN LNODES (5, N1).
C  THE ANGLE STATUS OF LNODES IS STORED IN (6, N1),
C      1 = ROW END ONLY
C      2 = ROW END OR SIDE
C      3 = ROW SIDE ONLY
C      4 = ROW SIDE OR ROW CORNER
C      5 = ROW CORNER ONLY
C      6 = ROW CORNER OR REVERSAL
C      7 = ROW REVERSAL ONLY
C  THE NUMBER OF NODES TO THE NEXT CORNER IS STORED IN (7, N1).
C  THE DEPTH OF THE ROW OF THIS NODE IS STORED IN (8, N1)
C
C  START AT THE FIRST BEGINNING (CORNER) NODE
C  AND GO TO THE END (CORNER) NODE
C
      N1 = NBEGIN
      N0 = LNODES (2, N1)
      NADJ1 = N0
      N2 = LNODES (3, N1)
      N3 = LNODES (3, N2)
C
C  SET UP FOR SMOOTHING AROUND THE BEGINNING AND ENDING OF THE ROW
C
      LNODES (4, NBEGIN) = - IABS (LNODES (4, NBEGIN))
      LNODES (4, NEND) = - IABS (LNODES (4, NEND))
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NBEGIN, ERR)
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NEND, ERR)
C
C  FIRST CHECK TO SEE IF THE ROW IS CIRCULAR - BEGINS AND ENDS AT N1
C
      IF (NEND .EQ. NBEGIN) THEN
C
C  A CIRCLUAR ROW THAT ENDS AT A ROW END (TEAR DROP SHAPE) NEEDS STARTED
C
         IF (LNODES (1, N1) .EQ. 1) THEN
            NEND = N0
            NADJ1 = NNN + 1
C
C  A CIRCLUAR ROW THAT HAS NO END, CORNER, OR REVERSAL (CIRCLE SHAPE)
C  NEEDS STARTED
C
         ELSEIF (LNODES (1, N1) .EQ. 3) THEN
            CALL EXTND1 (MXND, XN, YN, ANGLE, N0, N1, N2, X(1), Y(1),
     &         DIST1)
            CALL EXTND1 (MXND, XN, YN, ANGLE, N1, N2, N3, X(2), Y(2),
     &         DIST2)
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD2ND TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
            CALL ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &         BNSIZE, LNODES, X(1), Y(1), X(2), Y(2), DIST1, DIST2,
     &         NNN, LLL, KKK, N1, N2, NLOOP (1), SIZEIT, ERR, NOROOM,
     &         XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK,
     &         NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &         REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            N0 = NNN
            NADJ1 = NNN - 1
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N1, NNN-1)
               CALL D2NODE (MXND, XN, YN, NNN-1, NNN)
               CALL D2NODE (MXND, XN, YN, NNN, N2)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
            N1 = N2
C
C  A CIRCLUAR ROW THAT ENDS AT A ROW CORNER NEEDS STARTED
C
         ELSEIF (LNODES (1, N1) .EQ. 5) THEN
            AMBIG = .TRUE.
            LAMBIG = LNODES (5, N1)
            CALL EXTND3 (MXND, XN, YN, ANGLE (N1), N0, N1, N2, X, Y,
     &         DIST)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADD3ND TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
            CALL ADD3ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &         BNSIZE, LNODES, X, Y, DIST, NNN, LLL, KKK, N1, NLOOP (1),
     &         SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG,
     &         LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &         REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &         EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
C
            NADJ1 = NNN - 2
            LNODES (1, N1) = 1
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N1, NNN-2)
               CALL D2NODE (MXND, XN, YN, NNN-2, NNN-1)
               CALL D2NODE (MXND, XN, YN, NNN-1, NNN)
               CALL D2NODE (MXND, XN, YN, NNN, N1)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
            N0 = NNN
            N3 = LNODES (2, N2)
            GOTO 110
C
C  A CIRCLUAR ROW THAT ENDS AT A ROW REVERSAL NEEDS STARTED
C
         ELSEIF (LNODES (1, N1) .EQ. 7) THEN
            CALL EXTND5 (MXND, XN, YN, ANGLE (N2), N1, N2, N3, X, Y,
     &         DIST)
            CALL ADD3ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &         BNSIZE, LNODES, X, Y, DIST, NNN, LLL, KKK, N1, NLOOP (1),
     &         SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG,
     &         LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &         REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &         EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N1, NNN-2)
               CALL D2NODE (MXND, XN, YN, NNN-2, NNN-1)
               CALL D2NODE (MXND, XN, YN, NNN-1, NNN)
               CALL D2NODE (MXND, XN, YN, NNN, N1)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
C
            N0 = NNN
            CALL ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &         BNSIZE, LNODES, X(4), Y(4), X(5), Y(5), DIST, DIST, NNN,
     &         LLL, KKK, NNN, N1, NLOOP (1), SIZEIT, ERR, NOROOM, XNOLD,
     &         YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &         NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            LNODES (1, N1) = 1
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N0, NNN-1)
               CALL D2NODE (MXND, XN, YN, NNN-1, NNN)
               CALL D2NODE (MXND, XN, YN, NNN, N1)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
C
         ELSE
            CALL MESAGE ('PROBLEMS IN ADDROW - THE CIRCLUAR ROW')
            CALL MESAGE ('HAS NO END POINT CLASSIFICATION')
            ERR = .TRUE.
            GOTO 120
         ENDIF
      ENDIF
C
      KOUNT = 0
  100 CONTINUE
C
      KOUNT = KOUNT + 1
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      N3 = LNODES (3, N2)
  110 CONTINUE
C
C  NOW ADD THE NEXT ELEMENT(S) TO THE ROW
C
      IF (N2 .NE. NEND) THEN
C
C  IF N2 IS A ROW SIDE NODE, THEN PROJECT AT THE MIDLINE ANGLE
C
         IF (LNODES (1, N2) .EQ. 3) THEN
            CALL EXTND1 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO ADDNOD TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
            CALL ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, ANGLE,
     &         BNSIZE, LNODES, X(1), Y(1), DIST, NNN, KKK, LLL,
     &         N0, N1, N2, AMBIG, LAMBIG, SIZEIT, ERR, NOROOM, XNOLD,
     &         YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &         NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            AMBIG = .FALSE.
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N0, NNN)
               CALL D2NODE (MXND, XN, YN, N2, NNN)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
            N1 = N2
            GOTO 100
C
C  IF N2 IS A ROW CORNER NODE,
C  THEN PROJECT TWO NODES (1 & 3) AT AVERAGE ANGLES AND ANOTHER (2)
C  AS AN ISOPARAMETRIC TYPE
C
         ELSEIF (LNODES (1, N2) .EQ. 5) THEN
            AHOLD = ANGLE (N2)
            IF (NADJ1 .EQ. N1) NADJ1 = N0
            CALL EXTND3 (MXND, XN, YN, AHOLD, N1, N2, N3, X, Y, DIST)
C
            CALL ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, ANGLE,
     &         BNSIZE, LNODES, X(1), Y(1), DIST, NNN, KKK, LLL,
     &         N0, N1, N2, AMBIG, LAMBIG, SIZEIT, ERR, NOROOM, XNOLD,
     &         YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &         NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            AMBIG = .FALSE.
            NEW1 = NNN
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N0, NEW1)
               CALL D2NODE (MXND, XN, YN, N2, NEW1)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
C
C  TRY A PINCH TO MAKE SURE THAT CONTINUING DOES NOT OVERLAY
C  AN ALREADY GOOD CLOSURE
C
            CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &         NNNOLD, LLLOLD, ERR, NOROOM)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            LLLOLD = LLL
            NNNOLD = NNN
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP (1), XMIN, XMAX,
     &         YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
            CALL GETIME (TIME1)
C
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
               IF (VIDEO) CALL SNAPIT (2)
            ENDIF
C
            CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         NLOOP (1), ANGLE, LNODES, NADJ1, LLL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
            IF (ERR) GOTO 120
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL PINCH (MXND, MXCORN, MLN, NUID, XN, YN, ZN, LXK, KXL,
     &         NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1, NLOOP (1),
     &         KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &         XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &         LCORN, NCORN, IDUM1, IDUM2, GRAPH, VIDEO, KREG, NOROOM,
     &         ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 120
            IF (DONE) GOTO 120
            CALL GETIME (TIME1)
C
C  TRY A COLAPS TO MAKE SURE THAT CONTINUING DOES NOT COMPLICATE
C  AN ALREADY OVERLAPPED SIDE
C
            IF (NLOOP (1) .GT. 6) THEN
               NLOLD = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               CALL COLAPS (MXND, MXCORN, MLN, MXLOOP, NUID, XN, YN,
     &            ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1,
     &            KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &            XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &            LCORN, NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO,
     &            KREG, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               IF (DONE) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  CHECK TO SEE IF ANY OF THE CONCURRENT PERIMETERS OVERLAP
C
            IF (LINKPR (2, KPERIM) .NE. 0) THEN
               LINKPR (3, KPERIM) = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               IDUM = NADJ1
               CALL PCROSS (MXND, MXCORN, MLN, MXLOOP, MAXPRM,
     &            NUID, XN, YN, ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES,
     &            BNSIZE, LINKPR, KPERIM, IDUM, NADJ1, NEW1, KKKOLD,
     &            LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN, LCORN,
     &            NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO, KREG,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  THE ROW CAN BE CONTINUED IF THE NEW NODE AND THE ENDING NODE
C  IS STILL ON THE BOUNDARY AFTER THE PINCH, COLAPS AND PCROSS
C
            IF ( (NLOLD .EQ. NLOOP(1)) .AND.
     &         (IABS (LNODES (4, NEW1)) .EQ. 1) .AND.
     &         (IABS (LNODES (4, NEND)) .EQ. 1) .AND.
     &         (LNODES (2, N2) .EQ. NEW1) .AND.
     &         (LXN (1, NEW1) .GT. 0) .AND.
     &         (ANGLE (N2) .GT. 2.3561945) ) THEN
               CALL EXTND3 (MXND, XN, YN, AHOLD, N1, N2, N3, X, Y, DIST)
               CALL ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &            BNSIZE, LNODES, X(2), Y(2), X(3), Y(3), DIST, DIST,
     &            NNN, LLL, KKK, NEW1, N2, NLOOP (1), SIZEIT, ERR,
     &            NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &            MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &            REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((ERR) .OR. (NOROOM)) GOTO 120
               NEW2 = NNN - 1
               NEW3 = NNN
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL D2NODE (MXND, XN, YN, NEW1, NEW2)
                  CALL D2NODE (MXND, XN, YN, NEW2, NEW3)
                  CALL D2NODE (MXND, XN, YN, N2, NEW3)
                  CALL SFLUSH
                  IF (VIDEO) THEN
                     CALL SNAPIT (1)
                  ENDIF
               ENDIF
C
               N1 = N2
               GOTO 100
            ELSE
               NBEGIN = NADJ1
               NADJ1 = 0
               NADJ2 = 0
               GOTO 120
            ENDIF
C
C  IF N2 IS A ROW REVERSAL NODE,
C  THEN PROJECT THREE NODES (1, 3, & 5) AT AVERAGE ANGLES AND
C  TWO OTHERS (2 & 4) AS AN ISOPARAMETRIC TYPE
C
         ELSEIF (LNODES (1, N2) .EQ. 7) THEN
            AHOLD = ANGLE (N2)
            CALL EXTND5 (MXND, XN, YN, AHOLD, N1, N2, N3, X, Y, DIST)
C
            CALL ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN, ANGLE,
     &         BNSIZE, LNODES, X(1), Y(1), DIST, NNN, KKK, LLL,
     &         N0, N1, N2, AMBIG, LAMBIG, SIZEIT, ERR, NOROOM, XNOLD,
     &         YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD,
     &         NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN)
            IF ((ERR) .OR. (NOROOM)) GOTO 120
            AMBIG = .FALSE.
            NEW1 = NNN
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL D2NODE (MXND, XN, YN, N0, NEW1)
               CALL D2NODE (MXND, XN, YN, N2, NEW1)
               CALL SFLUSH
               IF (VIDEO) THEN
                  CALL SNAPIT (1)
               ENDIF
            ENDIF
C
C  TRY A PINCH TO MAKE SURE THAT CONTINUING DOES NOT OVERLAY
C  AN ALREADY GOOD CLOSURE
C
            CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &         NNNOLD, LLLOLD, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) GOTO 120
            LLLOLD = LLL
            NNNOLD = NNN
            CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         NLOOP (1), ANGLE, LNODES, NADJ1, LLL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
            IF (ERR) GOTO 120
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP (1), XMIN, XMAX,
     &         YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
C
            CALL GETIME (TIME1)
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
               IF (VIDEO) CALL SNAPIT (2)
            ENDIF
C
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL PINCH (MXND, MXCORN, MLN, NUID, XN, YN, ZN, LXK, KXL,
     &         NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1, NLOOP (1),
     &         KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &         XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &         LCORN, NCORN, IDUM1, IDUM2, GRAPH, VIDEO, KREG, NOROOM,
     &         ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 120
            IF (DONE) GOTO 120
            CALL GETIME (TIME1)
C
C  TRY A COLAPS TO MAKE SURE THAT CONTINUING DOES NOT COMPLICATE
C  AN ALREADY OVERLAPPED SIDE
C
            IF (NLOOP (1) .GT. 6) THEN
               NLOLD = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               CALL COLAPS (MXND, MXCORN, MLN, MXLOOP, NUID, XN, YN,
     &            ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1,
     &            KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &            XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &            LCORN, NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO,
     &            KREG, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               IF (DONE) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  CHECK TO SEE IF ANY OF THE CONCURRENT PERIMETERS OVERLAP
C
            IF (LINKPR (2, KPERIM) .NE. 0) THEN
               LINKPR (3, KPERIM) = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               IDUM = NADJ1
               CALL PCROSS (MXND, MXCORN, MLN, MXLOOP, MAXPRM,
     &            NUID, XN, YN, ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES,
     &            BNSIZE, LINKPR, KPERIM, IDUM, NADJ1, NEW1, KKKOLD,
     &            LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN, LCORN,
     &            NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO, KREG,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  THE ROW CAN BE CONTINUED IF THE NEW NODE AND THE ENDING NODE
C  IS STILL ON THE BOUNDARY AFTER THE PINCH, COLAPS AND PCROSS
C  AND IF THE ANGLE AT THE NODE IS SOMEWHAT AS IT SHOULD BE
C
            IF ( (NLOLD .EQ. NLOOP(1)) .AND.
     &         (IABS (LNODES (4, NEW1)) .EQ. 1) .AND.
     &         (LNODES (2, N2) .EQ. NEW1) .AND.
     &         (IABS (LNODES (4, NEND)) .EQ. 1) .AND.
     &         (LXN (1, NEW1) .GT. 0) .AND.
     &         (ANGLE (N2) .GT. 3.9269908 ) ) THEN
               CALL EXTND5 (MXND, XN, YN, AHOLD, N1, N2, N3, X, Y, DIST)
               CALL ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &            BNSIZE, LNODES, X(2), Y(2), X(3), Y(3), DIST, DIST,
     &            NNN, LLL, KKK, NEW1, N2, NLOOP (1), SIZEIT, ERR,
     &            NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &            MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &            REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((ERR) .OR. (NOROOM)) GOTO 120
               NEW2 = NNN - 1
               NEW3 = NNN
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL D2NODE (MXND, XN, YN, NEW1, NEW2)
                  CALL D2NODE (MXND, XN, YN, NEW2, NEW3)
                  CALL D2NODE (MXND, XN, YN, N2, NEW3)
                  CALL SFLUSH
                  IF (VIDEO) THEN
                     CALL SNAPIT (1)
                  ENDIF
               ENDIF
            ELSE
               NBEGIN = NADJ1
               NADJ1 = 0
               NADJ2 = 0
               GOTO 120
            ENDIF
C
C  TRY A PINCH TO MAKE SURE THAT CONTINUING DOES NOT OVERLAY
C  AN ALREADY GOOD CLOSURE
C
            CALL FIXLXN (MXND, LXN, NXL, NUID, NAVAIL, IAVAIL, NNN, LLL,
     &         NNNOLD, LLLOLD, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) GOTO 120
            LLLOLD = LLL
            NNNOLD = NNN
            CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         NLOOP (1), ANGLE, LNODES, NADJ1, LLL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
            IF (ERR) GOTO 120
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP (1), XMIN, XMAX,
     &         YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
C
            CALL GETIME (TIME1)
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
               IF (VIDEO) CALL SNAPIT (2)
            ENDIF
C
            CALL GETIME (TIME2)
            TIMEA = TIMEA + TIME2 - TIME1
            CALL PINCH (MXND, MXCORN, MLN, NUID, XN, YN, ZN, LXK, KXL,
     &         NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1, NLOOP (1),
     &         KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &         XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &         LCORN, NCORN, IDUM1, IDUM2, GRAPH, VIDEO, KREG, NOROOM,
     &         ERR)
            IF ((NOROOM) .OR. (ERR)) GOTO 120
            IF (DONE) GOTO 120
            CALL GETIME (TIME1)
C
C  TRY A COLAPS TO MAKE SURE THAT CONTINUING DOES NOT COMPLICATE
C  AN ALREADY OVERLAPPED SIDE
C
            IF (NLOOP (1) .GT. 6) THEN
               NLOLD = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               CALL COLAPS (MXND, MXCORN, MLN, MXLOOP, NUID, XN, YN,
     &            ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES, BNSIZE, NADJ1,
     &            KKKOLD, LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN,
     &            XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN,
     &            LCORN, NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO,
     &            KREG, NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               IF (DONE) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  CHECK TO SEE IF ANY OF THE CONCURRENT PERIMETERS OVERLAP
C
            IF (LINKPR (2, KPERIM) .NE. 0) THEN
               LINKPR (3, KPERIM) = NLOOP (1)
               CALL GETIME (TIME2)
               TIMEA = TIMEA + TIME2 - TIME1
               IDUM = NADJ1
               CALL PCROSS (MXND, MXCORN, MLN, MXLOOP, MAXPRM,
     &            NUID, XN, YN, ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES,
     &            BNSIZE, LINKPR, KPERIM, IDUM, NADJ1, NEW3, KKKOLD,
     &            LLLOLD, NNNOLD, IAVAIL, NAVAIL, DONE, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL, KKK, NNN, LCORN,
     &            NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO, KREG,
     &            NOROOM, ERR)
               IF ((NOROOM) .OR. (ERR)) GOTO 120
               CALL GETIME (TIME1)
            ENDIF
C
C  THE ROW CAN BE CONTINUED IF THE NEW NODE AND THE ENDING NODE
C  IS STILL ON THE BOUNDARY AFTER THE PINCH, COLAPS AND PCROSS
C
            IF ( (NLOLD .EQ. NLOOP(1)) .AND.
     &         (IABS (LNODES (4, NEW3)) .EQ. 1) .AND.
     &         (LNODES (2, N2) .EQ. NEW3) .AND.
     &         (IABS (LNODES (4, NEND)) .EQ. 1) .AND.
     &         (LXN (1, NEW3) .GT. 0) .AND.
     &         (ANGLE (N2) .GT. 2.3561945) ) THEN
               CALL EXTND5 (MXND, XN, YN, AHOLD, N1, N2, N3, X, Y, DIST)
               CALL ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &            BNSIZE, LNODES, X(4), Y(4), X(5), Y(5), DIST, DIST,
     &            NNN, LLL, KKK, NEW3, N2, NLOOP (1), SIZEIT, ERR,
     &            NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &            MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &            REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
               IF ((ERR) .OR. (NOROOM)) GOTO 120
               NEW4 = NNN - 1
               IF ((GRAPH) .OR. (VIDEO)) THEN
                  CALL D2NODE (MXND, XN, YN, NEW3, NEW4)
                  CALL D2NODE (MXND, XN, YN, NEW4, NNN)
                  CALL D2NODE (MXND, XN, YN, N2, NNN)
                  CALL SFLUSH
                  IF (VIDEO) THEN
                     CALL SNAPIT (1)
                  ENDIF
               ENDIF
            ELSE
               NBEGIN = NADJ1
               NADJ1 = 0
               NADJ2 = 0
               GOTO 120
            ENDIF
C
            N1 = N2
            GOTO 100
C
         ELSE
            CALL MESAGE ('PROBLEMS IN ADDROW - THE NEXT NODE IN THE')
            CALL MESAGE ('ROW DOES NOT HAVE THE RIGHT CLASSIFICATION')
            ERR = .TRUE.
            GOTO 120
         ENDIF
C
         IF (KOUNT .GT. NLOOP (1)) THEN
            CALL MESAGE ('PROBLEMS IN ADDROW - THE ROW DOESN''T STOP')
            ERR = .TRUE.
            GOTO 120
         ENDIF
C
C  CLOSE OFF THE END OF THE ROW - THE END NODE IS THE CONNECTION
C
      ELSE
         NADJ2 = LNODES (3, N2)
         CALL CONNOD (MXND, MLN, XN, YN, NUID, LXK, KXL, NXL, LXN,
     &      ANGLE, LNODES, NNN, KKK, LLL, NNNOLD, LLLOLD, N0, N1, N2,
     &      LNODES (3, N2), N1, NLOOP (1), IAVAIL, NAVAIL, GRAPH,
     &      VIDEO, NOROOM, ERR)
         IF ((NOROOM) .OR. (ERR)) GOTO 120
C
C  NOW SMOOTH THE CURRENT ROW INTO THE MESH
C
         CALL GETIME (TIME2)
         TIMEA = TIMEA + TIME2 - TIME1
         CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &      LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP (1), XMIN, XMAX,
     &      YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
         CALL GETIME (TIME1)
C
C  CALCULATE NEW ANGLES
C
         CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &      NLOOP (1), ANGLE, LNODES, N3, LLL, XMIN, XMAX, YMIN, YMAX,
     &      ZMIN, ZMAX, DEV1, KREG, ERR)
         IF (ERR) GOTO 120
C
C  PLOT THE NEW MESH
C
         IF (VIDEO) THEN
            CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
            CALL SNAPIT (2)
         ENDIF
      ENDIF
      NBEGIN = N3
C
  120 CONTINUE
      CALL GETIME (TIME2)
      TIMEA = TIMEA + TIME2 - TIME1
C
      RETURN
C
      END
