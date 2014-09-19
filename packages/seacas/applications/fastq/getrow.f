C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: getrow.f,v 1.3 1998/07/14 18:19:05 gdsjaar Exp $
C $Log: getrow.f,v $
C Revision 1.3  1998/07/14 18:19:05  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1995/06/28 19:21:18  gdsjaar
C Applied fixes found in memo dated May 13, 1991. The bug shows itself
C for rare cases of semicircular regions being paved.
C
c Revision 1.1.1.1  1990/11/30  11:08:41  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:08:39  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]GETROW.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETROW (MXND, MXCORN, MXPICK, MLN, NUID, LXK, KXL,
     &   NXL, LXN, LNODES, NCORN, LCORN, BNSIZE, ANGLE, XN, YN, ZN,
     &   ICOMB, ITYPE, NLOOP, NBEGIN, NEND, IAVAIL, NAVAIL, LLL, KKK,
     &   NNN, GRAPH, VIDEO, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1,
     &   KREG, SIZEIT, NEXTPR, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE GETROW = GETS THE CURRENT ROW TO START ON
C
C***********************************************************************
C
      DIMENSION LNODES(MLN, MXND), LCORN(MXCORN), ANGLE(MXND)
      DIMENSION BNSIZE(2, MXND)
      DIMENSION NUID(MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION NXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION ICOMB(MXCORN, MXPICK), ITYPE(MXPICK)
      DIMENSION ITEST(5), LTEST(5)
      DIMENSION IPINCH(4), JPINCH(4)
      DIMENSION XN(MXND), YN(MXND), ZN(MXND), X(1), Y(1)
C
      CHARACTER*3 DEV1
C
      LOGICAL POSBL2, POSBL3, POSBL4
      LOGICAL FOUND2, FOUND3, FOUND4
      LOGICAL GRAPH, ONLYC, CORNP, REDO2, REDO3, PPOSBL, VIDEO
      LOGICAL SIDPIN, ROWCHN, ERR, SIZEIT, NOROOM
C
      ERR = .FALSE.
      ONLYC = .FALSE.
      NPIN2 = 0
C
C  RESET EVERYTHING TO BE FREE
C
  100 CONTINUE
      IF (GRAPH) CALL LCOLOR ('YELOW')
      CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 0, ERR)
      IF (ERR) GOTO 270
C
C  GET THE CURRENT CORNERS
C
      CALL GETCRN (MXND, MXCORN, MLN, LNODES, NCORN, LCORN,
     &   ANGLE, XN, YN, LXN, NLOOP, NBEGIN, ONLYC, PPOSBL, GRAPH, ERR)
      IF (ERR) GOTO 270
C
C  GET ALL THE COMBINATIONS - NPICK IS THE NUMBER OF COMBINATIONS
C
      IF (PPOSBL) THEN
         CALL COMSRT (MXND, MXCORN, MXPICK, MLN, LNODES, LCORN, NCORN,
     &      ICOMB, ITYPE, NPICK)
      ELSE
         NPICK = 0
      ENDIF
C
C  NOW CHECK FOR THE STANDARD PRIMITIVE TYPES
C
      FOUND4 = .FALSE.
      FOUND3 = .FALSE.
      FOUND2 = .FALSE.
      REDO2 = .TRUE.
      REDO3 = .TRUE.
      SIDPIN = .FALSE.
C
C  SKIP THE PRIMITIVES IF SIZEIT IS IN EFFECT
C
      IF ((SIZEIT) .OR. (NEXTPR .NE. 0)) GOTO 130
C
C  SET UP THE MINIMUM ACCEPTABLE QUALITIES
C
      BEST2 = 3.
      BEST3 = 3.
      BEST4 = 4.
C
      DO 120 I = 1, NPICK
C
C  NOW GET THE BEST RECTANGLE COMBINATION
C
         IF (ITYPE(I) .LT. 5) THEN
            IF (ITYPE(I) .EQ. 4) THEN
               CALL QUAL4 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES,
     &            ICOMB(1, I), ANGLE, LXN, ITEST, LTEST, QUAL, POSBL4,
     &            ERR)
               IF (ERR) GOTO 270
C
C  GET THE RECTANGLE INTERPRETATION
C
               IF (POSBL4) THEN
                  IF (QUAL .LT. BEST4) THEN
                     IS2C = 0
                     IBEST4 = I
                     BEST4 = QUAL
                     FOUND4 = .TRUE.
                  ENDIF
               ENDIF
C
C  NOW GET THE BEST TRIANGLE COMBINATION
C
            ELSEIF (ITYPE(I) .EQ. 3) THEN
               CALL QUAL3 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES,
     &            ICOMB(1, I), ANGLE, LXN, ITEST, LTEST, QUAL, POSBL3,
     &            POSBL4, ERR)
               IF (ERR) GOTO 270
C
C  GET THE PURE TRIANGLE INTERPRETATION
C
               IF (POSBL3) THEN
                  IF (QUAL .LT. BEST3) THEN
                     IBEST3 = I
                     BEST3 = QUAL
                     FOUND3 = .TRUE.
                     REDO3 = .TRUE.
                  ENDIF
C
C  TRY A CHANGE TO A RECTANGLE - OR RATHER SET THE RIGHT ROW FOR
C  THE PROPER CONCLUSION OF A TRIANGLE
C
               ELSEIF (POSBL4) THEN
C
C  MAKE SURE THAT THE RESULTING CORNER IS NOT ON A BAD SIDE
C
                  CALL CH3TO4 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES,
     &               ICOMB(1, I), ANGLE, ITEST, LTEST, QUAL, POSBL4,
     &               ICHNG)
                  IF ((LXN(2, ICHNG) .GT. 0) .OR.
     &               ((LXN(2, ICHNG) .LT. 0) .AND.
     &               (CORNP(ANGLE(ICHNG))))) THEN
C
C  SEE IF WE SHOULD KEEP IT BASED ON QUALITY
C
                     IF (QUAL .LT. BEST3) THEN
                        IBEST3 = I
                        BEST3 = QUAL
                        FOUND3 = .TRUE.
                        REDO3 = .FALSE.
                        NBEG34 = NBEGIN
                        CALL KEEP3 (ITEST, LTEST, NBEG34, NEND34)
                     ENDIF
                  ENDIF
C                  IF ((POSBL4) .AND.
C     &               (QUAL .LT. BEST4)) THEN
C                     IS2C = II
C                     IBEST4 = I
C                     BEST4 = QUAL
C                     FOUND4 = .TRUE.
C                  ENDIF
               ENDIF
C
C  NOW GET THE BEST SEMICIRCLE COMBINATION
C
            ELSEIF (ITYPE(I) .EQ. 2) THEN
               CALL QUAL2N (MXND, MXCORN, MLN, NCORN, LCORN, LNODES,
     &            ICOMB(1, I), BNSIZE, ANGLE, LXN, ITEST, LTEST, QUAL,
     &            POSBL2, POSBL3, ROWCHN, SIDPIN, ISTART, IEND, JPINCH,
     &            NPIN2, ERR)
               IF (ERR) GOTO 270
C               CALL QUAL2 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES,
C     &            ICOMB(1, I), BNSIZE, ANGLE, LXN, ITEST, LTEST, QUAL,
C     &            POSBL2, POSBL3, ROWCHN, ISTART, IEND)
C
C  GET THE SEMICIRCLE INTERPRETATION
C
               IF (POSBL2) THEN
                  IF (QUAL .LT. BEST2) THEN
                     IS2C = 0
                     IBEST2 = I
                     BEST2 = QUAL
                     FOUND2 = .TRUE.
                     IF (SIDPIN) THEN
                        NPINCH = NPIN2
                        DO 110 IN = 1, NPINCH
                           IPINCH(IN) = JPINCH(IN)
  110                   CONTINUE
                     ELSE
                        SIDPIN = .FALSE.
                        IF (ROWCHN) THEN
                           NBEG24 = ISTART
                           NEND24 = IEND
                           REDO2 = .FALSE.
                        ELSE
                           REDO2 = .TRUE.
                        ENDIF
                     ENDIF
                  ENDIF
C
C  TRY A CHANGE TO A RECTANGLE - OR RATHER SET THE RIGHT ROW FOR
C  THE PROPER CONCLUSION OF A TRIANGLE
C
               ELSEIF (POSBL3) THEN
C
C  MAKE SURE THAT THE RESULTING CORNER IS NOT ON A BAD SIDE
C
                  IF ((LXN(2, IEND) .GT. 0) .OR.
     &               ((LXN(2, IEND) .LT. 0) .AND.
     &               (CORNP (ANGLE(IEND))))) THEN
C
C  SEE IF WE SHOULD KEEP IT BASED ON QUALITY
C
                     IF (QUAL .LT. BEST3) THEN
                        IBEST3 = I
                        BEST3 = QUAL
                        FOUND3 = .TRUE.
                        REDO3 = .FALSE.
                        NBEG34 = ISTART
                        NEND34 = IEND
                     ENDIF
                  ENDIF
               ENDIF
C
            ENDIF
         ENDIF
  120 CONTINUE
  130 CONTINUE
C
C  FOR NOW, THE RECTANGLE WILL ALWAYS WIN, ETC.
C
C  TAKE THE RECTANGLE
C
      IF (FOUND4) THEN
         CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 3, ERR)
         IF (ERR) GOTO 270
         DO 140 I = 1, NCORN
            IF (ICOMB(I, IBEST4) .EQ. 1) THEN
               LNODES(1, LCORN(I)) = 1
               IF (GRAPH) THEN
                  ISQR = LCORN(I)
                  X(1) = XN(ISQR)
                  Y(1) = YN(ISQR)
                  CALL SYMBOL (1, X, Y, 'SQUARE')
                  CALL SFLUSH
               ENDIF
            ENDIF
  140    CONTINUE
         IF (IS2C .GT. 0) THEN
            NBEGIN = IS2C
            LNODES(1, NBEGIN) = 1
            IF (GRAPH) THEN
               ISQR = NBEGIN
               X(1) = XN(ISQR)
               Y(1) = YN(ISQR)
               CALL SYMBOL (1, X, Y, 'SQUARE')
               CALL SFLUSH
            ENDIF
         ELSE
            DO 150 I = 1, NCORN
               IF (ICOMB (I, IBEST4) .EQ. 1) THEN
                  NBEGIN = LCORN(I)
                  II = I
                  GOTO 160
               ENDIF
  150       CONTINUE
            NBEGIN = LCORN(1)
  160       CONTINUE
         ENDIF
C
C  TAKE THE TRIANGLE
C
      ELSEIF (FOUND3) THEN
         CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 3, ERR)
         IF (ERR) GOTO 270
         DO 170 I = 1, NCORN
            IF (ICOMB(I, IBEST3) .EQ. 1) THEN
               LNODES(1, LCORN(I)) = 1
               IF (GRAPH) THEN
                  ISQR = LCORN(I)
                  X(1) = XN(ISQR)
                  Y(1) = YN(ISQR)
                  CALL SYMBOL (1, X, Y, 'SQUARE')
                  CALL SFLUSH
               ENDIF
            ENDIF
  170    CONTINUE
         IF (REDO3) THEN
            DO 180 I = 1, NCORN
               IF (ICOMB(I, IBEST3) .EQ. 1) THEN
                  NBEGIN = LCORN(I)
                  II = I
                  GOTO 190
               ENDIF
  180       CONTINUE
            NBEGIN = LCORN(1)
            II = IBEST3
         ELSE
            NBEGIN = NBEG34
            NEND = NEND34
            LNODES(1, NBEGIN) = 1
            LNODES(1, NEND) = 1
            GOTO 260
         ENDIF
  190    CONTINUE
C
C  OTHERWISE TAKE THE SEMICIRCLE
C
      ELSEIF (FOUND2) THEN
C
C  IF THE BEST SEMICIRCLE MUST BE TUCKED, THEN DO SO AND THEN
C  REDO THE WHOLE SORTING - A RECTANGLE SHOULD RESULT
C
         IF (SIDPIN) THEN
            DO 200 I = 1, NPINCH
               NBEGIN = LNODES(2, IPINCH(I))
C
C  MAR   K THE SMOOTHING
C
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, IPINCH(I)), ERR)
               IF (ERR) GOTO 270
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (2, LNODES (2, IPINCH(I))), ERR)
               IF (ERR) GOTO 270
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, IPINCH(I)), ERR)
               IF (ERR) GOTO 270
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, IPINCH(I))), ERR)
               IF (ERR) GOTO 270
               CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &            LNODES (3, LNODES (3, LNODES (3, IPINCH(I)))), ERR)
               IF (ERR) GOTO 270
C
               CALL TUCK (MXND, MLN, NUID, XN, YN, LXK, KXL, NXL, LXN,
     &            LNODES, IAVAIL, NAVAIL, LLL, KKK, NNN, IPINCH(I),
     &            NLOOP, GRAPH, NOROOM, ERR)
               IF (ERR) GOTO 270
               IF (VIDEO) THEN
                  CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &               YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                  CALL SNAPIT (1)
               ENDIF
  200       CONTINUE
            NNN2 = 1
            CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, DEV1, KREG)
            CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &         NLOOP, ANGLE, LNODES, NBEGIN, LLL, XMIN, XMAX, YMIN,
     &         YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
            IF (ERR) GOTO 270
            IF ((GRAPH) .OR. (VIDEO)) THEN
               CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX, YMIN,
     &            YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
               IF (VIDEO) CALL SNAPIT (1)
            ENDIF
            GOTO 100
         ELSE
            CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 3, ERR)
            IF (ERR) GOTO 270
            DO 210 I = 1, NCORN
               IF (ICOMB (I, IBEST2) .EQ. 1) THEN
                  LNODES(1, LCORN(I)) = 1
                  IF (GRAPH) THEN
                     ISQR = LCORN(I)
                     X(1) = XN(ISQR)
                     Y(1) = YN(ISQR)
                     CALL SYMBOL (1, X, Y, 'SQUARE')
                     CALL SFLUSH
                  ENDIF
               ENDIF
  210       CONTINUE
            IF (REDO2) THEN
               DO 220 I = 1, NCORN
                  IF (ICOMB(I, IBEST2) .EQ. 1) THEN
                     NBEGIN = LCORN(I)
                     II = I
                     GOTO 230
                  ENDIF
  220          CONTINUE
               NBEGIN = LCORN(1)
            ELSE
               NBEGIN = NBEG24
               NEND = NEND24
               LNODES(1, NBEGIN) = 1
               LNODES(1, NEND) = 1
               GOTO 260
            ENDIF
  230       CONTINUE
         ENDIF
C
C  CHECK FOR A ONE SIDED SEMICIRCLE
C
      ELSEIF ( (NCORN .EQ. 2) .AND. ((LNODES (7, LCORN(1)) .EQ. 1) .OR.
     &   (LNODES (7, LCORN(2)) .EQ. 1)) ) THEN
         IF (LNODES (7, LCORN(1)) .EQ. 1) THEN
            NBEGIN = LCORN (1)
            NEND = LCORN (2)
         ELSE
            NBEGIN = LCORN (2)
            NEND = LCORN (1)
         ENDIF
         GOTO 260
C
C  OTHERWISE, THE DEFAULT IS TO JUST START AT THE NEXT CORNER
C
      ELSEIF (NCORN .GT. 0) THEN
         CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 3, ERR)
         IF (ERR) GOTO 270
         II = 0
         NEND = 0
         INODE = NBEGIN
         DO 240 I = 1, NLOOP + 1
            CALL NDSTAT (INODE, LXN(1, INODE), ANGLE(INODE), ISTAT)
            LNODES(1, INODE) = ISTAT
C
C  SAVE THE FIRST NATURAL CORNER AS THE START
C
            IF ((II .EQ. 0) .AND. (ISTAT .EQ. 1)) THEN
               NBEGIN = INODE
               II = 1
C
C  A ROW END HAS BEEN FOUND
C
            ELSEIF (ISTAT .EQ. 1) THEN
               NEND = INODE
               GOTO 260
            ENDIF
            INODE = LNODES(3, INODE)
  240    CONTINUE
C
C  THE ROW IS A CLOSED LOOP BACK TO THE SAME CORNER
C
         IF ((II .NE. 0) .AND. (NEND .EQ. 0)) THEN
            NEND = NBEGIN
C
C  THE ROW DOESN'T CONTAIN ANY TRUE CORNERS - TREAT IT AS A CIRCLE
C
         ELSE
            CALL SETCIR (MXND, MLN, NLOOP, LNODES, NBEGIN, ERR)
            NEND = NBEGIN
         ENDIF
         GOTO 260
C
C  NO CORNERS - JUST SET EVERYTHING TO BE A SIDE
C
      ELSE
         CALL SETLOP (MXND, MLN, NLOOP, LNODES, NBEGIN, 3, ERR)
         IF (ERR) GOTO 270
         NEND = NBEGIN
         GOTO 260
      ENDIF
C
C  FIND THE NEXT NATURAL CORNER
C
      NEND = NBEGIN
      JJ = 0
      DO 250 I = 1, NCORN
         INODE = LCORN(I)
         IF ((LNODES(1, INODE) .EQ. 1) .AND.
     &      (INODE .NE. NBEGIN)) THEN
            IF (((JJ .EQ. 0) .AND. (I .LT. II)) .OR.
     &         ((JJ .LT. II) .AND. (I .GT. II))) THEN
               JJ = I
               NEND = INODE
            ENDIF
         ENDIF
  250 CONTINUE
C
  260 CONTINUE
C
      IF (GRAPH) THEN
C  5 IS PINK; 4 IS BLUE; 3 IS YELLOW; 0 IS BLACK ; 7 IS WHITE; 1 IS RED
         CALL LCOLOR ('PINK ')
         X(1) = XN(NBEGIN)
         Y(1) = YN(NBEGIN)
         CALL SYMBOL (1, X, Y, 'SQUARE')
         X(1) = XN(NEND)
         Y(1) = YN(NEND)
         CALL SYMBOL (1, X, Y, 'SQUARE')
         CALL SFLUSH
         CALL LCOLOR ('WHITE')
      ENDIF
C
  270 CONTINUE
C
      RETURN
C
      END
