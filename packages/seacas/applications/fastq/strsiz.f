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

C $Id: strsiz.f,v 1.4 1999/06/21 22:43:41 gdsjaar Exp $
C $Log: strsiz.f,v $
C Revision 1.4  1999/06/21 22:43:41  gdsjaar
C Fixed more uninitialized variables; one was causing core dump on g77
C compiled executable.
C
C VERSN was not consistently defined -- now 10 characters everywhere
C
C Updated so full version string output
C
C Added capability to debug memory using unit specified in EXT99
C variable. Similar to STRTUP in SUPLIB
C
C Cleaned up some other code
C
C Upped version
C
C Revision 1.3  1998/07/14 18:20:05  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/04/10 19:56:59  gdsjaar
C Fixed some logical variables
C
c Revision 1.1.1.1  1990/11/30  11:16:50  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:16:49  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]STRSIZ.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO STRSIZ TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
      SUBROUTINE STRSIZ (MAXNP, X, Y, NINT, N, XEND, YEND, XDIFF, YDIFF,
     &   D, ERR, TEST, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &   MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &   REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, GRAPH, DXMAX)
C***********************************************************************
C
C  SUBROUTINE STRSIZ = GETS INTERVALS ON A SRAIGHT LINE BASED ON ERROR
C                      SIZE
C
C***********************************************************************
C
      DIMENSION X (MAXNP), Y (MAXNP)
C
      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD), NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)
C
      LOGICAL GRAPH, REMESH, TEST, ERR, SGRAPH, MOVED
C
      IF (GRAPH) THEN
         SGRAPH = .TRUE.
      ELSE
         SGRAPH = .FALSE.
      ENDIF
C
      ITERAT = 100
      EPS = .01
C
      DSTNOW = 0.
      INTNOW = 0
      SIZNOW = 0.0
      IF (GRAPH) THEN
         CALL SYMBOL (1, X(1), Y(1), 'DIAMND')
         CALL PLTFLU
      ENDIF
  100 CONTINUE
      INTNOW = INTNOW + 1
      IF (DSTNOW + (SIZNOW * 1.3) .GT. D)THEN
C
C  THE END OF THE LINE (OR CLOSE ENOUGH) HAS BEEN REACHED
C
C  IF WE ARE TESTING OR THE INTERVALS MATCH, THEN SIMPLY FINISH THE
C  LINE.
C
         IF ((TEST) .OR. (INTNOW .EQ. NINT)) THEN
            NINT = INTNOW
            N = NINT + 1
            X(N) = XEND
            Y(N) = YEND
            IF (GRAPH) THEN
               CALL SYMBOL (1, X(INTNOW), Y(INTNOW), 'DIAMND')
               CALL MPD2VC (1, X(INTNOW), Y(INTNOW),
     &            X(N), Y(N))
               CALL SYMBOL (1, X(N), Y(N), 'DIAMND')
               CALL PLTFLU
            ENDIF
         ELSE
C
C  OTHERWISE, MAKE SURE THE INTERVALS ARE ALRIGHT AND ADD THE EXTRA ONE
C
            EPS = .001
            IF (INTNOW + 1 .NE. NINT) THEN
               CALL MESAGE ('** PROBLEMS WITH INTNOW '//
     &            'IN PLINE **')
               ERR = .TRUE.
               GOTO 160
            ENDIF
            X(INTNOW + 1) = (X(INTNOW) + XEND) * .5
            Y(INTNOW + 1) = (Y(INTNOW) + YEND) * .5
            N = NINT + 1
            X(N) = XEND
            Y(N) = YEND
            IF (GRAPH) THEN
               CALL SYMBOL (1, X(INTNOW + 1), Y(INTNOW + 1),
     &            'DIAMND')
               CALL SYMBOL (1, X(N), Y(N),
     &            'DIAMND')
               CALL MPD2VC (1, X(INTNOW), Y(INTNOW),
     &            X(INTNOW+1), Y(INTNOW+1))
               CALL MPD2VC (1, X(INTNOW+1), Y(INTNOW+1),
     &            X(N), Y(N))
               CALL PLTFLU
            ENDIF
         ENDIF
      ELSE
C
C  NOT TO THE END YET
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO CALL TO GETSIZ TO PASS MINIMUM
CC**              ELEMENT SIZE (SIZMIN) AND GETSIZ PARAMETERS OF
CC**              EMIN AND EMAX
C
         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &      MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &      X(INTNOW), Y(INTNOW), S1)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 8/2/90
CC* MODIFICATION: ADDED A SIZE ADJUSTMENT BASED ON THE REQUIRED VALUE
CC*               AT THE END OF THE SEGMENT AND AT THE AVERAGE OF THE
CC*               SEGMENTS - THE 2ND AND 3RD CALL TO GETSIZ.
         XNEW1 = X(1) + (((DSTNOW + S1) / D) * XDIFF)
         YNEW1 = Y(1) + (((DSTNOW + S1) / D) * YDIFF)
         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &      MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &      XNEW1, YNEW1, S2)
         XNEW2 = X(1) + (((DSTNOW + ((S1+S2) * .5)) / D) * XDIFF)
         YNEW2 = Y(1) + (((DSTNOW + ((S1+S2) * .5)) / D) * YDIFF)
         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &      MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &      XNEW2, YNEW2, S3)
C
         SIZNOW = (((S1 + S2) * .5) + S3) * .5
C
         DSTNOW = DSTNOW + SIZNOW
         X(INTNOW + 1) = X(1) + ((DSTNOW / D) * XDIFF)
         Y(INTNOW + 1) = Y(1) + ((DSTNOW / D) * YDIFF)
         IF (GRAPH) THEN
            CALL SYMBOL (1, X(INTNOW + 1), Y(INTNOW + 1),
     &         'DIAMND')
            CALL MPD2VC (1, X(INTNOW), Y(INTNOW),
     &         X(INTNOW+1), Y(INTNOW+1))
            CALL PLTFLU
         ENDIF
         GOTO 100
      ENDIF
C
C  ERASE THE NODES FOR REPLOTTING IF NEEDED
C
      IF ((.NOT. SGRAPH) .AND. (GRAPH)) THEN
         DO 110 J = 2, NINT
            CALL LCOLOR ('BLACK')
            CALL MPD2VC (1, X(J), Y(J), X(J+1), Y(J+1))
            CALL MPD2VC (1, X(J), Y(J), X(J-1), Y(J-1))
            CALL SYMBOL (1, X(J), Y(J), 'DIAMND')
            CALL LCOLOR ('WHITE')
            CALL PLTFLU
  110    CONTINUE
      ENDIF
C
C  NOW SMOOTH THE NODES ALONG THE LINE
C
      DO 130 I = 1, ITERAT
         MOVED = .FALSE.
         DO 120 J = 2, NINT
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG,
     &         LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK,
     &         REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN, X(J-1), Y(J-1), SIZE1)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG,
     &         LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK,
     &         REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN, X(J), Y(J), SIZE2)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG,
     &         LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK,
     &         REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &         IDIVIS, SIZMIN, EMAX, EMIN, X(J+1), Y(J+1), SIZE3)
            DIST1 = SQRT ( ((X(J-1) - X(J)) **2) +
     &         ((Y(J-1) - Y(J)) **2) )
            DIST2 = SQRT ( ((X(J) - X(J+1)) **2) +
     &         ((Y(J) - Y(J+1)) **2) )
            DTOTAL = DIST1 + DIST2
            RATIO = DIST1 / DTOTAL
            DRATIO = ((SIZE1 + SIZE2) * .5) /
     &         ( ((SIZE1 + SIZE2) * .5) +
     &         ((SIZE2 + SIZE3) * .5) )
            TRATIO = (RATIO + DRATIO) * .5
            IF (SGRAPH) THEN
               CALL LCOLOR ('BLACK')
               CALL MPD2VC (1, X(J), Y(J), X(J+1), Y(J+1))
               CALL MPD2VC (1, X(J), Y(J), X(J-1), Y(J-1))
               CALL SYMBOL (1, X(J), Y(J), 'DIAMND')
               CALL LCOLOR ('WHITE')
               CALL PLTFLU
            ENDIF
            X(J) = (TRATIO * (X(J+1) - X(J-1))) + X(J-1)
            Y(J) = (TRATIO * (Y(J+1) - Y(J-1))) + Y(J-1)
            IF (SGRAPH) THEN
               CALL MPD2VC (1, X(J), Y(J), X(J+1), Y(J+1))
               CALL MPD2VC (1, X(J), Y(J), X(J-1), Y(J-1))
               CALL SYMBOL (1, X(J), Y(J), 'DIAMND')
               CALL PLTFLU
            ENDIF
            DX1 = DIST1 / (.5 * (SIZE1 + SIZE2))
            DX2 = DIST2 / (.5 * (SIZE2 + SIZE3))
            IF (J .EQ. 2) THEN
               DXMAX = AMAX1 (DX1, DX2)
            ELSE
               DXMAX = AMAX1 (DXMAX, DX1, DX2)
            ENDIF
            DT = ABS((TRATIO * DTOTAL) - DIST1)
            IF (DT/DTOTAL .GT. EPS) MOVED = .TRUE.
  120    CONTINUE
         IF (.NOT. MOVED) GOTO 140
  130 CONTINUE
  140 CONTINUE
      IF ((.NOT. SGRAPH) .AND. (GRAPH)) THEN
         DO 150 J = 2, NINT
            CALL MPD2VC (1, X(J), Y(J), X(J+1), Y(J+1))
            CALL MPD2VC (1, X(J), Y(J), X(J-1), Y(J-1))
            CALL SYMBOL (1, X(J), Y(J), 'DIAMND')
            CALL PLTFLU
  150    CONTINUE
      ENDIF
C
  160 CONTINUE
C
      RETURN
C
      END
