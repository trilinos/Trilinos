C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CRCSIZ (MAXNP, X, Y, NINT, N, XEND, YEND, XCEN, YCEN,
     &   THETA1, THETA2, TANG, AA, BB, ERR, TEST, XNOLD, YNOLD, NXKOLD,
     &   LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &   GRAPH)
C***********************************************************************

C  SUBROUTINE CRCSIZ = GETS INTERVALS ON AN ARC LINE BASED ON ERROR
C                      SIZE

C***********************************************************************

      DIMENSION X (MAXNP), Y (MAXNP)

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD), NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL GRAPH, REMESH, TEST, ERR, SGRAPH, MOVED

      IF (GRAPH) THEN
         SGRAPH = .TRUE.
      ELSE
         SGRAPH = .FALSE.
      ENDIF

      ITERAT = 100
      EPS = .01

      DELANG = 0.
      ANGNOW = 0.
      INTNOW = 0
      IF (GRAPH) THEN
         CALL SYMBOL (1, X(1), Y(1), 'DIAMND')
         CALL PLTFLU
      ENDIF
  100 CONTINUE
      INTNOW = INTNOW + 1
      IF ( ((TANG .GT. 0.) .AND.
     &   (THETA1 + ANGNOW + (DELANG * 1.3) .GT. TANG))
     &   .OR. ((TANG .LT. 0.) .AND.
     &   (THETA1 + ANGNOW + (DELANG * 1.3) .LT. TANG))
     &   ) THEN

C  THE END OF THE LINE (OR CLOSE ENOUGH) HAS BEEN REACHED

C  IF WE ARE TESTING OR THE INTERVALS MATCH, THEN SIMPLY FINISH THE
C  LINE.

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

C  OTHERWISE, MAKE SURE THE INTERVALS ARE ALRIGHT AND ADD THE EXTRA ONE

            EPS = .001
            IF (INTNOW + 1 .NE. NINT) THEN
               CALL MESAGE ('** PROBLEMS WITH INTNOW '//
     &            'IN PLINE **')
               ERR = .TRUE.
               GOTO 160
            ENDIF
            ANG = THETA1 + ANGNOW +
     &         ((TANG - (THETA1 + ANGNOW)) * .5)
            RADIUS = BB * EXP (AA * ANG)
            X (INTNOW + 1) = XCEN + COS (ANG) * RADIUS
            Y (INTNOW + 1) = YCEN + SIN (ANG) * RADIUS
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

C  NOT TO THE END YET

         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG,
     &      BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &      REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &      EMIN, X(INTNOW), Y(INTNOW), S1)
         DELANG = S1 / (BB * EXP (AA * (THETA1 + ANGNOW)))
         IF (TANG .LT. 0.) DELANG = - DELANG
         ANG1 = ANGNOW + DELANG
         RAD1 = BB * EXP (AA * (THETA1 + ANG1))
         XNEW1 = XCEN + COS (THETA1 + ANG1) * RAD1
         YNEW1 = YCEN + SIN (THETA1 + ANG1) * RAD1
         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG,
     &      BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &      REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &      EMIN, XNEW1, YNEW1, S2)
         DELANG = ((S1 + S2) * .5) / (BB * EXP (AA * (THETA1 + ANGNOW)))
         IF (TANG .LT. 0.) DELANG = - DELANG
         ANG1 = ANGNOW + DELANG
         RAD1 = BB * EXP (AA * (THETA1 + ANG1))
         XNEW1 = XCEN + COS (THETA1 + ANG1) * RAD1
         YNEW1 = YCEN + SIN (THETA1 + ANG1) * RAD1
         CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG,
     &      BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &      REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &      EMIN, XNEW1, YNEW1, S3)
         SIZNOW = (((S1 + S2) * .5) + S3) * .5
         DELANG = SIZNOW / (BB * EXP (AA * (THETA1 + ANGNOW)))
         IF (TANG .LT. 0.) DELANG = - DELANG
         ANGNOW = ANGNOW + DELANG
         RADIUS = BB * EXP (AA * (THETA1 + ANGNOW))
         X (INTNOW + 1) = XCEN + COS (THETA1 + ANGNOW) * RADIUS
         Y (INTNOW + 1) = YCEN + SIN (THETA1 + ANGNOW) * RADIUS
         IF (GRAPH) THEN
            CALL SYMBOL (1, X(INTNOW + 1), Y(INTNOW + 1),
     &         'DIAMND')
            CALL MPD2VC (1, X(INTNOW), Y(INTNOW),
     &         X(INTNOW+1), Y(INTNOW+1))
            CALL PLTFLU
         ENDIF
         GOTO 100
      ENDIF

C  ERASE THE LINES FOR SMOOTHING IF NEEDED

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

C  NOW SMOOTH THE NODES ALONG THE LINE

      DO 130 I = 1, ITERAT
         MOVED = .FALSE.
         ANGNOW = 0.
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
            SIZNOW = DTOTAL * TRATIO
            DELANG = SIZNOW / (BB * EXP (AA * (THETA1 + ANGNOW)))
            IF (TANG .LT. 0.) DELANG = -DELANG
            ANGNOW = ANGNOW + DELANG
            ANG = THETA1 + ANGNOW
            RADIUS = BB * EXP (AA * ANG)
            X (J) = XCEN + COS (ANG) * RADIUS
            Y (J) = YCEN + SIN (ANG) * RADIUS
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

  160 CONTINUE

      RETURN

      END
