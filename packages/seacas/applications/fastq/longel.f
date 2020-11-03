C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LONGEL (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, NAVAIL, IAVAIL, NODE, KELEM, ANG, TOLER,
     &   N1, N2, KREG, XMIN, XMAX, YMIN, YMAX, KKK, LLL, DONE, GRAPH,
     &   VIDEO, NOROOM, ERR, KKKADD)
C***********************************************************************

C  SUBROUTINE LONGEL = AN ELONGATED ELEMENT OVER 150 DEGREES GETS A
C                      3 ELEMENT REPLACEMENT FOR THE TWO ELEMENTS THERE

C***********************************************************************

      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), NUID(MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION NODES(4)

      LOGICAL NOROOM, ERR, DONE, GRAPH, CCW, VIDEO

      CCW = .TRUE.

C  SEE IF THE ANGLE IS WITHIN BOUNDS

      IF (ANG .GT. TOLER) THEN
         CALL GNXKA (MXND, XN, YN, KELEM, NODES, AREA, LXK, NXL, CCW)
         NODE2 = NODES(1) + NODES(2) + NODES(3) + NODES(4) - NODE
     &      - N1 - N2
         D12 = SQRT ( ((XN (NODE2) - XN (N1)) ** 2) +
     &      ((YN (NODE2) - YN (N1)) ** 2) )
         D22 = SQRT ( ((XN (NODE2) - XN (N2)) ** 2) +
     &      ((YN (NODE2) - YN (N2)) ** 2) )
         DN1 = SQRT ( ((XN (NODE) - XN (N1)) ** 2) +
     &      ((YN (NODE) - YN (N1)) ** 2) )
         DN2 = SQRT ( ((XN (NODE) - XN (N2)) ** 2) +
     &      ((YN (NODE) - YN (N2)) ** 2) )
         DMIN = ((DN1 + DN2) * .5) * 1.7
         DMAX = (DN1 + DN2) * .5

C  SEE IF IT IS A LONG LEGGED BEAST

         IF ((D12 .GT. DMIN) .OR. (D22 .GT. DMIN)) THEN

C  FIND L1, L2, L3, AND L4

            DO 100 I = 1, 4
               LTEST = LXK (I, KELEM)
               IF ( ((NXL (1, LTEST) .EQ. NODE) .AND.
     &            (NXL (2, LTEST) .EQ. N1)) .OR.
     &            ((NXL (2, LTEST) .EQ. NODE) .AND.
     &            (NXL (1, LTEST) .EQ. N1)) ) THEN
                  L1 = LTEST
                  GOTO 110
               ENDIF
  100       CONTINUE
            CALL MESAGE('** PROBLEMS IN LONGEL FINDING L1 **')
            ERR = .TRUE.
            GOTO 250
  110       CONTINUE

            DO 120 I = 1, 4
               LTEST = LXK (I, KELEM)
               IF ( ((NXL (1, LTEST) .EQ. NODE) .AND.
     &            (NXL (2, LTEST) .EQ. N2)) .OR.
     &            ((NXL (2, LTEST) .EQ. NODE) .AND.
     &            (NXL (1, LTEST) .EQ. N2)) ) THEN
                  L2 = LTEST
                  GOTO 130
               ENDIF
  120       CONTINUE
            CALL MESAGE('** PROBLEMS IN LONGEL FINDING L2 **')
            ERR = .TRUE.
            GOTO 250
  130       CONTINUE

            DO 140 I = 1, 4
               LTEST = LXK (I, KELEM)
               IF ( ((NXL (1, LTEST) .EQ. NODE2) .AND.
     &            (NXL (2, LTEST) .EQ. N1)) .OR.
     &            ((NXL (2, LTEST) .EQ. NODE2) .AND.
     &            (NXL (1, LTEST) .EQ. N1)) ) THEN
                  IF (D12 .GT. D22) THEN
                     L4 = LTEST
                  ELSE
                     L3 = LTEST
                  ENDIF
                  GOTO 150
               ENDIF
  140       CONTINUE
            CALL MESAGE('** PROBLEMS IN LONGEL FINDING L4/L3 **')
            ERR = .TRUE.
            GOTO 250
  150       CONTINUE

            DO 160 I = 1, 4
               LTEST = LXK (I, KELEM)
               IF ( ((NXL (1, LTEST) .EQ. NODE2) .AND.
     &            (NXL (2, LTEST) .EQ. N2)) .OR.
     &            ((NXL (2, LTEST) .EQ. NODE2) .AND.
     &            (NXL (1, LTEST) .EQ. N2)) ) THEN
                  IF (D12 .GT. D22) THEN
                     L3 = LTEST
                  ELSE
                     L4 = LTEST
                  ENDIF
                  GOTO 170
               ENDIF
  160       CONTINUE
            CALL MESAGE('** PROBLEMS IN LONGEL FINDING L3/L4 **')
            ERR = .TRUE.
            GOTO 250
  170       CONTINUE

C  NOW FIND KELEM2

            KELEM2 = KXL (1, L4) + KXL (2, L4) - KELEM
            IF (KELEM2 .EQ. 0) GOTO 250

C  NOW FIND NODE3 - THE NODE THAT WILL BE PART OF THE NEWLY
C  FORMED ELEMENT

            CALL GNXKA (MXND, XN, YN, KELEM2, NODES, AREA,
     &         LXK, NXL, CCW)

            IF (D12 .GT. D22) THEN
               DO 180 I = 1, 4
                  IF (NODES (I) .EQ. N1) THEN
                     IF (I .EQ. 1) THEN
                        NODE3 = NODES (4)
                     ELSE
                        NODE3 = NODES (I - 1)
                     ENDIF
                     GOTO 190
                  ENDIF
  180          CONTINUE
               CALL MESAGE('** PROBLEMS IN LONGEL FINDING NODE3/D11 **')
               ERR = .TRUE.
               GOTO 250
  190          CONTINUE

            ELSE
               DO 200 I = 1, 4
                  IF (NODES (I) .EQ. N2) THEN
                     IF (I .EQ. 4) THEN
                        NODE3 = NODES (1)
                     ELSE
                        NODE3 = NODES (I + 1)
                     ENDIF
                     GOTO 210
                  ENDIF
  200          CONTINUE
               CALL MESAGE('** PROBLEMS IN LONGEL FINDING NODE3/D11 **')
               ERR = .TRUE.
               GOTO 250
  210          CONTINUE
            ENDIF

C  NOW FIND L5

            DO 220 I = 1, 4
               LTEST = LXK (I, KELEM2)
               IF (D12 .GT. D22) THEN
                  IF ( ( (NXL (1, LTEST) .EQ. NODE3) .AND.
     &               (NXL (2, LTEST) .EQ. N1) ) .OR.
     &               ( (NXL (1, LTEST) .EQ. N1) .AND.
     &               (NXL (2, LTEST) .EQ. NODE3) ) ) THEN
                     L5 = LTEST
                     GOTO 230
                  ENDIF
               ELSE
                  IF ( ( (NXL (1, LTEST) .EQ. NODE3) .AND.
     &               (NXL (2, LTEST) .EQ. N2) ) .OR.
     &               ( (NXL (1, LTEST) .EQ. N2) .AND.
     &               (NXL (2, LTEST) .EQ. NODE3) ) ) THEN
                     L5 = LTEST
                     GOTO 230
                  ENDIF
               ENDIF
  220       CONTINUE
            CALL MESAGE('** PROBLEMS IN LONGEL FINDING L5 **')
            ERR = .TRUE.
            GOTO 250
  230       CONTINUE

C  NOW CHECK TO SEE IF IT MAKES SENSE TO ADD THE EXTRA ELEMENT TO
C  IMPROVE AN ELONGATED ONE

            DL3 = SQRT ( ((XN (NXL (1, L3)) - XN (NXL (2, L3))) ** 2) +
     &         ((YN (NXL (1, L3)) - YN (NXL (2, L3))) ** 2) )
            DL4 = SQRT ( ((XN (NXL (1, L4)) - XN (NXL (2, L4))) ** 2) +
     &         ((YN (NXL (1, L4)) - YN (NXL (2, L4))) ** 2) )
            DL5 = SQRT ( ((XN (NXL (1, L5)) - XN (NXL (2, L5))) ** 2) +
     &         ((YN (NXL (1, L5)) - YN (NXL (2, L5))) ** 2) )
            IF ((DL3 .GT. DMAX) .OR. (DL5 .GT. DL4)) GOTO 250

C  ADD THE EXTRA ELEMENT

            IF (GRAPH) THEN
               CALL LCOLOR ('PINK ')
               DO 240 IL = 1, 4
                  IL1 = NXL (1, LXK (IL, KELEM))
                  IL2 = NXL (2, LXK (IL, KELEM))
                  CALL D2NODE (MXND, XN, YN, IL1, IL2)
                  IL1 = NXL (1, LXK (IL, KELEM2))
                  IL2 = NXL (2, LXK (IL, KELEM2))
                  CALL D2NODE (MXND, XN, YN, IL1, IL2)
  240          CONTINUE
               CALL LCOLOR ('WHITE')
               CALL SFLUSH
            ENDIF
            IF ((D12 .GT. D22) .AND. (LXN (4, N1) .NE. 0)) THEN
               CALL UNDELM (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL,
     &            NXL, LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, NODE2, NODE3,
     &            N1, NODE, L5, L1, L4, KELEM, KELEM2, NOROOM, ERR,
     &            GRAPH, VIDEO)
               IF ((ERR) .OR. (DONE)) GOTO 250
               DONE = .TRUE.
            ELSEIF (LXN (4, N2) .NE. 0) THEN
               CALL UNDELM (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL,
     &            NXL, LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, NODE2, NODE,
     &            N2, NODE3, L2, L5, L4, KELEM2, KELEM, NOROOM, ERR,
     &            GRAPH, VIDEO)
               IF ((ERR) .OR. (DONE)) GOTO 250
               DONE = .TRUE.
            ENDIF
            KKKADD = KKK
         ENDIF
      ENDIF

  250 CONTINUE
      RETURN

      END
