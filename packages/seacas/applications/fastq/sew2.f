C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SEW2 (MXND, MLN, NUID, LXK, KXL, NXL, LXN, LNODES,
     &   IAVAIL, NAVAIL, LLL, KKK, NNN, I1, I2, J1, J2, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE SEW2 = COLLAPSES A LOOP INTO TWO POSSIBLE LOOPS

C***********************************************************************

      DIMENSION NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20)

      LOGICAL ERR, NOROOM

      ERR = .FALSE.

C  GET THE APPROPRIATE LINES AND NODES TO BE DELETED

      IF ((LXN (2, J1) .LT. 0) .AND. (LXN (2, J2) .LT. 0)) THEN
         LSTAY = LNODES (5, J1)
         LGONE = LNODES (5, I1)
         NGONE1 = I1
         NGONE2 = I2
         NSTAY1 = J2
         NSTAY2 = J1

      ELSEIF (LXN (2, J1) .LT. 0) THEN
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = J2
         NGONE2 = I2
         NSTAY1 = I1
         NSTAY2 = J1

      ELSEIF (LXN (2, J2) .LT. 0) THEN
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = I1
         NGONE2 = J1
         NSTAY1 = J2
         NSTAY2 = I2

      ELSE
         LSTAY = LNODES (5, I1)
         LGONE = LNODES (5, J1)
         NGONE1 = J2
         NGONE2 = J1
         NSTAY1 = I1
         NSTAY2 = I2
      ENDIF

      KOLD = KXL (1, LGONE)
      KNEW = KXL (1, LSTAY)

C  DELETE THE OLD LINE AND REDO LINK ARRAYS

      IF (KNEW .EQ. 0) THEN
         KXL (1, LSTAY) = KOLD
         KXL (2, LSTAY) = 0
      ELSE
         KXL (1, LSTAY) = KNEW
         KXL (2, LSTAY) = KOLD
      ENDIF

      KXL (1, LGONE) = 0
      KXL (2, LGONE) = 0
      NXL (1, LGONE) = 0
      NXL (2, LGONE) = 0

C  FIX THE LINES PER ELEMENT ARRAY FOR THE ONE ELEMENT CHANGING

      IF (KOLD .GT. 0) THEN
         DO 100 II = 1, 4
            IF (LXK (II, KOLD) .EQ. LGONE) THEN
               LXK (II, KOLD) = LSTAY
               GOTO 110
            ENDIF
  100    CONTINUE

         CALL MESAGE ('** PROBLEMS IN SEW2 FIXING THE CHANGING'//
     &      'ELEMENT **')
         ERR = .TRUE.
         GOTO 180

  110    CONTINUE
      ENDIF

C  RECONNECT ALL LINES CONNECTING TO NGONE2 TO NSTAY2

      CALL GETLXN (MXND, LXN, NGONE2, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN SEW2 FINDING LXN FOR NGONE2 **')
         GOTO 180
      ENDIF
      DO 120 II = 1, NL
         LL = L1LIST (II)
         IF (NXL (1, LL) .EQ. NGONE2) THEN
            NXL (1, LL) = NSTAY2
         ELSEIF (NXL (2, LL) .EQ. NGONE2) THEN
            NXL (2, LL) = NSTAY2
         ENDIF
  120 CONTINUE

C  FIX LXN ARRAY
C  UNHOOK LGONE FROM NGONE2 OR NSTAY2 AS NEEDED

      IF (LGONE .EQ. LNODES (5, J1)) THEN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, J1,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE2 LINES **')
            GOTO 180
         ENDIF

      ELSE
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, I2,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE2 LINES **')
            GOTO 180
         ENDIF

      ENDIF

C  ADD ALL LINES STILL HOOKED TO NGONE2 TO THE LIST OF LINES FOR NSTAY2

      DO 130 II = 1, NL
         LL = L1LIST (II)
         IF (LL .NE. LGONE) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NSTAY2, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 ADDING NSTAY2 '//
     &            'LINES **')
               GOTO 180
            ENDIF
         ENDIF
  130 CONTINUE

C  DELETE NGONE2 (UNHOOK EVERYTHING FROM IT)

      DO 140 II = 1, 3
         LXN (II, NGONE2) = 0
  140 CONTINUE
      LXN (4, NGONE2) = IAVAIL
      IAVAIL = NGONE2
      NAVAIL = NAVAIL+1
      NUID (NGONE2) = 0

C  RECONNECT ALL LINES CONNECTING TO NGONE1 TO NSTAY1

      CALL GETLXN (MXND, LXN, NGONE1, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN SEW2 GETTING NGONE1 LINES **')
         GOTO 180
      ENDIF
      DO 150 II = 1, NL
         LL = L1LIST (II)
         IF (NXL (1, LL) .EQ. NGONE1) THEN
            NXL (1, LL) = NSTAY1
         ELSEIF (NXL (2, LL) .EQ. NGONE1) THEN
            NXL (2, LL) = NSTAY1
         ENDIF
  150 CONTINUE

C  FIX LXN ARRAY
C  UNHOOK LGONE FROM NGONE1 OR NSTAY1 AS APPROPRIATE

      IF (LGONE .EQ. LNODES (5, I1)) THEN
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, I1,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE1 LINES **')
            GOTO 180
         ENDIF

      ELSE
         CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, J2,
     &      LGONE, NNN, ERR, NOROOM)
         IF ((NOROOM) .OR. (ERR)) THEN
            CALL MESAGE ('** PROBLEMS IN SEW2 DELETING NGONE1 LINES **')
            GOTO 180
         ENDIF

      ENDIF

C  ADD ALL LINES STILL HOOKED TO NGONE1 TO THE LIST OF LINES FOR NSTAY1

      DO 160 II = 1, NL
         LL = L1LIST (II)
         IF (LL .NE. LGONE) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NSTAY1, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN SEW2 ADDING NSTAY1'//
     &            ' LINES **')
               GOTO 180
            ENDIF
         ENDIF
  160 CONTINUE

C  DELETE NGONE1 (UNHOOK EVERYTHING FROM IT)

      DO 170 II = 1, 3
         LXN (II, NGONE1) = 0
  170 CONTINUE
      LXN (4, NGONE1) = IAVAIL
      IAVAIL = NGONE1
      NAVAIL = NAVAIL+1
      NUID (NGONE1) = 0

C  NOW FIX THE LNODES ARRAY FOR BOTH OF THE LOOPS

      IF (NGONE2 .EQ. J1) THEN
         LNODES (2, NSTAY2) = LNODES (2, NGONE2)
         LNODES (3, LNODES (2, NGONE2)) = NSTAY2
      ELSE
         LNODES (2, LNODES (3, NGONE2)) = NSTAY2
         LNODES (3, NSTAY2) = LNODES (3, NGONE2)
         LNODES (5, NSTAY2) = LNODES (5, NGONE2)
      ENDIF
      IF (NGONE1 .EQ. J2) THEN
         LNODES (2, LNODES (3, NGONE1)) = NSTAY1
         LNODES (3, NSTAY1) = LNODES (3, NGONE1)
         LNODES (5, NSTAY1) = LNODES (5, NGONE1)
      ELSE
         LNODES (2, NSTAY1) = LNODES (2, NGONE1)
         LNODES (3, LNODES (2, NGONE1)) = NSTAY1
      ENDIF

      I1 = NSTAY1
      I2 = NSTAY2
      J1 = NGONE1
      J2 = NGONE2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   I1, ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, I1), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, LNODES (3, I1)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, I1), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, LNODES (2, I1)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   I2, ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, I2), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (3, LNODES (3, I2)), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, I2), ERR)
      IF (ERR) GOTO 180
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   LNODES (2, LNODES (2, I2)), ERR)
      IF (ERR) GOTO 180

  180 CONTINUE

      RETURN

      END
