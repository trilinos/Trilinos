C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TUCK (MXND, MLN, NUID, XN, YN, LXK, KXL, NXL, LXN,
     &   LNODES, IAVAIL, NAVAIL, LLL, KKK, NNN, N1, NLOOP, GRAPH,
     &   NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE TUCK = COLLAPSES TWO SIDE LINES INTO A ROW END NODE.
C                      THIS IS REFERRED TO AS A TUCK.

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20)

      LOGICAL GRAPH, ERR, NOROOM

      ERR = .FALSE.

C  CHECK TO MAKE SURE THAT THE NODE STILL EXISTS

      IF (LXN (1, N1) .LE. 0) THEN
         ERR = .TRUE.
         CALL MESAGE ('** PROBLEMS IN TUCK - N1 DOES NOT EXIST **')
         GOTO 290
      ENDIF

C  GET ALL THE DEFINITIONS IN ORDER

      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      L1 = LNODES (5, N0)
      L2 = LNODES (5, N1)
      KOLD = KXL (1, L1)
      KL2 = KXL (1, L2)

C  FIND L5 AND NC2

      DO 100 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF (LTEST .NE. L1) THEN
            IF (NXL (1, LTEST) .EQ. N1) THEN
               L5 = LTEST
               NC2 = NXL (2, LTEST)
               GOTO 110
            ELSEIF (NXL (2, LTEST) .EQ. N1) THEN
               L5 = LTEST
               NC2 = NXL (1, LTEST)
               GOTO 110
            ENDIF
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L5 AND NC2 **')
      ERR = .TRUE.
      GOTO 290
  110 CONTINUE

C  FIND L4 AND NC1

      DO 120 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF ( (LTEST .NE. L1) .AND. (LTEST .NE. L5) ) THEN
            IF (NXL (1, LTEST) .EQ. N0) THEN
               L4 = LTEST
               NC1 = NXL (2, LTEST)
               GOTO 130
            ELSEIF (NXL (2, LTEST) .EQ. N0) THEN
               L4 = LTEST
               NC1 = NXL (1, LTEST)
               GOTO 130
            ENDIF
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L4 AND NC1 **')
      ERR = .TRUE.
      GOTO 290
  130 CONTINUE

C  FIND L3

      DO 140 I = 1, 4
         LTEST = LXK (I, KOLD)
         IF ( (LTEST .NE. L1) .AND. (LTEST .NE. L5) .AND.
     &      (LTEST .NE. L4) ) THEN
            L3 = LTEST
            GOTO 150
         ENDIF
  140 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L3 **')
      ERR = .TRUE.
      GOTO 290
  150 CONTINUE

C  FIND THE ELEMENT KL5

      IF (KXL (1, L5) .EQ. KOLD) THEN
         KL5 = KXL (2, L5)
      ELSEIF (KXL (2, L5) .EQ. KOLD) THEN
         KL5 = KXL (1, L5)
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK FINDING KL5 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF

C  NOW THAT ALL THE NECESSARY VARIABLES HAVE BEEN DEFINED,
C  START BY DELETING LINE L1, L2, AND L5

      KXL (1, L1) = 0
      KXL (2, L1) = 0
      NXL (1, L1) = 0
      NXL (2, L1) = 0
      KXL (1, L2) = 0
      KXL (2, L2) = 0
      NXL (1, L2) = 0
      NXL (2, L2) = 0
      KXL (1, L5) = 0
      KXL (2, L5) = 0
      NXL (1, L5) = 0
      NXL (2, L5) = 0

C  NOW FIX THE KXL ARRAY FOR LINE L3 HAVING KL5 INSTEAD OF KOLD

      IF (KXL (1, L3) .EQ. KOLD) THEN
         KXL (1, L3) = KL5
      ELSEIF (KXL (2, L3) .EQ. KOLD) THEN
         KXL (2, L3) = KL5
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK REPLACING KOLD FOR L3 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF

C  NOW FIX THE KXL ARRAY FOR LINE L3 HAVING KL5 INSTEAD OF KOLD

      IF (KXL (1, L4) .EQ. KOLD) THEN
         KXL (1, L4) = KL2
      ELSEIF (KXL (2, L4) .EQ. KOLD) THEN
         KXL (2, L4) = KL2
      ELSE
         CALL MESAGE ('** PROBLEMS IN TUCK REPLACING KOLD FOR L4 **')
         ERR = .TRUE.
         GOTO 290
      ENDIF

C  FIX THE LINES PER ELEMENT ARRAY FOR ELEMENT KL5 TO REFLECT
C  L3 TAKING L5'S PLACE

      DO 160 I = 1, 4
         IF (LXK (I, KL5) .EQ. L5) THEN
            LXK (I, KL5) = L3
            GOTO 170
         ENDIF
  160 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L5 IN KL5 **')
      ERR = .TRUE.
      GOTO 290
  170 CONTINUE

C  FIX THE LINES PER ELEMENT ARRAY FOR ELEMENT KL2 TO REFLECT
C  L4 TAKING L2'S PLACE

      DO 180 I = 1, 4
         IF (LXK (I, KL2) .EQ. L2) THEN
            LXK (I, KL2) = L4
            GOTO 190
         ENDIF
  180 CONTINUE
      CALL MESAGE ('** PROBLEMS IN TUCK FINDING L2 IN KL2 **')
      ERR = .TRUE.
      GOTO 290
  190 CONTINUE

C  RECONNECT ALL LINES CONNECTED TO N1 TO NC1 EXCEPT L5 AND L2

      CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK GETTING N1 LINES **')
         GOTO 290
      ENDIF
      IF (GRAPH) CALL LCOLOR ('BLACK')
      DO 200 I = 1, NL
         LL = L1LIST (I)
         IF ((GRAPH) .AND. (NXL (1, LL) .GT. 0) .AND.
     &      (NXL (2, LL) .GT. 0) )
     &      CALL D2NODE (MXND, XN, YN, NXL (1, LL), NXL (2, LL))
         IF (NXL (1, LL) .EQ. N1) THEN
            NXL (1, LL) = NC1
         ELSEIF (NXL (2, LL) .EQ. N1) THEN
            NXL (2, LL) = NC1
         ENDIF
  200 CONTINUE
      IF (GRAPH) THEN
         CALL LCOLOR ('WHITE')
         CALL SFLUSH
      ENDIF

C  FIX LXN ARRAY
C  UNHOOK L1, L2 AND L5 FROM N1

      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L1, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N1 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L2, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L2 FROM N1 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N1,
     &   L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L5 FROM N1 **')
         GOTO 290
      ENDIF

C  ADD ALL LINES STILL HOOKED TO N1 TO THE LIST OF LINES FOR NC1

      DO 210 I = 1, NL
         LL = L1LIST (I)
         IF ((LL .NE. L2) .AND. (LL .NE. L5) .AND. (LL .NE. L1)) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         NC1, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN TUCK HOOKING N1'' LINES'//
     &            ' TO NC1 **')
               GOTO 290
            ENDIF
         ENDIF
  210 CONTINUE

C  DELETE N1

      DO 220 I = 1, 3
         LXN (I, N1) = 0
  220 CONTINUE
      LXN (4, N1) = IAVAIL
      IAVAIL = N1
      NAVAIL = NAVAIL+1
      NUID (N1) = 0

C  RECONNECT ALL LINES CONNECTED TO N2 TO N0 (EXCEPT L2)

      CALL GETLXN (MXND, LXN, N2, L1LIST, NL, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK GETTING N2 LINES **')
         GOTO 290
      ENDIF
      IF (GRAPH) CALL LCOLOR ('BLACK')
      DO 230 I = 1, NL
         LL = L1LIST (I)
         IF ((GRAPH) .AND. (NXL (1, LL) .GT. 0) .AND.
     &      (NXL (2, LL) .GT. 0) )
     &      CALL D2NODE (MXND, XN, YN, NXL (1, LL), NXL (2, LL))
         IF (NXL (1, LL) .EQ. N2) THEN
            NXL (1, LL) = N0
         ELSEIF (NXL (2, LL) .EQ. N2) THEN
            NXL (2, LL) = N0
         ENDIF
  230 CONTINUE
      IF (GRAPH) THEN
         CALL LCOLOR ('WHITE')
         CALL SFLUSH
      ENDIF

C  FIX LXN ARRAY
C  UNHOOK L2 FROM N2, L1 FROM N0, AND L5 FROM NC2

      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N2,
     &   L2, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L2 FROM N2 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N0,
     &   L1, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N0 **')
         GOTO 290
      ENDIF
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NC2,
     &   L5, NNN, ERR, NOROOM)
      IF ((NOROOM) .OR. (ERR)) THEN
         CALL MESAGE ('** PROBLEMS IN TUCK UNHOOKING L1 FROM N0 **')
         GOTO 290
      ENDIF

C  ADD ALL LINES STILL HOOKED TO N2 TO THE LIST OF LINES FOR N0

      DO 240 I = 1, NL
         LL = L1LIST (I)
         IF (LL .NE. L2) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL,
     &         N0, LL, NNN, ERR, NOROOM)
            IF ((NOROOM) .OR. (ERR)) THEN
               CALL MESAGE ('** PROBLEMS IN TUCK HOOKING N2'' LINES'//
     &            ' TO N0 **')
               GOTO 290
            ENDIF
         ENDIF
  240 CONTINUE

C  DELETE N2

      DO 250 I = 1, 3
         LXN (I, N2) = 0
  250 CONTINUE
      LXN (4, N2) = IAVAIL
      IAVAIL = N2
      NAVAIL = NAVAIL+1
      NUID (N2) = 0

C  NOW DELETE THE OLD ELEMENT

      DO 260 I = 1, 4
         LXK (I, KOLD) = 0
  260 CONTINUE

C  NOW FIX THE LNODES ARRAY

      LNODES (3, N0) = LNODES (3, N2)
      LNODES (2, LNODES (3, N2) ) = N0
      LNODES (5, N0) = LNODES (5, N2)

      NLOOP = NLOOP - 2
      ERR = .FALSE.

C  NOW REDRAW THE ELEMENTS

      IF (GRAPH) THEN
         CALL LCOLOR ('BLACK')
         CALL D2NODE (MXND, XN, YN, N0, N1)
         CALL D2NODE (MXND, XN, YN, NC2, N1)
         CALL D2NODE (MXND, XN, YN, N2, N1)
         CALL LCOLOR ('WHITE')
         CALL GETLXN (MXND, LXN, N0, L1LIST, NL, ERR)
         IF (ERR) GOTO 290
         DO 270 II = 1, NL
            IDRAW = L1LIST (II)
            NODE1 = NXL (1, IDRAW)
            NODE2 = NXL (2, IDRAW)
            CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  270    CONTINUE
         CALL GETLXN (MXND, LXN, NC1, L1LIST, NL, ERR)
         IF (ERR) GOTO 290
         DO 280 II = 1, NL
            IDRAW = L1LIST (II)
            NODE1 = NXL (1, IDRAW)
            NODE2 = NXL (2, IDRAW)
            CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  280    CONTINUE
         CALL SFLUSH
      ENDIF

C  FLAG NODES FOR SMOOTHING

      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NC1, ERR)
      IF (ERR) GOTO 290
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NC2, ERR)
      IF (ERR) GOTO 290

  290 CONTINUE

      RETURN

      END
