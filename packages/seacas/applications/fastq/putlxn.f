C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PUTLXN (MXND, NL, LXN, NUID, NODE, LINES, NAVAIL,
     &   IAVAIL, NNN, ERR, NOROOM)
C***********************************************************************

C  SUBROUTINE PUTLXN = DEFINE THE LINES FOR NODE AS  (LINES (I), I=1, NL)

C***********************************************************************

C  NOTE:
C     SAME CONTINUATION ENTRIES ARE USED AS ALREADY IN USE
C     FOR THIS NODE.
C     THIS NODE WILL BE FLAGGED AS A BOUNDARY NODE IF
C     LXN (2, NODE) .LT. 0   (NOT IF LINES (2) .LT. 0)

C***********************************************************************

      DIMENSION LINES (NL), LXN (4, MXND), NUID (MXND)

      LOGICAL ERR, BOUND, NOROOM

      BOUND = .FALSE.

      IF (LXN (2, NODE) .LT. 0)BOUND = .TRUE.
      NN = NODE
      NDONE = 0

C  FILL IN NEXT 3  (4 LAST TIME) NODES

  100 CONTINUE
      N4 = LXN (4, NN)
      NR = MIN0 (4, NL - NDONE)
      DO 110 I = 1, NR
         J = NDONE + I
         LXN (I, NN) = LINES (J)
  110 CONTINUE

C  CLEAR REMAINING PORTION

      IF (NR .LT. 4) THEN
         NZ = NR + 1
         DO 120 I = NZ, 4
            LXN (I, NN) = 0
  120    CONTINUE
      ENDIF

C  TAG BOUNDARY NODES

      IF (BOUND)LXN (2, NN) =  - LXN (2, NN)

C  TAG CONTINUATIONS

      IF (NDONE .GT. 1)LXN (1, NN) =  - LXN (1, NN)
      IF (NDONE + 4 .GE. NL) THEN

C  COLLECT GARBAGE

  130    CONTINUE
         IF (N4 .GE. 0) THEN

C  UPDATE NNN

  140       CONTINUE
            IF (LXN (1, NNN) .NE. 0) THEN
               RETURN
            ELSE
               NNN = NNN - 1
               GOTO 140
            ENDIF
         ENDIF

         NR =  - N4
         N4 = LXN (4, NR)
         LXN (1, NR) = 0
         LXN (2, NR) = 0
         LXN (3, NR) = 0
         LXN (4, NR) = IAVAIL
         IAVAIL = NR
         NAVAIL = NAVAIL + 1
         GOTO 130
      ENDIF

C  NEED ANOTHER LINE IN THE TABLE

      NDONE = NDONE + 3
      NEXTR = IABS (N4)
      IF (N4 .LT. 0) THEN
         LXN (4, NN) =  - NEXTR
         NN = NEXTR
         GOTO 100
      ENDIF

C  RESERVE A NEW LINE IN LXN TABLE

      IF (NAVAIL .LT. 1) THEN
         WRITE ( * , 10000)NODE
         ERR = .TRUE.
         NOROOM = .TRUE.
         RETURN
      ENDIF
      NEW = IAVAIL
      IF (NEW .GT. NNN)NNN = NEW
      IAVAIL = LXN (4, IAVAIL)
      NAVAIL = NAVAIL - 1
      LXN (4, NEW) = 0
      NUID (NEW) = 0
      NEXTR = NEW

C  INSERT LINK TO NEXT LINE

      LXN (4, NN) =  - NEXTR
      NN = NEXTR
      GOTO 100

10000 FORMAT (' NODE TABLE OVERFLOW IN PUTLXN - NODE = ', I5)

      END
