C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE NCKLCE(MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &   NNN, NNNOLD, LLL, NAVAIL, IAVAIL, EPS, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE NCKLCE = INSERTS AN EXRTRA RING OF ELEMENTS JUST INSIDE
C                      THE REGION BOUNDARY

C***********************************************************************

C  NOTE:
C     ONLY ARRAYS LXK, NXL, XN, AND YN ARE INPUT TO NCKLCE.
C     ARRAYS KXL AND LXN ARE RECREATED BY SUBROUTINE CROSS AFTER
C     LXK AND NXL ARE MODIFIED.

C***********************************************************************

      DIMENSION LINES(20)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), NXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), NUID(MXND)

      LOGICAL ERR, NOROOM

      NNNX = NNN
      LLLX = LLL
      KKKX = KKK
      ERR = .TRUE.
      NOROOM = .TRUE.

C  COUNT BOUNDARY NODES TO CHECK FOR IMPENDING OVERFLOW

      NUMB = 0
      DO 100 I = NNNOLD + 1, NNN
         IF ((LXN(2, I) .LT. 0) .AND. (LXN(1, I) .GT. 0)) THEN
            NUMB = NUMB + 1
         END IF
  100 CONTINUE
      IF (NNN + NUMB .GT. MXND) THEN
         CALL MESAGE ('NODE TABLE OVERFLOW IN NCKLCE')
         RETURN
      ELSE IF (LLL + 2*NUMB .GT. MXND*3) THEN
         CALL MESAGE ('LINE TABLE OVERFLOW IN NCKLCE')
         RETURN
      ELSE IF (KKK + NUMB .GT. MXND) THEN
         CALL MESAGE ('ELEMENT TABLE OVERFLOW IN NCKLCE')
         RETURN
      END IF
      NOROOM = .FALSE.

C     FIND FIRST BOUNDARY NODE

      DO 110 I = NNNOLD + 1, NNN
         IF (LXN(1, I) .GT. 0) THEN
            NBDY1 = I
            IF (LXN(2, I) .LT. 0) GO TO 120
         END IF
  110 CONTINUE
      NODE = -1
      WRITE (*, 10000) NODE
      RETURN
  120 CONTINUE
      NOLD = -1
      NODE = NBDY1

C  FIND NEXT NODE ON THE BOUNDARY
C  LOOK AT ALL NEIGHBORING NODES

  130 CONTINUE
      CALL GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
      IF (ERR) RETURN

      DO 140 IL = 1, NL
         L = LINES(IL)
         IM = NXL(1, L) + NXL(2, L) - NODE

C  DISALLOW PREVIOUS NODE AND NON-BOUNDARY LINES

         IF ((IM .NE. NOLD) .AND. (KXL(2, L) .LE. 0)) THEN
            NNXT = IM
            GO TO 150
         END IF
  140 CONTINUE
      ERR = .TRUE.
      WRITE (*, 10000) NODE
      RETURN
  150 CONTINUE

C  CREATE A NEW NODE *ON-TOP-OF* OLD BOUNDARY NODE

      NNN = NNN + 1
      XN(NNN) = XN(NODE)
      YN(NNN) = YN(NODE)
      NUID(NNN) = NUID(NODE)
      NUID(NODE) = 0
      DO 160 I = 1, 4
         LXN(I, NNN) = 0
  160 CONTINUE

C  CREATE TWO NEW LINES -- ONE CONNECTING TO THIS NODE, AND
C  ONE ON *TOP* OF THE NEW ELEMENT

      LLL = LLL + 1
      NXL(1, LLL) = NODE
      NXL(2, LLL) = NNN
      KXL(1, LLL) = 0
      KXL(2, LLL) = 0
      LLL = LLL + 1
      NXL(1, LLL) = NNN
      NXL(2, LLL) = NNN + 1
      KXL(1, LLL) = 0
      KXL(2, LLL) = 0
      IF (NNXT .EQ. NBDY1)NXL(2, LLL) = NNNX + 1

C     CREATE A NEW ELEMENT

      KKK = KKK + 1
      LXK(1, KKK) = L
      LXK(2, KKK) = LLL - 1
      LXK(3, KKK) = LLL
      LXK(4, KKK) = LLL + 1
      IF (NNXT .EQ. NBDY1)LXK(4, KKK) = LLLX + 1

C  CHECK FOR COMPLETION OF LOOP AROUND BOUNDARY

      IF (NNXT .NE. NBDY1) THEN
         NOLD = NODE
         NODE = NNXT
         GO TO 130
      END IF

C     RE-SETUP AVAILABLE LXN-SPACE LINKS

      IOLD = 0
      NAVAIL = 0
      DO 170 I = 1, NNNX
         IF (LXN(1, I) .EQ. 0) THEN
            IF (IOLD .LE. 0) THEN
               IAVAIL = I
            ELSE
               LXN(4, IOLD) = I
            END IF
            IOLD = I
            NAVAIL = NAVAIL + 1
         END IF
  170 CONTINUE
      IF (IOLD .LE. 0) THEN
         IAVAIL = NNN + 1
      ELSE
         LXN(4, IOLD) = NNN + 1
      END IF
      NAVAIL = NAVAIL + (MXND - NNN)
      IF (NNN .LT. MXND) THEN
         NNN1 = NNN + 1
         DO 180 I = NNN1, MXND
            LXN(1, I) = 0
            LXN(2, I) = 0
            LXN(3, I) = 0
            LXN(4, I) = I + 1
  180    CONTINUE
      END IF

C     COMPLETE KXL AND LXN ARRAYS

      CALL CCROSS(4, KKK, 2, LLL, LXK, KXL, KKKX + 1, LLLX + 1, NOROOM,
     &   ERR)
      IF (NOROOM)RETURN
      IF (ERR) THEN
         CALL MESAGE ('ERROR IN NCKLCE - LXK TABLE GENERATION')
         RETURN
      END IF

      LLLX1 = LLLX + 1
      DO 200 L = LLLX1, LLL
         DO 190 I = 1, 2
            CALL ADDLXN(MXND, LXN, NUID, NAVAIL, IAVAIL, NXL(I, L), L,
     &         NNN, ERR, NOROOM)
            IF (ERR) THEN
               CALL MESAGE ('ERROR IN NCKLCE - NXL TABLE GENERATION')
               RETURN
            END IF
  190    CONTINUE
  200 CONTINUE

C  USE SMOGS TO REPOSITION THE OLD BOUNDARY NODES

      DO 210 I = NNNOLD + 1, NNN
         LXN(2, I) = -LXN(2, I)
  210 CONTINUE
      RONE = 1.
      CALL SMOGS(MXND, XN, YN, NXL, LXN, NNN, NNNOLD, 3, EPS, RONE)

C  FLAG NEW BOUNDARY NODES (ONLY)

      DO 220 I = NNNOLD + 1, NNN
         LXN(2, I) = IABS(LXN(2, I))
  220 CONTINUE
      NNNX1 = NNNX + 1
      DO 230 I = NNNX1, NNN
         LXN(2, I) = -LXN(2, I)
  230 CONTINUE
      RETURN

10000 FORMAT(' IN NCKLCE, THE PERIMETER IS NOT CONTINUOUS AT NODE', I5)

      END
