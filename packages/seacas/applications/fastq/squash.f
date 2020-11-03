C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SQUASH (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &   KKKOLD, NNN, NAVAIL, IAVAIL, ASMALL, DONE, NOROOM, ERR)
C***********************************************************************

C  SUBROUTINE SQUASH = TRIES TO SQUASH THE MOST SHARPELY RHOMBICAL
C                      ELEMENT IN THE MESH

C***********************************************************************

C  VARIABLES:
C     ASMALL  =  AN ANGLE CONSIDERED TO BE TOO SMALL (DEGREES)
C     DONE    =  .TRUE. IF A SQUASH WAS DONE

C***********************************************************************

      DIMENSION NODES(4), LINES(4), L1LIST(20)
      DIMENSION RHOM(6), KR(6)
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), NUID(MXND)

      LOGICAL ERR, NOROOM, DONE, CCW

C  CHECK FOR IMPENDING OVERFLOW

      NOROOM = .FALSE.
      ERR = .FALSE.
      DONE = .FALSE.
      IF (NNN .GE. MXND - 1) THEN
         CALL MESAGE ('INSUFFICIENT STORAGE IN SQUASH')
         NOROOM = .TRUE.
         ERR = .TRUE.
         RETURN
      END IF

C  INITIALIZATION

      PI = ATAN2(0.0, -1.0)
      IF (ASMALL .GT. 80.0)ASMALL = 80.0
      ALOP = ASMALL*PI/180.
      DMAX = SIN(0.5*ALOP)/COS(0.5*ALOP)
      DMAX = DMAX**2

C  CLEAR RECORDS

      DO 100 I = 1, 6
         RHOM(I) = 10.0
         KR(I) = 0
  100 CONTINUE

C  ELEMENT LOOP

      DO 130 K = KKKOLD + 1, KKK
         IF (LXK(1, K) .GT. 0) THEN
            CCW = .FALSE.
            CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
            N1 = NODES(1)
            N2 = NODES(2)
            N3 = NODES(3)
            N4 = NODES(4)
            DIAG1 = (XN(N3) - XN(N1))**2 + (YN(N3) - YN(N1))**2
            DIAG2 = (XN(N4) - XN(N2))**2 + (YN(N4) - YN(N2))**2
            IF (DIAG2 .NE. 0.0) THEN
               DRAT = DIAG1/DIAG2
               IF (DRAT .GT. 1.0) DRAT = 1.0/DRAT
            ELSE
               DRAT = 0.
            END IF

C  RECORD WORST FIVE

            IF (DRAT .LE. DMAX) THEN
               KR(6) = K
               RHOM(6) = DRAT
               DO 110 I = 1, 5
                  J = 7 - I
                  IF (RHOM(J) .GE. RHOM(J - 1)) GO TO 120
                  TEMP = RHOM(J - 1)
                  RHOM(J - 1) = RHOM(J)
                  RHOM(J) = TEMP
                  KTEMP = KR(J - 1)
                  KR(J - 1) = KR(J)
                  KR(J) = KTEMP
  110          CONTINUE
  120          CONTINUE
            END IF
         END IF
  130 CONTINUE

C  SQUASH WORST ELEMENT NOT VIOLATING THE BOUNDARY

      DO 140 IK = 1, 5
         K = KR(IK)
         IF (K .GT. 0) THEN
            CCW = .TRUE.
            CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
            N1 = NODES(1)
            N2 = NODES(2)
            N3 = NODES(3)
            N4 = NODES(4)
            DIAG1 = (XN(N3) - XN(N1))**2 + (YN(N3) - YN(N1))**2
            DIAG2 = (XN(N4) - XN(N2))**2 + (YN(N4) - YN(N2))**2

C  ARRANGE NODES SO THE SHORT DIAGONAL IS FROM 1ST TO 3RD NODES
C  AND INSURE THAT THE NODE TO BE DELETED IS NOT A BOUNDARY NODE

            IF (DIAG1 .GT. DIAG2) CALL NXKORD (NODES, N2)
            N1 = NODES(1)
            N3 = NODES(3)
            IF (LXN(2, N1) .LE. 0) CALL NXKORD (NODES, N3)
            N1 = NODES(1)
            IF (LXN(2, N1) .GT. 0) GO TO 150
         END IF
  140 CONTINUE
      RETURN

C  PREPARE FOR THE SQUASH OF ELEMENT K

  150 CONTINUE
      DONE = .TRUE.
      N1 = NODES(1)
      N2 = NODES(2)
      N3 = NODES(3)
      N4 = NODES(4)

C  FIND THE LINES ASSOCIATED WITH THE ELEMENT TO BE DELETED

      DO 160 I = 1, 4
         J = I + 1
         IF (J .GE. 5) J = 1
         CALL FNDLNK (MXND, LXK, NXL, K, NODES(I), NODES(J), LINES(I),
     &      ERR)
         IF (ERR) RETURN
         IF (LINES(I) .EQ. 0) THEN
            ERR = .TRUE.
            RETURN
         END IF
  160 CONTINUE

C  FIND ELEMENTS ON OTHER SIDES OF THE LINES
C  K2 AND K3 ARE NEVER NEEDED

      L1 = LINES(1)
      L2 = LINES(2)
      L3 = LINES(3)
      L4 = LINES(4)
      K1 = KXL(1, L1) + KXL(2, L1) - K
      K4 = KXL(1, L4) + KXL(2, L4) - K

C  FIX LXK ARRAY
C  DISCARD L1 FOR L2 IN K1

      DO 170 I = 1, 4
         IF (LXK(I, K1) .EQ. L1) THEN
            LXK(I, K1) = L2
            GO TO 180
         END IF
  170 CONTINUE
      WRITE(*, 10000)K1, L1
      ERR = .TRUE.
      RETURN
  180 CONTINUE

C  DISCARD L4 FOR L3 IN K4

      DO 190 I = 1, 4
         IF (LXK(I, K4) .EQ. L4) THEN
            LXK(I, K4) = L3
            GO TO 200
         END IF
  190 CONTINUE
      WRITE(*, 10000)K1, L1
      ERR = .TRUE.
      RETURN
  200 CONTINUE

C  DELETE ELEMENT K

      DO 210 I = 1, 4
         LXK(I, K) = 0
  210 CONTINUE

C  FIX KXL ARRAY
C  DISCARD K FOR K1 WITH L2

      IF (KXL(1, L2) .EQ. K) THEN
         KXL(1, L2) = K1
      ELSE IF (KXL(2, L2) .EQ. K) THEN
         KXL(2, L2) = K1
      END IF

C  DISCARD K FOR K4 WITH L3

      IF (KXL(1, L3) .EQ. K) THEN
         KXL(1, L3) = K4
      ELSE IF (KXL(2, L3) .EQ. K) THEN
         KXL(2, L3) = K4
      END IF

C  DELETE L1 AND L4

      KXL(1, L1) = 0
      KXL(2, L1) = 0
      KXL(1, L4) = 0
      KXL(2, L4) = 0

C  FIX NXL ARRAY
C  DELETE L1 AND L4

      NXL(1, L1) = 0
      NXL(2, L1) = 0
      NXL(1, L4) = 0
      NXL(2, L4) = 0

C  RECONNECT ALL LINES CONNECTING TO NODE 1 TO NODE 3

      CALL GETLXN (MXND, LXN, N1, L1LIST, NL, ERR)
      IF (ERR) RETURN
      DO 220 I = 1, NL
         LL = L1LIST(I)
         IF (NXL(1, LL) .EQ. N1) THEN
            NXL(1, LL) = N3
         ELSE IF (NXL(2, LL) .EQ. N1) THEN
            NXL(2, LL) = N3
         END IF
  220 CONTINUE

C  FIX LXN ARRAY
C  UNHOOK L1 FROM N2 AND L4 FROM N4

      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N2, L1, NNN, ERR,
     &   NOROOM)
      IF (ERR) RETURN
      CALL DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N4, L4, NNN, ERR,
     &   NOROOM)
      IF (ERR) RETURN

C  ADD ALL LINES HOOKED TO N3 TO THE LIST OF LINES FOR N3

      DO 230 I = 1, NL
         LL = L1LIST(I)
         IF ((LL .NE. L1) .AND. (LL .NE. L4)) THEN
            CALL ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, N3, LL, NNN,
     &         ERR, NOROOM)
            IF (ERR) RETURN
         END IF
  230 CONTINUE

C  DELETE N1 (UNHOOK EVERYTHING FROM IT)

      DO 240 I = 1, 3
         LXN(I, N1) = 0
  240 CONTINUE
      LXN(4, N1) = IAVAIL
      IAVAIL = N1
      NAVAIL = NAVAIL + 1

C  FIX XN AND YN ARRAYS
C  DEFINE POSITION OF N3

      IF (LXN(2, N3) .GT. 0) THEN
         XN(N3) = 0.5*(XN(N1) + XN(N3))
         YN(N3) = 0.5*(YN(N1) + YN(N3))
      END IF
      NUID(N1) = 0
      RETURN

10000 FORMAT(' IN SQUASH,  ELEMENT', I5, ' DOES NOT CONTAIN LINE', I5)
      END
