C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADD3ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &   BNSIZE, LNODES, X, Y, DIST, NNN, LLL, KKK, N1, NLOOP, SIZEIT,
     &   ERR, NOROOM, XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &   MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &   REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE ADD2ND = ADDS A NEW ELEMENT JUTTING OUT FROM AN EXISTING
C                      LINE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION X(3), Y(3)

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL SIZEIT, ERR, NOROOM

      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X(1)
      YN (NNN) = Y(1)
      NODE1 = NNN
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X(2)
      YN (NNN) = Y(2)
      NODE2 = NNN
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X(3)
      YN (NNN) = Y(3)
      NODE3 = NNN

C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE

      IF (LXN (2, N1) .LT. 0) THEN
         BNSIZE (1, NODE1) = DIST
         BNSIZE (1, NODE2) = DIST
         BNSIZE (1, NODE3) = DIST
         BNSIZE (2, NODE1) = 1.
         BNSIZE (2, NODE2) = 1.
         BNSIZE (2, NODE3) = 1.
      ELSE
         IF (SIZEIT) THEN

C**               LOCATION SIZE AND PROJECTING FROM LOCATION SIZE.

            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X(1), Y(1),
     &         SIZE1)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XN(N1),
     &         YN(N1), SIZE2)
            SIZNEW = AMIN1 (SIZE1, SIZE2)
         ELSE
            SIZNEW = BNSIZE (1, N1)
         ENDIF
         BNSIZE (1, NODE1) = SIZNEW
         BNSIZE (1, NODE2) = SIZNEW
         BNSIZE (1, NODE3) = SIZNEW
         IF ((BNSIZE (1, N1) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NODE1) = 1.
            BNSIZE (2, NODE2) = 1.
            BNSIZE (2, NODE3) = 1.
         ELSE
            BNSIZE (2, NODE1) = DIST / SIZNEW
            BNSIZE (2, NODE2) = DIST / SIZNEW
            BNSIZE (2, NODE3) = DIST / SIZNEW
         ENDIF
      ENDIF

C  MAKE NXL ARRAY
C  ADD THE FOUR NEW LINES

      LLL = LLL+1
      L1 = LLL
      NXL (1, L1) = N1
      NXL (2, L1) = NODE1

      LLL = LLL+1
      L2 = LLL
      NXL (1, L2) = NODE1
      NXL (2, L2) = NODE2

      LLL = LLL+1
      L3 = LLL
      NXL (1, L3) = NODE2
      NXL (2, L3) = NODE3

      LLL = LLL+1
      L4 = LLL
      NXL (1, L4) = NODE3
      NXL (2, L4) = N1

C  MAKE THE NEW ELEMENT

      KKK = KKK+1
      LXK (1, KKK) = L4
      LXK (2, KKK) = L3
      LXK (3, KKK) = L2
      LXK (4, KKK) = L1
      CALL ADDKXL (MXND, KXL, KKK, L1)
      CALL ADDKXL (MXND, KXL, KKK, L2)
      CALL ADDKXL (MXND, KXL, KKK, L3)
      CALL ADDKXL (MXND, KXL, KKK, L4)

C  ZERO OUT THE LXN ARRAY

      DO 100 I = 1, 4
         LXN (I, NODE1) = 0
         LXN (I, NODE2) = 0
         LXN (I, NODE3) = 0
  100 CONTINUE

C  REDO THE LNODES ARRAY

      LNODES (1, NODE1) = 0
      LNODES (1, NODE2) = 0
      LNODES (1, NODE3) = 0
      LNODES (1, N1) = 0

      LNODES (2, NODE1) = N1
      LNODES (2, NODE2) = NODE1
      LNODES (2, NODE3) = NODE2

      LNODES (3, N1) = NODE1
      LNODES (3, NODE1) = NODE2
      LNODES (3, NODE2) = NODE3
      LNODES (3, NODE3) = N1

      LNODES (4, NODE1) = - 1
      LNODES (4, NODE2) = - 1
      LNODES (4, NODE3) = - 1

      LNODES (5, N1) = L1
      LNODES (5, NODE1) = L2
      LNODES (5, NODE2) = L3
      LNODES (5, NODE3) = L4

      LNODES (8, NODE1) = LNODES (8, N1) + 1
      LNODES (8, NODE2) = LNODES (8, N1) + 1
      LNODES (8, NODE3) = LNODES (8, N1) + 1

      NLOOP = NLOOP + 4

      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)

  110 CONTINUE
      RETURN

      END
