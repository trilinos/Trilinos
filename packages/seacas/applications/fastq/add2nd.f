C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADD2ND (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &   BNSIZE, LNODES, X1, Y1, X2, Y2, DIST1, DIST2, NNN, LLL, KKK,
     &   N1, N2, NLOOP, SIZEIT, ERR, NOROOM, XNOLD, YNOLD, NXKOLD,
     &   LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD, NNXK, REMESH,
     &   REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE ADD2ND = ADDS A NEW ELEMENT JUTTING OUT FROM AN EXISTING
C                      LINE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND), BNSIZE (2, MXND)

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL SIZEIT, ERR, NOROOM

      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X1
      YN (NNN) = Y1
      NODE1 = NNN
      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = X2
      YN (NNN) = Y2
      NODE2 = NNN

C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE

      IF (LXN (2, N1) .LT. 0) THEN
         BNSIZE (1, NODE1) = DIST1
         BNSIZE (2, NODE1) = 1.
      ELSE
         IF (SIZEIT) THEN

C**               LOCATION SIZE AND PROJECTING FROM LOCATION SIZE.

            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X1, Y1,
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
         IF ((BNSIZE (1, N1) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NODE1) = 1.
         ELSE
            BNSIZE (2, NODE1) = DIST1 / SIZNEW
         ENDIF
      ENDIF
      IF (LXN (2, N2) .LT. 0) THEN
         BNSIZE (1, NODE2) = DIST2
         BNSIZE (2, NODE2) = 1.
      ELSE
         IF (SIZEIT) THEN
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, X2, Y2,
     &         SIZNEW)
         ELSE
            SIZNEW = BNSIZE (1, N2)
         ENDIF
         BNSIZE (1, NODE2) = SIZNEW
         IF ((BNSIZE (1, N2) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NODE2) = 1.
         ELSE
            BNSIZE (2, NODE2) = DIST2 / SIZNEW
         ENDIF
      ENDIF

C  MAKE NXL ARRAY
C  ADD THE THREE NEW LINES

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
      NXL (2, L3) = N2

C  MAKE THE NEW ELEMENT

      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N1)
      LXK (2, KKK) = L3
      LXK (3, KKK) = L2
      LXK (4, KKK) = L1
      CALL ADDKXL (MXND, KXL, KKK, L1)
      CALL ADDKXL (MXND, KXL, KKK, L2)
      CALL ADDKXL (MXND, KXL, KKK, L3)
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))

C  ZERO OUT THE LXN ARRAY

      DO 100 I = 1, 4
         LXN (I, NODE1) = 0
         LXN (I, NODE2) = 0
  100 CONTINUE

C  REDO THE LNODES ARRAY

      LNODES (1, NODE1) = 0
      LNODES (1, NODE2) = 0
      LNODES (1, N1) = 0
      LNODES (1, N2) = 0

      LNODES (2, NODE1) = N1
      LNODES (2, NODE2) = NODE1
      LNODES (2, N2) = NODE2

      LNODES (3, N1) = NODE1
      LNODES (3, NODE1) = NODE2
      LNODES (3, NODE2) = N2

      LNODES (4, NODE1) = - 1
      LNODES (4, NODE2) = - 1

      LNODES (5, N1) = L1
      LNODES (5, NODE1) = L2
      LNODES (5, NODE2) = L3

      LNODES (8, NODE1) = LNODES (8, N1) + 1
      LNODES (8, NODE2) = LNODES (8, N2) + 1

      NLOOP = NLOOP + 2

      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N2, ERR)

  110 CONTINUE
      RETURN

      END
