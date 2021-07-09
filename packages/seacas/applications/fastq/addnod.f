C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDNOD (MXND, MLN, XN, YN, LXK, KXL, NXL, LXN,
     &   ANGLE, BNSIZE, LNODES, XNEW, YNEW, DIST, NNN, KKK, LLL,
     &   N0, N1, N2, AMBIG, LAMBIG, SIZEIT, ERR, NOROOM, XNOLD,
     &   YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR, MLINK, NPNOLD, NPEOLD,
     &   NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN,
     &   EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE ADDNOD = ADDS A NEW ELEMENT TO A NEW NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)

      DIMENSION XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL AMBIG, SIZEIT, ERR, NOROOM

      NNN = NNN+1
      IF (NNN .GT. MXND) THEN
         NOROOM = .TRUE.
         GOTO 110
      ENDIF
      XN (NNN) = XNEW
      YN (NNN) = YNEW

C  PUT THE BEGINNING BOUNDARY DISTANCE IN PLACE

      IF (LXN (2, N2) .LT. 0) THEN
         BNSIZE (1, NNN) = DIST
         BNSIZE (2, NNN) = 1.
      ELSE
         IF (SIZEIT) THEN

C**               LOCATION SIZE AND PROJECTING FROM LOCATION SIZE.

            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XNEW, YNEW,
     &         SIZE1)
            CALL GETSIZ (XNOLD, YNOLD, NXKOLD, LINKEG, LISTEG, BMESUR,
     &         MLINK, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &         REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, XN(N2),
     &         YN(N2), SIZE2)
            SIZNEW = AMIN1 (SIZE1, SIZE2)
         ELSE
            SIZNEW = BNSIZE (1, N2)
         ENDIF
         BNSIZE (1, NNN) = SIZNEW
         IF ((BNSIZE (1, N2) .EQ. 0.) .OR. (SIZEIT)) THEN
            BNSIZE (2, NNN) = 1.
         ELSE
            BNSIZE (2, NNN) = DIST / SIZNEW
         ENDIF
      ENDIF

C  MAKE LXN, NXL, KXL, AND LXK ARRAYS
C  FIRST, ADD THE NEW NODE'S LINES

      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = N0
      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = N2
      DO 100 I = 1, 4
         LXN (I, NNN) = 0
  100 CONTINUE

C  MAKE THE NEW ELEMENT

      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N0)
      IF (AMBIG) THEN
         LXK (2, KKK) = LAMBIG
      ELSE
         LXK (2, KKK) = LNODES (5, N1)
      ENDIF
      LXK (3, KKK) = LLL
      LXK (4, KKK) = LLL-1
      CALL ADDKXL (MXND, KXL, KKK, LLL)
      CALL ADDKXL (MXND, KXL, KKK, LLL-1)
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
      IF (AMBIG) THEN
         CALL ADDKXL (MXND, KXL, KKK, LAMBIG)
      ELSE
         CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
      ENDIF

C  REDO THE LNODES ARRAY

      LNODES (2, N2) = NNN
      LNODES (3, N0) = NNN
      LNODES (1, N0) = 0
      LNODES (1, NNN) = 0
      LNODES (1, N2) = 0
      LNODES (2, NNN) = N0
      LNODES (3, NNN) = N2
      LNODES (5, NNN) = LLL
      LNODES (5, N0) = LLL-1
      LNODES (4, NNN) = - 1
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, N1, ERR)
      IF (ERR) GOTO 110
      IF (.NOT. AMBIG) LNODES (4, N1) = - 2
      LNODES (8, NNN) = LNODES (8, N2) + 1

  110 CONTINUE
      RETURN

      END
