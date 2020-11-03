C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDJUT (MXND, XN, YN, LXK, KXL, NXL, LXN,
     +   ANGLE, LNODES, XNEW, YNEW, NNN, LLL, NOLD, NLOOP, JUTTED)
C***********************************************************************

C  SUBROUTINE ADDJUT = ADDS A NEW NODE JUTTING OUT FROM AN EXISTING
C                      NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (7, MXND)

      LOGICAL JUTTED

      NNN = NNN+1
      XN (NNN) = XNEW
      YN (NNN) = YNEW

C  MAKE LXN AND NXL ARRAYS

C  ADD THE NEW NODE'S LINES

      LLL = LLL+1
      NXL (1, LLL) = NNN
      NXL (2, LLL) = NOLD

      DO 100 I = 1, 4
         LXN (I, NNN) = 0
  100 CONTINUE

      KXL (1, LLL) = 0
      KXL (2, LLL) = 0

C  REDO THE LNODES ARRAY

      LNODES (1, NNN) = 0
      LNODES (2, NNN) = NOLD
      LNODES (3, NNN) = NOLD
      LNODES (4, NNN) = - 1
      LNODES (5, NNN) = LLL

      LNODES (1, NOLD) = 0
      LNODES (3, NOLD) = NNN
      LNODES (5, NOLD) = LLL

      NLOOP = NLOOP + 2
      JUTTED = .TRUE.

      RETURN

      END
