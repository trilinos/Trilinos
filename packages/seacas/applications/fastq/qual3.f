C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QUAL3 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   ANGLE, LXN, ITEST, LTEST, QUAL, POSBL3, POSBL4, ERR)
C***********************************************************************

C  SUBROTINE QUAL3 = CHECKS THE QUALITY OF A TRIANGLE INTERPRETATION

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION ICOMB (MXCORN), ITEST (3), LTEST (3), LXN (4, MXND)

      LOGICAL ERR, POSBL3, POSBL4

      REAL NICKS, NICKC

      ERR = .FALSE.

C  ASSUME PERFECT QUALITY

      QUAL = 0.
      POSBL3 = .FALSE.
      POSBL4 = .FALSE.

C  FIRST GET THE INTERVAL LENGTHS TO THE CHOSEN CORNERS

      ILEN = 3
      CALL SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES, ICOMB,
     &   ITEST, LTEST, ERR)
      IF (ERR) GOTO 110

C  SEE IF A TRIANGLE INTERPRETATION IS POSSIBLE WITH THESE INTERVALS

      MMAX = MAX0 (LTEST(1), LTEST(2), LTEST(3))
      IF (LTEST(1) .EQ. MMAX) THEN
         IBIG = LTEST(1)
         ISUM = LTEST(2) + LTEST(3)
      ELSEIF (LTEST(2) .EQ. MMAX) THEN
         IBIG = LTEST(2)
         ISUM = LTEST(3) + LTEST(1)
      ELSE
         IBIG = LTEST(3)
         ISUM = LTEST(1) + LTEST(2)
      ENDIF

      IDIFF = ISUM - IBIG

      IF (IDIFF .LT. 0) THEN
         RETURN
      ELSEIF (IDIFF .LT. 2) THEN
         POSBL4 = .TRUE.
      ELSE
         POSBL3 = .TRUE.
      ENDIF

C  NOW ADD UP THE NICKS FOR BAD ANGLES AND BAD CONNECTIVITY

      DO 100 I = 1, NCORN
         NODE = LCORN (I)
         IF (ICOMB (I) .EQ. 1) THEN
            QUAL = QUAL + NICKC (ANGLE (NODE), LXN (1, NODE))
         ELSE
            QUAL = QUAL + NICKS (ANGLE (NODE), LXN (1, NODE))
         ENDIF
  100 CONTINUE

      IF (POSBL4) QUAL = QUAL * 1.3
      MMIN = MIN0 (LTEST(1), LTEST(2), LTEST(3))
      IF (MMIN.EQ.1) QUAL = QUAL * 1.3

  110 CONTINUE

      RETURN

      END
