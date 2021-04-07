C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QUAL4 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   ANGLE, LXN, ITEST, LTEST, QUAL, POSBL4, ERR)
C***********************************************************************

C  SUBROTINE QUAL4 = CHECKS THE QUALITY OF A RECTANGLE INTERPRETATION

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION ICOMB (MXCORN), ITEST (4), LTEST (4), LXN (4, MXND)

      LOGICAL ERR, POSBL4

      REAL NICKS, NICKC

      ERR = .FALSE.

C  ASSUME PERFECT QUALITY

      QUAL = 0.
      POSBL4 = .FALSE.

C  FIRST GET THE INTERVAL LENGTHS TO THE CHOSEN CORNERS

      ILEN = 4
      CALL SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES, ICOMB,
     &   ITEST, LTEST, ERR)
      IF (ERR) GOTO 110

C  GET SEE IF A RECTANGLE INTERPRETATION IS POSSIBLE WITH
C  THESE INTERVALS

      IF ( (LTEST(1) .EQ. LTEST(3)) .AND.
     &   (LTEST(2) .EQ. LTEST(4)) ) THEN
         POSBL4 = .TRUE.
      ELSE
         RETURN
      ENDIF

C  NOT ADD UP THE NICKS FOR BAD ANGLES

      DO 100 I =1, NCORN
         NODE = LCORN (I)
         IF (ICOMB (I) .EQ. 1) THEN
            QUAL = QUAL + (.8 * NICKC (ANGLE (NODE), LXN (1, NODE)) )
         ELSE
            QUAL = QUAL + (.8 * NICKS (ANGLE (NODE), LXN (1, NODE)) )
         ENDIF
  100 CONTINUE

  110 CONTINUE

      RETURN

      END
