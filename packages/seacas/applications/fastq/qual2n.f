C    Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QUAL2N (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   BNSIZE, ANGLE, LXN, ITEST, LTEST, QUAL, POSBL2, POSBL3, ROWCHN,
     &   SIDPIN, ISTART, IEND, IPINCH, NPINCH, ERR)
C***********************************************************************

C  SUBROTINE QUAL2 = CHECKS THE QUALITY OF A SEMICIRCLE INTERPRETATION

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION BNSIZE (2, MXND), IPINCH(4)
      DIMENSION ICOMB (MXCORN), ITEST (2), LTEST (2), LXN (4, MXND)

      LOGICAL ERR, POSBL2, POSBL3, ROWCHN, SHRUNK, SIDPIN

      REAL NICKS, NICKC

      ERR = .FALSE.

C ... See note below regarding bug...
      ISTEP = 0

C  ASSUME PERFECT QUALITY

      QUAL = 0.
      POSBL2 = .FALSE.
      POSBL3 = .FALSE.
      ROWCHN = .FALSE.
      SIDPIN = .FALSE.

C  FIRST GET THE INTERVAL LENGTHS TO THE CHOSEN CORNERS

      ILEN = 2
      CALL SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES, ICOMB,
     &   ITEST, LTEST, ERR)
      IF (ERR) GOTO 130

C  SEE IF A SEMICIRCLE INTERPRETATION IS POSSIBLE WITH
C  THESE INTERVALS

      IF ( (LTEST(1) .GE. 2) .AND. (LTEST(2) .GE. 2) ) THEN
         POSBL2 = .TRUE.
      ENDIF

C  NOT ADD UP THE NICKS FOR BAD ANGLES

      DO 100 I =1, NCORN
         NODE = LCORN (I)
         IF (ICOMB (I) .EQ. 1) THEN
            QUAL = QUAL + NICKC (ANGLE (NODE), LXN (1, NODE))
         ELSE
            QUAL = QUAL + NICKS (ANGLE (NODE), LXN (1, NODE))
         ENDIF
  100 CONTINUE

C  NOW SEE IF A TRIANGLE INTERPRETATION IS WARRANTED

      IF (LTEST (1) .GT. LTEST (2)) THEN
         I1 = ITEST (1)
         L1 = LTEST (1)
         I2 = ITEST (2)
         L2 = LTEST (2)
      ELSE
         I1 = ITEST (2)
         L1 = LTEST (2)
         I2 = ITEST (1)
         L2 = LTEST (1)
      ENDIF
      LDIF = (L1 - L2) / 2
      IF (LDIF .GT. L1 / 2) LDIF = L1 - LDIF

C  THIS TESTS THE FORCED RECTANGLE - THE TWO NEW ROWS MUST BE
C  ENDED AT CURRENT SIDE NODES

      IF (L1 .EQ. L2) THEN
         NCHG1 = JUMPLP (MXND, MLN, LNODES, I1, L1 / 2)
         NCHG2 = JUMPLP (MXND, MLN, LNODES, I2, L1 / 2)
         IF ( (SHRUNK (BNSIZE (2, NCHG1), LNODES(8, NCHG1))) .OR.
     &      (SHRUNK (BNSIZE (2, NCHG2), LNODES(8, NCHG2))) ) THEN
           IF (NPINCH .EQ. 0) THEN
             SIDPIN = .TRUE.
             POSBL2 = .TRUE.
             IPINCH(1) = NCHG1
             IPINCH(2) = NCHG2
             NPINCH = 2
             ISTART = I1
           ENDIF
         ENDIF

C  SEE IF THE ROW NEEDS ADJUSTED SO THAT A RECTANGLE REMAINS POSSIBLE
C  WITH A SIGNIFICANTLY REDUCED ELEMENT SIZE ON THE LONG SIDE

      ELSE
C ... There is a bug here since ISTEP is not defined
C     Since it has been 'kindof' working for several years,
C     we assume that ISTEP=0 will work the first time through
C     and we add that initialization to the beginning of this routine
         NCHG1 = JUMPLP (MXND, MLN, LNODES, I1, ISTEP)
         NCHG2 = JUMPLP (MXND, MLN, LNODES, I1, L1 - ISTEP)
         IF ( (SHRUNK (BNSIZE (2, NCHG1), LNODES(8, NCHG1))) .AND.
     &      (SHRUNK (BNSIZE (2, NCHG2), LNODES(8, NCHG2))) ) THEN
            ROWCHN = .TRUE.
            ISTART = I2
            IEND = I1

C  CHECK THE SIZE REDUCTIONS AND TRIANGLE INTERPRETATION

         ELSE
            DO 110 ISTEP = LDIF + 1, L1 / 2 - 1
               NCHG1 = JUMPLP (MXND, MLN, LNODES, I1, ISTEP)
               NCHG2 = JUMPLP (MXND, MLN, LNODES, I1, L1 - ISTEP)
               IF ( (SHRUNK (BNSIZE (2, NCHG1), LNODES(8, NCHG1))) .AND.
     &            (SHRUNK (BNSIZE (2, NCHG2), LNODES(8, NCHG2))) ) THEN
                  ROWCHN = .TRUE.
                  ISTART = I2
                  IEND = I1
                  GOTO 120
               ENDIF

  110       CONTINUE
  120       CONTINUE
         ENDIF
      ENDIF

  130 CONTINUE

      RETURN

      END
