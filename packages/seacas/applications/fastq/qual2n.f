C $Id: qual2n.f,v 1.3 1998/11/24 20:45:08 gdsjaar Exp $
C $Log: qual2n.f,v $
C Revision 1.3  1998/11/24 20:45:08  gdsjaar
C Added code to avoid array bound read errors and uninitialized
C variables. In some cases, the correct fix was difficult to determine,
C so added something that looked like it made sense...
C
C This fixes problems with very slow run times on g77-compiled code. It
C was taking an uninitialized variable to be INT_MAX instead of zero
C which resulted in lots of iterations through a loop. This variable was
C initialized to zero since that is what it was being set to on the sun
C and when compiled with fort77 (f2c-based).  Gives the exact same mesh
C on linux and sun for several test cases.
C
C Revision 1.2  1995/06/28 19:21:20  gdsjaar
C Applied fixes found in memo dated May 13, 1991. The bug shows itself
C for rare cases of semicircular regions being paved.
C
c Revision 1.1.1.1  1990/11/30  11:14:15  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:14:14  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]QUAL2N.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE QUAL2N (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   BNSIZE, ANGLE, LXN, ITEST, LTEST, QUAL, POSBL2, POSBL3, ROWCHN,
     &   SIDPIN, ISTART, IEND, IPINCH, NPINCH, ERR)
C***********************************************************************
C
C  SUBROTINE QUAL2 = CHECKS THE QUALITY OF A SEMICIRCLE INTERPRETATION
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION BNSIZE (2, MXND), IPINCH(4)
      DIMENSION ICOMB (MXCORN), ITEST (2), LTEST (2), LXN (4, MXND)
C
      LOGICAL ERR, POSBL2, POSBL3, ROWCHN, SHRUNK, SIDPIN
C
      REAL NICKS, NICKC
C
      ERR = .FALSE.

C ... See note below regarding bug...
      ISTEP = 0
C
C  ASSUME PERFECT QUALITY
C
      QUAL = 0.
      POSBL2 = .FALSE.
      POSBL3 = .FALSE.
      ROWCHN = .FALSE.
      SIDPIN = .FALSE.
C
C  FIRST GET THE INTERVAL LENGTHS TO THE CHOSEN CORNERS
C
      ILEN = 2
      CALL SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES, ICOMB,
     &   ITEST, LTEST, ERR)
      IF (ERR) GOTO 130
C
C  SEE IF A SEMICIRCLE INTERPRETATION IS POSSIBLE WITH
C  THESE INTERVALS
C
      IF ( (LTEST(1) .GE. 2) .AND. (LTEST(2) .GE. 2) ) THEN
         POSBL2 = .TRUE.
      ENDIF
C
C  NOT ADD UP THE NICKS FOR BAD ANGLES
C
      DO 100 I =1, NCORN
         NODE = LCORN (I)
         IF (ICOMB (I) .EQ. 1) THEN
            QUAL = QUAL + NICKC (ANGLE (NODE), LXN (1, NODE))
         ELSE
            QUAL = QUAL + NICKS (ANGLE (NODE), LXN (1, NODE))
         ENDIF
  100 CONTINUE
C
C  NOW SEE IF A TRIANGLE INTERPRETATION IS WARRANTED
C
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
C
C  THIS TESTS THE FORCED RETANGLE - THE TWO NEW ROWS MUST BE
C  ENDED AT CURRENT SIDE NODES
C
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
C
C  SEE IF THE ROW NEEDS ADJUSTED SO THAT A RECTANGLE REMAINS POSSIBLE
C  WITH A SIGNIFICANTLY REDUCED ELEMENT SIZE ON THE LONG SIDE
C
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
C
C  CHECK THE SIZE REDUCTIONS AND TRIANGLE INTERPRETATION
C
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
C
  110       CONTINUE
  120       CONTINUE
         ENDIF
      ENDIF
C
  130 CONTINUE
C
      RETURN
C
      END
