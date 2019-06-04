C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C

C $Id: qual2.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: qual2.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:14:12  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:14:11  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]QUAL2.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE QUAL2 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   BNSIZE, ANGLE, LXN, ITEST, LTEST, QUAL, POSBL2, POSBL3, ROWCHN,
     &   ISTART, IEND)
C***********************************************************************
C
C  SUBROTINE QUAL2 = CHECKS THE QUALITY OF A SEMICIRCLE INTERPRETATION
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION BNSIZE (2, MXND)
      DIMENSION ICOMB (MXCORN), ITEST (2), LTEST (2), LXN (4, MXND)
C
      LOGICAL ERR, POSBL2, POSBL3, ROWCHN, SHRUNK
C
      REAL NICKS, NICKC

C ... See note below regarding bug...
      ISTEP = 0
C
C  ASSUME PERFECT QUALITY
C
      QUAL = 0.
      POSBL2 = .FALSE.
      POSBL3 = .FALSE.
      ROWCHN = .FALSE.
C
C  FIRST GET THE INTERVAL LENGTHS TO THE CHOSEN CORNERS
C
      ILEN = 2
      CALL SPACED (MXND, MXCORN, MLN, ILEN, NCORN, LCORN, LNODES, ICOMB,
     &   ITEST, LTEST, ERR)
      IF (ERR) RETURN
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
C  THIS TESTS THE FORCED TRIANGLE - THE NEW ROW MUST BE
C  ENDED AT A CURRENT SIDE NODE
C
      IF (L1 .EQ. L2) THEN
         NCHG1 = JUMPLP (MXND, MLN, LNODES, I1, LDIF)
         NCHG2 = JUMPLP (MXND, MLN, LNODES, I2, LDIF)
         QUAL13 = QUAL + NICKC (ANGLE (NCHG1), LXN (1, NCHG1))
         QUAL23 = QUAL + NICKC (ANGLE (NCHG2), LXN (1, NCHG2))
         IF ( (SHRUNK (BNSIZE (2, NCHG1), LNODES(8, NCHG1))) .AND.
     &      (SHRUNK (BNSIZE (2, NCHG2), LNODES(8, NCHG2))) ) THEN
            POSBL3 = .TRUE.
            POSBL2 = .FALSE.
            IF (QUAL13 .LT. QUAL23) THEN
               IEND = NCHG1
               ISTART = I1
            ELSE
               IEND = NCHG2
               ISTART = I2
            ENDIF
         ELSEIF (SHRUNK (BNSIZE (2, NCHG1), LNODES(8, NCHG1))) THEN
            POSBL3 = .TRUE.
            POSBL2 = .FALSE.
            IEND = NCHG1
            ISTART = I1
         ELSEIF (SHRUNK (BNSIZE (2, NCHG2), LNODES(8, NCHG2))) THEN
            POSBL3 = .TRUE.
            POSBL2 = .FALSE.
            IEND = NCHG2
            ISTART = I2
         ELSE
            POSBL3 = .FALSE.
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
      RETURN
C
      END
