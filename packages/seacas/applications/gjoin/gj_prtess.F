C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      SUBROUTINE PRTESS (ISTAT, NESS1, NESS2, IDESS, NEESS, ISCR)
C=======================================================================
C $Id: prtess.F,v 1.2 2001/06/26 17:38:54 gdsjaar Exp $

C   --*** PRTESS *** (GJOIN) Print summary for element side sets
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTESS prints a summary for the element side sets, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each element side set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NESS1 - IN - the number of element side sets in the first database
C   --   NESS2 - IN - the number of element side sets in the second database
C   --   IDESS - IN - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   ISCR - SCRATCH - size = NESS1+NESS2

      INTEGER ISTAT(*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NUMESS = NESS1 + NESS2

      IF (NUMESS .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Element side sets:'

      WRITE (STRA, 10000) NUMESS
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 IESS = 1, NUMESS
         IF (ISTAT(IESS) .EQ. 0) THEN
            WRITE (STRA, 10000) IESS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IESS .GT. NESS1) CH = '*'
            WRITE (*, 10010) CH, IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS)
         ELSE IF (ISTAT(IESS) .GT. 0) THEN
            IF (ISTAT(IESS) .EQ. IESS) THEN
               CALL GETALL (IESS, NUMESS, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I5)') IDESS(IESS)
               CALL SQZSTR (STRB, LSTRB)
               DO 100 ISET = 1, NINSET
                  IF (ISET .EQ. 1) THEN
#if NeedsDoubleEscape
                     STR2 = '\\  '
#else
                     STR2 = '\ '
#endif
                  ELSE IF (ISET .EQ. NINSET) THEN
                     STR2 = '/ '
                  ELSE
                     STR2 = ' |'
                  END IF
                  I = ISCR(ISET)
                  WRITE (STRA, 10000) I
                  CALL PCKSTR (1, STRA)
                  CH = ' '
                  IF (I .GT. NESS1) CH = '*'
                  WRITE (*, 10010) CH, IDESS(I), STRA(:LSTRA),
     &               NEESS(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(IESS) .LT. 0) THEN
            WRITE (STRA, 10000) IESS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IESS .GT. NESS1) CH = '*'
            WRITE (*, 10010) CH, IDESS(IESS), STRA(:LSTRA),
     &         NEESS(IESS), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Set', I6, 1X, A, ':',
     &   I9, ' elements and sides', :, 3X, A, :, 3X, 5A)
      END
