C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C=======================================================================
      SUBROUTINE BANNER (NOUT, QAINFO, LINE1, LINE2, LINE3)
C=======================================================================
C$Id: banner.f,v 1.5 2009/03/25 12:46:01 gdsjaar Exp $

C   --*** BANNER *** (ETCLIB) Print program banner
C   --   Written by Amy Gilkey - revised 11/24/87
C   --
C   --BANNER prints the banner at the start of any program.  The banner
C   --is printed to the standard output device or an output file.
C   --
C   --Parameters:
C   --   NOUT - IN - the output file number; 0 if standard output device
C   --   QAINFO - IN - the current program QA information:
C   --      (1) = program name
C   --      (2) = revision date
C   --      (3) = version as "QA xx.xx" or "X  xx.xx" or "   xx.xx"
C   --      (4) = program name with version appended
C   --      (5) = date of current run
C   --      (6) = time of current run
C   --   LINE1, LINE2, LINE3 - IN - three-line program description;
C   --      first blank line causes rest to be skipped

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      INTEGER NOUT
      CHARACTER*(*) QAINFO(6)
      CHARACTER*(*) LINE1, LINE2, LINE3

      CHARACTER*80 BANR
      CHARACTER*40 BLANK
      SAVE BLANK

      DATA BLANK / ' ' /

      NCEN(LEN) = MAX (1, (80 - LEN + 1) / 2)

      IF (NOUT .LE. 0) THEN
         WRITE (*, *)
      ELSE
         WRITE (NOUT, *)
      END IF

      IF (QAINFO(3) .EQ. ' ') THEN
         WRITE (BANR, '(7A)')
     &      QAINFO(1), ' Test Version'
      ELSE
         WRITE (BANR, '(7A)')
     &      QAINFO(1), ' Version ', QAINFO(3)(2:8)
      END IF
      CALL SQZSTR (BANR, LBANR)
      IF (NOUT .LE. 0) THEN
         WRITE (*, 10000) BLANK(:NCEN(LBANR+8)),
     &      '*** ', BANR(:LBANR), ' ***'
      ELSE
         WRITE (NOUT, 10000) BLANK(:NCEN(LBANR+8)),
     &      '*** ', BANR(:LBANR), ' ***'
      END IF

      LREV = 8 + LENSTR(QAINFO(2))
      IF (NOUT .LE. 0) THEN
         WRITE (*, 10000) BLANK(:NCEN(LREV)), 'Revised ', QAINFO(2)
      ELSE
         WRITE (NOUT, 10000) BLANK(:NCEN(LREV)), 'Revised ', QAINFO(2)
      END IF

      IF (NOUT .LE. 0) THEN
         WRITE (*, *)
      ELSE
         WRITE (NOUT, *)
      END IF

      L = 0
      IF (LINE1 .NE. ' ') THEN
         L = LENSTR(LINE1)
         IF (NOUT .LE. 0) THEN
            WRITE (*, 10000) BLANK(:NCEN(L)), LINE1(:L)
         ELSE
            WRITE (NOUT, 10000) BLANK(:NCEN(L)), LINE1(:L)
         END IF
      END IF
      IF (LINE2 .NE. ' ') THEN
         L = LENSTR(LINE2)
         IF (NOUT .LE. 0) THEN
            WRITE (*, 10000) BLANK(:NCEN(L)), LINE2(:L)
         ELSE
            WRITE (NOUT, 10000) BLANK(:NCEN(L)), LINE2(:L)
         END IF
      END IF
      IF (LINE3 .NE. ' ') THEN
         L = LENSTR(LINE3)
         IF (NOUT .LE. 0) THEN
            WRITE (*, 10000) BLANK(:NCEN(L)), LINE3(:L)
         ELSE
            WRITE (NOUT, 10000) BLANK(:NCEN(L)), LINE3(:L)
         END IF
      END IF
      IF (L .NE. 0) THEN
         IF (NOUT .LE. 0) THEN
            WRITE (*, *)
         ELSE
            WRITE (NOUT, *)
         END IF
      END IF

      IF (QAINFO(5)(3:3) .NE. '/') THEN
         WRITE (BANR, 10010) QAINFO(5)(:4), QAINFO(5)(5:6),
     &        QAINFO(5)(7:8), QAINFO(6)
      ELSE
         WRITE (BANR, 10020) QAINFO(5), QAINFO(6)
      ENDIF         

      CALL SQZSTR (BANR, L)
      IF (NOUT .LE. 0) THEN
         WRITE (*, 10000) BLANK(:NCEN(L)), BANR(:L)
         WRITE (*, 10030)
         WRITE (*, *)
      ELSE
         WRITE (NOUT, 10000) BLANK(:NCEN(L)), BANR(:L)
         WRITE (NOUT, 10030)
         WRITE (NOUT, *)
      END IF


      RETURN
10000  FORMAT (8A)
10010  FORMAT ('Run on ', A4, '-', A2, '-', A2, ' at ', A8)
10020  FORMAT ('Run on ', A8, ' at ', A8)
10030  FORMAT (/,15x,
     *   '==== Email seacas-help@sandia.gov for support ====')
      END
