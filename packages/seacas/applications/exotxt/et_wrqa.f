C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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

C     $Id: wrqa.f,v 1.5 2007/10/17 18:46:10 gdsjaar Exp $
C=======================================================================
      SUBROUTINE WRQA (NTXT, NQAREC, QAREC, NINFO, INFO, QAINFO)
C=======================================================================

C   --*** WRQA *** (TXTEXO) Write QA and information records
C   --   Written by Amy Gilkey - revised 09/30/87
C   --
C   --WRQA writes the QA records and the information records.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   NQAREC - IN - the number of QA records
C   --   QAREC  - IN - the QA records containing:
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO  - IN - the number of information records
C   --   INFO   - IN - the information records
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      CHARACTER*32 QAREC(4,*)
      CHARACTER*32 QAINFO(6)
      CHARACTER*80 INFO(*)

      WRITE (NTXT, '(A)') '! QA Records'
      WRITE (NTXT, 10010) NQAREC+1, '! QA records'

      DO 100 IQA = 1, NQAREC
         WRITE (NTXT, '(A)') (QAREC(I,IQA), I=1,4)
  100 CONTINUE

C ... Add record for this code
      write (ntxt, '(A)') QAINFO(1), QAINFO(3), QAINFO(5), QAINFO(6)
      WRITE (NTXT, '(A)') '! Information Records'
      WRITE (NTXT, 10010) NINFO, '! information records'

      IF (NINFO .GT. 0) THEN
         DO 110 I = 1, NINFO
C ... Some codes are embedding newlines in the info records which screws up the text output           
C     Filter them (and other strange characters) out here...
           do 105 j=1, 80
             if (ichar(info(i)(j:j)) .lt.  32 .or.
     *           ichar(info(i)(j:j)) .gt. 126) then
               info(i)(j:j) = ' '
             end if
 105         continue
            WRITE (NTXT, 10000) INFO(I)
  110    CONTINUE
      END IF

      RETURN
10000  FORMAT (99 (A, :, 1X))
10010  FORMAT (1I10, 6X, A)
      END
