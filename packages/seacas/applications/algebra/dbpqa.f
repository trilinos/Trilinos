C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
C=======================================================================
      SUBROUTINE DBPQA (OPTION, NQAREC, QAREC, NINFO, INFO)
C=======================================================================

C   --*** DBPQA *** (EXOLIB) Print QA and information records
C   --   Written by Amy Gilkey - revised 02/02/88
C   --
C   --DBPQA displays the QA records and the information records.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --                 'Q' to print QA records
C   --                 'I' to print information records
C   --   NQAREC - IN - the number of QA records (if OPTION)
C   --   QAREC  - IN - the QA records containing: (if OPTION)
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO - IN - the number of information records (if OPTION)
C   --   INFO  - IN - the information records (if OPTION)

      include 'params.blk'

      CHARACTER*(*) OPTION
      INTEGER NQAREC
      CHARACTER*(MXSTLN) QAREC(4,*)
      INTEGER NINFO
      CHARACTER*(MXLNLN) INFO(*)

      LOGICAL ALL

      ALL = OPTION .EQ. '*'

      IF (ALL .OR. (INDEX (OPTION, 'Q') .GT. 0)) THEN
         IF (NQAREC .GT. 0) THEN
            WRITE (*, *)
            DO 100 IQA = 1, NQAREC
               call sqzstr(qarec(1, iqa), icode)
               call sqzstr(qarec(2, iqa), ivers)
               call sqzstr(qarec(3, iqa), idate)
               call sqzstr(qarec(4, iqa), itime)
               WRITE (*, 10000)qarec(1,iqa)(:icode),
     &         qarec(2, iqa)(:max(ivers,11)), qarec(3, iqa)(:idate),
     &         qarec(4, iqa)(:itime)
10000           FORMAT (1X, 'Code:  ', A, '  version  ', A,
     &                  '  on  ', A, '  at  ', A)
  100       CONTINUE
         END IF
      END IF

      IF (ALL .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         IF (NINFO .GT. 0) THEN
            WRITE (*, *)
            DO 110 I = 1, NINFO
               WRITE (*, '(1x,a)') INFO(I)(:LENSTR(INFO(I)))
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
