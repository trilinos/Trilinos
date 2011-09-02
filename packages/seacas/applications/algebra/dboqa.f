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
      SUBROUTINE DBOQA (NDB, QAINFO, NQAREC, QAREC, NINFO, INFO)
C=======================================================================

C   --*** DBOQA *** (EXOLIB) Write QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --   Modified by for EXODUSIIV2- 9/11/95
C   --
C   --Parameters:
C   --   NDB    - IN - the database number
C   --   QAINFO - IN - the QA record for the current run of algebra
C   --   NQAREC - IN - the number of QA records written only if >= 0
C   --   QAREC  - IN - the QA records containing:
C   --                 (1) - the analysis code name
C   --                 (2) - the analysis code QA descriptor
C   --                 (3) - the analysis date
C   --                 (4) - the analysis time
C   --   NINFO  - IN - the number of information records; written only if >= 0
C   --   INFO   - IN - the information records

      include 'params.blk'

      INTEGER NDB
      CHARACTER*(MXSTLN)  QAINFO(6)
      CHARACTER*(MXSTLN)  QAREC(4,*)
      CHARACTER*(MXLNLN)  INFO(*)

      nqarec = nqarec + 1
      QAREC(1,nqarec) = QAINFO(1)
      QAREC(2,nqarec) = QAINFO(3)
      QAREC(3,nqarec) = QAINFO(5)
      QAREC(4,nqarec) = QAINFO(6)

C     Write QA records to the file ndbout
C     There is always at least one QA record
      call expqa(ndb, nqarec, qarec, ierr)

C     Write information records to the file ndbout
      IF (NINFO .GT. 0) THEN
         call expinf(ndb, ninfo, info, ierr)
      END IF

      RETURN
      END
