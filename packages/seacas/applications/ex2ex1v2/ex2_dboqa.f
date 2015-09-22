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
      SUBROUTINE DBOQA (NDB, NQAREC, QAREC, NINFO, INFO)
C=======================================================================
C$Id: dboqa.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dboqa.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:37  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:36  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:17  gdsjaar
c Initial revision
c 

C   --*** DBOQA *** (EXOLIB) Write QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBOQA writes the QA records and the information records.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NQAREC - IN - the number of QA records (must be at least one);
C   --      written only if >= 0
C   --   QAREC - IN - the QA records containing:
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   NINFO - IN - the number of information records; written only if >= 0
C   --   INFO - IN - the information records
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      INTEGER NDB
      INTEGER NQAREC
      CHARACTER*8 QAREC(4,*)
      INTEGER NINFO
      CHARACTER*80 INFO(*)

      IF (NQAREC .LT. 0) GOTO 120

      WRITE (NDB) NQAREC

      DO 100 IQA = 1, NQAREC
         WRITE (NDB) (QAREC(I,IQA), I=1,4)
  100 CONTINUE

      IF (NINFO .LT. 0) GOTO 120

      WRITE (NDB) NINFO

      DO 110 I = 1, NINFO
         WRITE (NDB) INFO(I)
  110 CONTINUE

  120 CONTINUE
      RETURN
      END
