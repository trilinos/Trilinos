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

C $Id: rwqa.f,v 1.4 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RWQA (NTXT, NDB, C, QAINFO, *)
C=======================================================================

C   --*** RDQA *** (TXTEXO) Read QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --RDQA reads the QA records and the information records.
C   --
C   --Note that the number of QA records and information records to be read
C   --are read in this routine.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   MAXQA - IN - the maximum number of QA records to store
C   --   MAXINF - IN - the number of information records to store
C   --   NQAREC - OUT - the number of QA records; <0 if not read
C   --   QAREC - OUT - the QA records containing:
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   NINFO - OUT - the number of information records; <0 if not read
C   --   INFO - OUT - the information records
C   --   EXODUS - OUT - set false if GENESIS format, true if EXODUS so far
C   --   * - return statement if error encountered, including end-of-file;
C   --      NOT used if valid GENESIS file; message is printed
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      include 'exodusII.inc'
      CHARACTER*1 C(*)
      CHARACTER*(MXSTLN) QAINFO(6)
      
      READ (NTXT, *, END=150, ERR=150)
      READ (NTXT, *, END=150, ERR=150) NQAREC

      nqarec = nqarec + 1
      call mcrsrv('QAREC', kqarec, mxstln*4*nqarec)
      call mcstat (nerr, mem)
      if (nerr .ne. 0) then
        call memerr()
        return 1
      end if
      
      call rwqa1(ntxt, ndb, nqarec, c(kqarec), qainfo, *190)
      
      call mcdel('QAREC')
      call mcstat (nerr, mem)
      if (nerr .ne. 0) then
        call memerr()
        return 1
      end if

      READ (NTXT, *, END=170, ERR=170)
      READ (NTXT, *, END=170, ERR=170) NINFO

      call mcrsrv('INFO', kinfo, mxlnln*ninfo)
      call mcstat (nerr, mem)
      if (nerr .ne. 0) then
        call memerr()
        return 1
      end if

      call rwinfo(ntxt, ndb, ninfo, c(kinfo), *190)

      call mcdel('INFO')
      call mcstat (nerr, mem)
      if (nerr .ne. 0) then
        call memerr()
        return 1
      end if

      return
  150 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NUMBER OF QA RECORDS')
      GOTO 190
  170 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NUMBER OF INFORMATION RECORDS')
      GOTO 190
 190  continue
      return 1
      end

      subroutine rwqa1(ntxt, ndb, nqarec, qarec, qainfo, *)
      include 'exodusII.inc'
      character*(mxstln) qarec(4,*)
      character*(mxstln) qainfo(6)
      character*10 stra

      DO 100 IQA = 1, NQAREC-1
         READ (NTXT, '(A)', END=160, ERR=160)
     &      (QAREC(I,IQA), I=1,4)
  100 CONTINUE

C ... Add record for this code
      qarec(1,nqarec) = qainfo(1)
      qarec(2,nqarec) = qainfo(3)
      qarec(3,nqarec) = qainfo(5)
      qarec(4,nqarec) = qainfo(6)

      if (nqarec .gt. 0) then
        call expqa(ndb, nqarec, qarec, ierr)
        if (ierr .ne. 0) return 1
      end if
      return
  160 CONTINUE
      CALL INTSTR (1, 0, IQA, STRA, LSTRA)
      CALL PRTERR ('FATAL', 'Reading QA RECORD ' // STRA(:LSTRA))
      return 1
      end

      subroutine rwinfo(ntxt, ndb, ninfo, info, *)
      include 'exodusII.inc'
      character*(mxlnln) info(*)
      character*10 stra

      if (ninfo .le. 0) return
      DO 120 I = 1, NINFO
         READ (NTXT, '(A)', END=180, ERR=180) INFO(I)
  120 CONTINUE

      call expinf(ndb, ninfo, info, ierr)
      if (ierr .ne. 0) return 1
      return
  180 CONTINUE
      CALL INTSTR (1, 0, I, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading INFORMATION RECORD ' // STRA(:LSTRA))
      return 1
      end
