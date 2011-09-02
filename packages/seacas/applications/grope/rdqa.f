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
      SUBROUTINE RDQA (NDB, NQAREC, NINFO, KQAREC, KINFO, C)
C=======================================================================

C   --*** RDQA *** (GROPE) Read QA and information records
C   --
C   --RDQA reads the QA records and the information records.
C   --An error message is displayed if the end of file is read, unless
C   --end of GENESIS file.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NQAREC - OUT - the number of QA records
C   --   NINFO - OUT - the number of information records

      include 'exodusII.inc'

      CHARACTER*1 C(*)

      CHARACTER*80 ERRMSG

      NQAREC = 0
      NINFO = 0

      call exinq(ndb, EXQA,   nqarec, rdum, cdum, ierr)
      call exinq(ndb, EXINFO, ninfo,  rdum, cdum, ierr)
      call mcrsrv('QAREC',  kqarec, nqarec*4*MXSTLN)
      call mcrsrv('INFREC', kinfo,  ninfo*MXLNLN)
      CALL MCSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      if (nqarec .gt. 0) then
C     ... Wrapper to get strings the right length
         call exgqaw(ndb, c(kqarec), ierr)
         if (ierr .ne. 0) go to 180
      end if
      if (ninfo .gt. 0) then
C     ... Wrapper to get info record the right length
         call exginw(ndb, c(kinfo), ierr)
         if (ierr .ne. 0) go to 190
      end if

      RETURN

 160  CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'Memory Allocation'
      GOTO 210
 180  CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'QA RECORDS'
      GOTO 210
 190  CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM) 'INFORMATION RECORDS'
      GOTO 210
 210  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
 220  CONTINUE
      RETURN
      
10000 FORMAT (A)
      END
      
      subroutine exgqaw(ndb, qarec, ierr)
      include 'exodusII.inc'
      character*(mxstln) qarec(4, *)
      call exgqa(ndb, qarec, ierr)
      return
      end
      
      subroutine exginw(ndb, info, ierr)
      include 'exodusII.inc'
      character*(mxlnln) info(*)
      call exginf(ndb, info, ierr)
      return
      end
      
