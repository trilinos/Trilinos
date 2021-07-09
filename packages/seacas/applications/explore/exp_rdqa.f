C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDQA (NDB, NQAREC, NINFO, KQAREC, KINFO, C)
C=======================================================================

C   --*** RDQA *** (EXPLORE) Read QA and information records
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

