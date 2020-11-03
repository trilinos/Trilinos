C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
      character*32 stra

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
      character*32 stra

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
