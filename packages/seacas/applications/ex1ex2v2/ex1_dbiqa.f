C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIQA (NDB, OPTION,
     &   MAXQA, MAXINF, NQAREC, QAREC, NINFO, INFO,
     &   EXODUS, *)
C=======================================================================

C   --*** DBIQA *** (EXOLIB) Read QA and information records
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBIQA reads the QA records and the information records.
C   --
C   --Note that the number of QA records and information records to be read
C   --are read in this routine.
C   --
C   --This routine calls DBIV0.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'Q' to store QA records
C   --      'I' to store information records
C   --   MAXQA - IN - the maximum number of QA records to store
C   --   MAXINF - IN - the number of information records to store
C   --   NQAREC - OUT - the number of QA records; <0 if end-of-file
C   --   QAREC - OUT - the QA records containing: (if OPTION)
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   NINFO - OUT - the number of information records; <0 if end-of-file
C   --   INFO - OUT - the information records (if OPTION)
C   --   EXODUS - OUT - false if GENESIS file, true if EXODUS file so far
C   --   * - return statement if error encountered, including end-of-file;
C   --      NOT used if valid GENESIS file; message is printed
C   --
C   --Database must be positioned at start of QA records upon entry;
C   --upon exit positioned at end of information records.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER MAXQA, MAXINF
      INTEGER NQAREC
      CHARACTER*8 QAREC(4,*)
      INTEGER NINFO
      CHARACTER*80 INFO(*)
      LOGICAL EXODUS

      CHARACTER*80 ERRMSG

      EXODUS = .FALSE.
      NQAREC = -999
      NINFO = -999

      READ (NDB, END=160, ERR=170, IOSTAT=IERR) NQAREC

C ... If too many records, skip early records and store last written ones.
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'Q') .GT. 0)) THEN
         iskp = max(0, nqarec - maxqa)
         DO 110 IQA = 1, iskp
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
  110    CONTINUE
         DO 100 IQA = 1, MIN (MAXQA, NQAREC)
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
     &         (QAREC(I,IQA), I=1,4)
  100    CONTINUE
      ELSE
         DO 120 IQA = 1, NQAREC
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
  120    CONTINUE
      END IF

      READ (NDB, END=160, ERR=190, IOSTAT=IERR) NINFO

C ... If too many records, skip early records and store last written ones.
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         iskp = max(0, ninfo - maxinf)
         DO 140 I = 1, iskp
            READ (NDB, END=200, ERR=200, IOSTAT=IERR)
  140    CONTINUE
         DO 130 I = 1, MIN (MAXINF, NINFO)
            READ (NDB, END=200, ERR=200, IOSTAT=IERR) INFO(I)
  130    CONTINUE
      ELSE
         DO 150 I = 1, NINFO
            READ (NDB, END=200, ERR=200, IOSTAT=IERR)
  150    CONTINUE
      END IF

      EXODUS = .TRUE.

      CALL DBIV0 (NQAREC, NINFO)

  160 CONTINUE
      RETURN

  170 CONTINUE
      ERRMSG = 'NUMBER OF QA RECORDS'
      GOTO 210
  180 CONTINUE
      ERRMSG = 'QA RECORDS'
      GOTO 210
  190 CONTINUE
      ERRMSG = 'NUMBER OF INFORMATION RECORDS'
      GOTO 210
  200 CONTINUE
      ERRMSG = 'INFORMATION RECORDS'
      GOTO 210
  210 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
