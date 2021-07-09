C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE OPNFIL (IUNIT, INOUT, FFORM, IDAPAR, IERR)
C=======================================================================
C   --*** OPNFIL *** (ETCLIB) Open a file
C   --   Written by Amy Gilkey - revised 04/20/88
C   --
C   --OPNFIL opens an input or output file.  The file name is retrieved
C   --with a call to EXNAME.  The file is rewound after open.
C   --
C   --Parameters:
C   --   IUNIT - IN - the file unit number
C   --   INOUT - IN -
C   --      'I' if input file (status = "old")
C   --      'O' if output file (status = "new")
C   --      'S' if scratch file (status = "scratch")
C   --      'U' if opened with status = "unknown"
C   --   FFORM - IN - the file type:
C   --      'F' if formatted file
C   --      'U' if unformatted file
C   --      'L' if formatted, no carriage control file
C   --   IDAPAR - IN - direct access parameters:
C   --      (1) = number of words in record; 0 if not direct access
C   --      (2) = number of records
C   --   IERR - OUT - the open error status (0 = no error)

C   --Routines Called:
C   --   EXNAME - (SUPES) Get filename from unit
C   --   LENSTR - (STRLIB) Get string length

      CHARACTER INOUT
      CHARACTER FFORM
      INTEGER IDAPAR(2)

      CHARACTER*1024 FILNAM
      CHARACTER*11 FORM
      CHARACTER*7 STAT
      CHARACTER*8 CDUM

      IERR = 999
      CALL EXNAME (IUNIT, FILNAM, LNAM)

      IF (FFORM .EQ. 'U') THEN
         FORM = 'UNFORMATTED'
      ELSE
         FORM = 'FORMATTED'
      END IF
      LFORM = LENSTR (FORM)

      IF (INOUT .EQ. 'I') THEN
         STAT = 'OLD'
      ELSE IF (INOUT .EQ. 'O') THEN
         STAT = 'NEW'
      ELSE IF (INOUT .EQ. 'U') THEN
         STAT = 'UNKNOWN'
      ELSE IF (INOUT .EQ. 'S') THEN
         STAT = 'SCRATCH'
      ELSE
         IERR = -999
         GOTO 10
      END IF
      LSTAT = LENSTR (STAT)

      NINREC = IDAPAR(1)
      IF (NINREC .EQ. 0) THEN
         OPEN (UNIT=IUNIT, FILE=FILNAM(:LNAM), FORM=FORM(:LFORM),
     &      STATUS=STAT(:LSTAT), IOSTAT = IERR)
         IF (IERR .NE. 0 .AND. INOUT .EQ. 'I') THEN
            call screrr(iunit, filnam, lnam, '...', '...')
         END IF
      ELSE
         NRECS = IDAPAR(2)
         CALL EXPARM (CDUM, CDUM, IDUM, KSCU, KNSU, IDAU)
         IF (IDAU .EQ. 0) THEN
            NICE = 512
         ELSE
            NICE = (512 * KNSU - KNSU + 1) / KSCU
         END IF
         NBLKS = (NINREC + NICE-1) / NICE
         NINREC = NBLKS * NICE
         IF (IDAU .EQ. 0) THEN
            LREC = (NINREC * KSCU + KNSU-1) / KNSU
         ELSE
            LREC = NINREC
         END IF

         IF (INOUT .EQ. 'S') THEN
            STAT = 'unknown'
            LSTAT = LENSTR (STAT)
         END IF
         NSIZE = NINREC * NRECS
         OPEN (UNIT=IUNIT, FILE=FILNAM(:LNAM), FORM=FORM(:LFORM),
     &      ACCESS='DIRECT', RECL=LREC,
     &      STATUS=STAT(:LSTAT), ERR=10)

      END IF

      IERR = 0
 10   CONTINUE
      RETURN
      END
