C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE SETITL (TWODB)
C=======================================================================

C   --*** SETITL *** (GJOIN) Select database title
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --SETITL selects the database title from the existing titles and
C   --user-supplied instructions.
C   --
C   --Parameters:
C   --   TWODB - IN - true iff two databases

      include 'exodusII.inc'
      include 'gj_params.blk'
      include 'gj_titles.blk'
      include 'gj_filnum.blk'

      PARAMETER (MAXFLD=1)
      LOGICAL TWODB

      CHARACTER*8 WORD, VERB
      INTEGER INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER IFIELD(MAXFLD)
      REAL RFIELD(MAXFLD)

      CHARACTER*8 CMDTBL(8)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   '1       ', '2       ', 'CHANGE  ',
     2   'LIST    ', 'HELP    ', 'UP      ', 'EXIT    ',
     3   '        ' /

C   --Print the input database titles and the output database title.

   50 CONTINUE

      WRITE (*, *)
      IF ((.NOT. TWODB) .OR. (TITLE1 .EQ. TITLE2)) THEN
         WRITE (*, 55) 'Database title:'
         WRITE (*, 55) TITLE1(:LENSTR(TITLE1))
      ELSE
         WRITE (*, 55) 'Database titles:'
         WRITE (*, 55) TITLE1(:LENSTR(TITLE1))
         WRITE (*, 55) TITLE2(:LENSTR(TITLE2))
      END IF
      WRITE (*, 55) 'Output database title:'
      WRITE (*, 55) TITLE(:LENSTR(TITLE))
   55 FORMAT (1X, 5A)

  100 CONTINUE

      WRITE (*, *)
      CALL FREFLD (0, 0, 'TITLE> ', MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

      IF (IOSTAT .LT. 0) GOTO 110
      IF (NUMFLD .EQ. 0) GOTO 100
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD
      CFIELD(IFLD-1) = VERB
      INTYP(IFLD-1)  = 0
      CALL OUTLOG (KLOG, MIN(MAXFLD,NUMFLD), INTYP,
     $     CFIELD, IFIELD, RFIELD)

      IF (VERB .EQ. '1') THEN
         TITLE = TITLE1

         GOTO 50

      ELSE IF (VERB .EQ. '2') THEN
         IF (TWODB) THEN
            TITLE = TITLE2
         ELSE
            TITLE = TITLE1
         END IF

         GOTO 50

      ELSE IF (VERB .EQ. 'CHANGE') THEN
         CALL GETINP (0, 0, 'New title> ', TITLE, IOSTAT)
         INTYP(1) = 0
         CALL OUTLOG (KLOG, 1, INTYP, TITLE, IFIELD, RFIELD)

         GOTO 50

      ELSE IF (VERB .EQ. 'LIST') THEN
         GOTO 50

      ELSE IF (VERB .EQ. 'HELP') THEN
         WRITE (*, 10000)
10000     FORMAT (
     &      /,1X,'Valid Commands:',
     &      /,4X,'1  -  copy title from first database',
     &      /,4X,'2  -  copy title from second database (if any)',
     &      /,4X,'CHANGE  -  change title to user-specified title'
     &      /,4X,'LIST  -  list database titles'
     &      /,4X,'UP  -  go up a command level'
     &      )

      ELSE IF (VERB .EQ. 'UP' .OR. VERB .EQ. 'EXIT') THEN
         GOTO 110

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
      END IF

      GOTO 100

  110 CONTINUE
      RETURN
      END
