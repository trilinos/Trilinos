C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE SETSTA (PROMPT, ERRMSG, NITEM1, NITEM2, IDS, ISTAT, *)
C=======================================================================

C   --*** SETSTA *** (GJOIN) Set status of items
C   --   Written by Amy Gilkey - revised 02/23/88
C   --   Revised by Greg Sjaardema -
C   --      06/12/90 - Delete multiple occurrences of ID in list
C   --
C   --SETSTA sets the status of a list of items from user-supplied
C   --instructions.
C   --
C   --Parameters:
C   --   PROMPT - IN - the prompt string
C   --   ERRMSG - IN - the item type for error messages, end with 'number'
C   --   NITEM1 - IN - the number of items in the first set
C   --   NITEM2 - IN - the number of items in the second set
C   --   IDS - IN - the IDs of the items in both sets
C   --   ISTAT - IN/OUT - the status of each item (from both sets):
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   * - return statement for print items with status

      include 'gj_filnum.blk'
      PARAMETER (MAXFLD=80)

      CHARACTER*(*) PROMPT
      CHARACTER*(*) ERRMSG
      INTEGER ISTAT(*)
      INTEGER IDS(*)

      LOGICAL MATSTR, FFEXST, FOUND
      CHARACTER*8 WORD, VERB
      CHARACTER*5 STRA
      INTEGER INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD)
      INTEGER IFIELD(MAXFLD)
      REAL RFIELD(MAXFLD)

      CHARACTER*8 CMDTBL(11)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'ID      ', 'DELETE  ', 'COMBINE ', 'RESET   ', 'CHANGE  ',
     2   'LIST    ', 'HELP    ', 'UP      ', 'EXIT    ', 'INCREMEN',
     3   '        ' /

      MATCH = 0
      NITEMS = NITEM1 + NITEM2
      IF (NITEMS .LE. 0) RETURN

  100 CONTINUE

      WRITE (*, *)
      CALL FREFLD (0, 0, PROMPT, MAXFLD,
     &   IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 160
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

      IF (VERB .EQ. 'ID') THEN
        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    'block/set number', 0, ITEM, *150)
        IF ((ITEM .LE. 0) .OR. (ITEM .GT. NITEMS)) THEN
          CALL PRTERR ('CMDERR', 'Invalid block/set number')
          GOTO 150
        END IF

        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    'new ID', 0, ID, *150)
        IDS(ITEM) = ID

        RETURN 1

      ELSE IF (VERB .EQ. 'CHANGE') THEN
        MATCH = 3
 98     CONTINUE
        FOUND = .FALSE.

        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    'block/set ID', 0, ID, *150)
        IF (ID .LE. 0) RETURN 1

        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    'new ID', 0, IDNEW, *150)

        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, 'ALL', WORD)
          IF (MATSTR(WORD, 'FIRST', 1)) THEN
            MATCH = 1
          ELSE IF (MATSTR(WORD, 'SECOND', 1)) THEN
            MATCH = 2
          ELSE IF (MATSTR(WORD, 'BOTH', 1)) THEN
            MATCH = 3
          ELSE IF (MATSTR(WORD, 'ALL', 1)) THEN
            MATCH = 3
          ELSE
            CALL PRTERR ('CMDERR',
     *        'Invalid CHANGE identifier: ' // WORD)
            GOTO 98
          END IF
        END IF

        if (match .eq. 1 .or. match .eq. 3) then
          ITEM = LOCINT (ID, NITEM1, IDS)
          IF ((ITEM .GT. 0) .AND. (ITEM .LE. NITEM1)) THEN
            FOUND = .TRUE.
            IDS(ITEM) = IDNEW
          END IF
        end if

C ... Check for same ID in second set of IDs
        if (match .eq. 2 .or. match .eq. 3) then
          IF (NITEM2 .GT. 0) THEN
            ITEM = LOCINT (ID, NITEM2, IDS(NITEM1+1))
            IF ((ITEM .GT. 0) .AND. (ITEM .LE. NITEM2)) THEN
              FOUND = .TRUE.
              IDS(NITEM1 + ITEM) = IDNEW
            END IF
          END IF
        end if

        if (FOUND .EQV. .FALSE.) then
          WRITE (STRA, '(I5)') ID
          CALL SQZSTR (STRA, LSTRA)
          CALL PRTERR ('CMDERR',
     *      'Invalid block/set ID ' // STRA(:LSTRA))
        end if

        GO TO 98

      ELSE IF (VERB .EQ. 'DELETE') THEN
 99     CONTINUE
        CALL FFINTG (IFLD, INTYP, IFIELD,
     &    'block/set ID', 0, ID, *150)
         IF (ID .LE. 0) RETURN 1
         ITEM = LOCINT (ID, NITEMS, IDS)
         IF ((ITEM .LE. 0) .OR. (ITEM .GT. NITEMS)) THEN
           WRITE (STRA, '(I5)') ID
           CALL SQZSTR (STRA, LSTRA)
           CALL PRTERR ('CMDERR',
     *       'Invalid block/set ID ' // STRA(:LSTRA))
           GOTO 99
         END IF

         IOLD = ISTAT(ITEM)
         ISTAT(ITEM) = -1

C ... Check and Remove any combines with deleted ID

         IF (IOLD .GT. 0) THEN
           NCOMB = INTCNT (IOLD, ISTAT, NITEMS)
           IF (NCOMB .EQ. 1) THEN
             INEW = LOCINT (IOLD, NITEMS, ISTAT)
             ISTAT(INEW) = 0
           ELSE IF ((NCOMB .GT. 1) .AND. (IOLD .EQ. ITEM)) THEN
             INEW = LOCINT (IOLD, NITEMS, ISTAT)
             CALL CHGINT (IOLD, INEW, ISTAT, NITEMS)
           END IF
         END IF

C ... Check for same ID in second set of IDs

         IF (NITEM2 .GT. 0) THEN
            ITEM = LOCINT (ID, NITEM2, IDS(NITEM1+1))
            IF (ITEM .GT. 0 .AND. ITEM .LE. NITEM2) THEN
               ISTAT(NITEM1 + ITEM) = -1
            END IF
         END IF
         GO TO 99

      ELSE IF (VERB .EQ. 'COMBINE') THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'block/set ID', 0, ID1, *150)
         IFIRST = LOCINT (ID1, NITEMS, IDS)
         IF ((IFIRST .LE. 0) .OR. (IFIRST .GT. NITEMS)) THEN
            WRITE (STRA, '(I5)') ID
            CALL SQZSTR (STRA, LSTRA)
            CALL PRTERR ('CMDERR',
     &         'Invalid block/set ID ' // STRA(:LSTRA))
            GOTO 150
         END IF
  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'block/set ID', 0, ID2, *120)
            if (id1 .eq. id2 .and. ifirst .le. nitem1) then
C ... If the two ids match (combine 2 2) then user probably
C     has done a reset and now wants to recombine items from
C     each database.
              ITEM = LOCINT (ID2, NITEM2, IDS(NITEM1+1)) + NITEM1
            else
              ITEM = LOCINT (ID2, NITEMS, IDS)
            end if
            IF ((ITEM .LE. 0) .OR. (ITEM .GT. NITEMS)) THEN
               WRITE (STRA, '(I5)') ID
               CALL SQZSTR (STRA, LSTRA)
               CALL PRTERR ('CMDERR',
     &            'Invalid block/set ID ' // STRA(:LSTRA))
               GOTO 120
            END IF
            ISTAT(IFIRST) = IFIRST
            ISTAT(ITEM) = IFIRST
  120       CONTINUE
            GOTO 110
         END IF

         INEW = IFIRST
         DO 130 IOLD = 1, NITEMS
            IF (ISTAT(IOLD) .EQ. INEW) THEN
               CALL CHGINT (IOLD, INEW, ISTAT, NITEMS)
            END IF
  130    CONTINUE

         RETURN 1

      ELSE IF (VERB .EQ. 'RESET') THEN
         IF (.NOT. FFEXST (IFLD, INTYP)) THEN
C         --Reset all items
            DO 140 ITEM = 1, NITEMS
               ISTAT(ITEM) = 0
  140       CONTINUE

         ELSE
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'block/set ID', 0, ID, *150)
            ITEM = LOCINT (ID, NITEMS, IDS)
            IF ((ITEM .LE. 0) .OR. (ITEM .GT. NITEMS)) THEN
               CALL PRTERR ('CMDERR', 'Invalid block/set ID')
               GOTO 150
            END IF

            IOLD = ISTAT(ITEM)
            ISTAT(ITEM) = 0

            IF (IOLD .GT. 0) THEN
               NCOMB = INTCNT (IOLD, ISTAT, NITEMS)
               IF (NCOMB .EQ. 1) THEN
                  INEW = LOCINT (IOLD, NITEMS, ISTAT)
                  ISTAT(INEW) = 0
               ELSE IF ((NCOMB .GT. 1) .AND. (IOLD .EQ. ITEM)) THEN
                  INEW = LOCINT (IOLD, NITEMS, ISTAT)
                  CALL CHGINT (IOLD, INEW, ISTAT, NITEMS)
               END IF
            END IF
         END IF

         RETURN 1

C ----------------------------------------------------------------------
C SYNTAX: INCREMENT FIRST|SECOND|BOTH|ALL increment
      ELSE IF (VERB .EQ. 'INCREMEN') THEN
        DO 142 ITEM = 1, NITEMS
          ISTAT(ITEM) = 0
 142    CONTINUE

C ... See if incrementing 'first', 'second', 'both', or 'all'
        IF (FFEXST (IFLD, INTYP)) THEN
          CALL FFCHAR (IFLD, INTYP, CFIELD, 'ALL', WORD)
          IF (MATSTR(WORD, 'FIRST', 1)) THEN
            MATCH = 1
          ELSE IF (MATSTR(WORD, 'SECOND', 1)) THEN
            MATCH = 2
          ELSE IF (MATSTR(WORD, 'BOTH', 1)) THEN
            MATCH = 3
          ELSE IF (MATSTR(WORD, 'ALL', 1)) THEN
            MATCH = 3
          ELSE
            CALL PRTERR ('CMDERR',
     *        'Invalid INCREMENT identifier: ' // WORD)
            GOTO 149
          END IF
        END IF

C ... At this point, have type of match specified. Now get increment
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'block/set ID increment', 0, INC, *149)
         ELSE
           call prterr ('CMDERR', 'ID Increment not specified')
           go to 149
         END IF

C ... Do the increment
         if (match .eq. 1 .or. match .eq. 3) then
           do 143 item = 1, nitem1
             ids(item) = ids(item) + inc
 143       continue
         end if
         if (match .eq. 2 .or. match .eq. 3) then
           do 144 item = nitem1+1, nitems
             ids(item) = ids(item) + inc
 144       continue
         end if

C ... Redo the combinations if any match
      DO 145 ITEM = 1, NITEMS
         I = LOCINT (IDS(ITEM), NITEMS, IDS)
         IF (I .LT. ITEM) THEN
            ISTAT(ITEM) = I
            ISTAT(I) = I
         END IF
 145   CONTINUE

         return 1
 149    continue
      ELSE IF (VERB .EQ. 'HELP') THEN
         WRITE (*, 10000)
10000     FORMAT (
     &      /,1X,'Valid Commands:',
     &      /,4X,'ID n newid ',
     &     '-  change the block/set ID'
     &      /,4X,'      ("n" is the block/set NUMBER, not the ID)'
     &      /,4X,'CHANGE id idnew [FIRST|SECOND|BOTH]'
     &      /,4x,'      - change the ID from id to idnew'
     &      /,4X,'DELETE id  -  delete a block/set'
     &      /,4X,'COMBINE id1 id2 ...  -  combine blocks/sets'
     &      /,4X,'INCREMENT FIRST|SECOND|BOTH increment - increment ids'
     &      /,4X,'RESET id  -  reset the block/set'
     &      /,4X,'LIST  -  list information about the blocks/sets'
     &      /,4X,'UP/EXIT  -  go up a command level'
     &      )

      ELSE IF (VERB .EQ. 'LIST') THEN
         RETURN 1

      ELSE IF (VERB .EQ. 'UP' .OR. VERB .EQ. 'EXIT') THEN
         GOTO 160

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
      END IF

  150 CONTINUE
      GOTO 100

  160 CONTINUE
      RETURN
      END
