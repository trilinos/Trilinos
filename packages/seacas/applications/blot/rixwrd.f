C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RIXWRD (INLINE, IFLD, INTYP, CFIELD,
     &   SELMSG, LLIST, LIST, NUMSEL, IXSEL, *)
C=======================================================================

C   --*** RIXWRD *** (BLOT) Parse selection command
C   --   Written by Amy Gilkey - revised 05/17/88
C   --
C   --RIXWRD selects the items listed in the command.  If there are no
C   --fields, all items are selected.  If the first field is ADD, the
C   --items are added to the items already selected, otherwise
C   --only the listed items are selected.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD - IN/OUT - the number of the next entry to scan, incremented
C   --   INTYP - IN - the free-format field types
C   --      -1 = none, 0 = name, 1 = real, 2 = integer, 3 = character
C   --   CFIELD - IN - the scanned character line entries
C   --   SELMSG - IN - the type of item for error messages
C   --   LLIST - IN - the number of words in LIST
C   --   LIST - IN - the list of words which may be in the range
C   --   NUMSEL - OUT - the number of selected words
C   --   IXSEL - OUT - the indices of the selected words
C   --   * - return statement if error before any items selected

      include 'params.blk'
      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) SELMSG
      CHARACTER*(*) LIST(*)
      INTEGER IXSEL(*)

      LOGICAL FFEXST, FFMATC
      CHARACTER*80 ERRMSG
      CHARACTER*(MXNAME) WORD

      IF (.NOT. FFEXST (IFLD, INTYP)) THEN

C      --Select all items if no fields

         NUMSEL = LLIST
         DO 100 I = 1, LLIST
            IXSEL(I) = I
  100    CONTINUE

      ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN

C      --Select no items if OFF

         CALL FFADDC ('OFF', INLINE)
         NUMSEL = 0

      ELSE

C      --Reset to none selected unless ADD

         IF (FFMATC (IFLD, INTYP, CFIELD, 'ADD', 3)) THEN
            CALL FFADDC ('ADD', INLINE)
         ELSE
            NUMSEL = 0
         END IF

  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN

C         --Read word of range

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

C         --Find word in list and set range

            IX = LOCSTR (WORD, LLIST, LIST)
            IF (IX .LE. 0) THEN
               ERRMSG = '"' // WORD(:LENSTR(WORD)) //
     &            '" is an invalid ' // SELMSG // ', ignored'
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
            ELSE IF (LOCINT (IX, NUMSEL, IXSEL) .LE. 0) THEN
               CALL FFADDC (WORD, INLINE)
               NUMSEL = NUMSEL + 1
               IXSEL(NUMSEL) = IX
            END IF

            GOTO 110
         END IF

         IF (NUMSEL .EQ. 0) THEN
            ERRMSG = 'No ' // SELMSG // 's are selected'
            CALL PRTERR ('CMDWARN', ERRMSG(:LENSTR(ERRMSG)))
         END IF
      END IF

      RETURN
      END
