C     Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
C     of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C     NTESS, the U.S. Government retains certain rights in this software.
C     
C     See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RMIXINT (INLINE, IFLD, INTYP, CFIELD, IFIELD,
     &     SELMSG, MAXSEL, NUMSEL, IXSEL, MAP, *)
C=======================================================================

C     --*** RIXINT *** (BLOT) Parse selection command
C     --   Written by Amy Gilkey - revised 05/20/88
C     --
C     --RMIXINT selects the items listed in the command.  If there are no
C     --fields, all the items are selected.  If the first field is ADD, the
c     --items are added to the items already selected, otherwise only the
C     --listed items are selected.
C     --
C     --Parameters:
C     --   INLINE - IN/OUT - the parsed input line for the log file
C     --   IFLD - IN/OUT - the free-field reader index
C     --   INTYP - IN - the free-field reader field types
C     --   CFIELD - IN - the character fields
C     --   IFIELD - IN - the integer fields
C     --   SELMSG - IN - the type of item for error messages
C     --   MAXSEL - IN - the number of the maximum selected items
C     --   NUMSEL - IN/OUT - the number of selected items; set to zero
C     --      upon entry unless ADD is the first field
C     --   IXSEL - IN/OUT - the selected items
C     --   MAP - IN - map local ids to global ids.  User specifies global,
C     ixsel stores local. Map size is MAXSEL
C     --   * - return statement if error before any items selected

      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER*(*) SELMSG
      INTEGER IXSEL(*)
      INTEGER MAP(*)

      LOGICAL FFEXST, FFNUMB, FFMATC
      CHARACTER*80 ERRMSG
      INTEGER IRNG(3)
      CHARACTER*32 ISTR

      IF (.NOT. (FFEXST (IFLD, INTYP))) THEN

C     --Select all items if no fields

         NUMSEL = MAXSEL
         DO I = 1, MAXSEL
            IXSEL(I) = I
         end do

      ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN

C     --Select no items if OFF

         CALL FFADDC ('OFF', INLINE)
         NUMSEL = 0

      ELSE

C     --Reset to none selected unless ADD

         IF (FFMATC (IFLD, INTYP, CFIELD, 'ADD', 3)) THEN
            CALL FFADDC ('ADD', INLINE)
         ELSE
            IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
               ERRMSG =
     &              'Expected "OFF" or "ADD" or ' // SELMSG // ' range'
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 140
            END IF
            NUMSEL = 0
         END IF

 110     CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN

C     --Scan numeric range

            CALL FFVRNG (IFLD, INTYP, CFIELD, IFIELD,
     &           SELMSG, -MAXSEL, IRNG, *130)
            CALL FFADDV (IRNG, INLINE)

C     --Store the range selected

            DO ID = IRNG(1), IRNG(2), IRNG(3)
C     ... See if a node/element with id 'I' exists in map.
C     Store the index of where this occurs.
               I = LOCINT(ID, MAXSEL, MAP)

               if (I .eq. 0) then
                  call intstr(1, -1, ID, ISTR, LSTR);
                  ERRMSG = 'No ' // SELMSG // ' with id equal '
     *                 // ISTR(:LSTR) // ' found.'
               ELSE
                  IF (LOCINT (I, NUMSEL, IXSEL) .LE. 0) THEN
                     IF (NUMSEL .GE. MAXSEL) THEN
                        ERRMSG = 'Too many ' // SELMSG // 's selected'
                        CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
                        GOTO 130
                     END IF

                     NUMSEL = NUMSEL + 1
                     IXSEL(NUMSEL) = I
                  END IF
               END IF
            end do

            GOTO 110
         END IF

 130     CONTINUE
         IF (NUMSEL .EQ. 0) THEN
            ERRMSG = 'No ' // SELMSG // 's are selected'
            CALL PRTERR ('CMDWARN', ERRMSG(:LENSTR(ERRMSG)))
         END IF
      END IF

      RETURN

 140  CONTINUE
      RETURN 1
      END
