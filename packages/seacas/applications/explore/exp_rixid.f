C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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
      SUBROUTINE RIXID (INLINE, IFLD, INTYP, CFIELD, IFIELD,
     &   SELMSG, MAXSEL, IDSEL, NUMSEL, IXSEL, *)
C=======================================================================

C   --*** RIXID *** (BLOT) Parse selection command
C   --
C   --RIXID selects the items listed in the command.  If there are no
C   --fields, all the items are selected.  If the first field is ADD,
C   --the items are added to the items already selected, otherwise
C   --only the listed items are selected.
C   --
C   --Parameters:
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD - IN/OUT - the free-field reader index
C   --   INTYP - IN - the free-field reader field types
C   --   CFIELD - IN - the character fields
C   --   IFIELD - IN - the integer fields
C   --   SELMSG - IN - the type of item for error messages
C   --   MAXSEL - IN - the number of the maximum selected items
C   --   IDSEL - IN - the IDs of the items to be selected
C   --   NUMSEL - IN/OUT - the number of selected items; set to zero
C   --      upon entry unless ADD is the first field
C   --   IXSEL - IN/OUT - the selected items
C   --   * - return statement if error before any items selected

      include 'exodusII.inc'
      CHARACTER*(*) INLINE
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER IFIELD(*)
      CHARACTER*(*) SELMSG
      INTEGER IDSEL(*)
      INTEGER IXSEL(*)

      LOGICAL FFEXST, FFNUMB, FFMATC
      CHARACTER*80 ERRMSG
      CHARACTER*(MXNAME) WORD

      IF (.NOT. (FFEXST (IFLD, INTYP))) THEN

C      --Select all items if no fields

         NUMSEL = MAXSEL
         DO 100 I = 1, MAXSEL
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
            IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
               ERRMSG =
     &            'Expected "OFF" or "ADD" or ' // SELMSG // ' range'
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
               GOTO 130
            END IF
            NUMSEL = 0
         END IF

  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN

C         --Scan ID

            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         SELMSG, 0, ID, *120)

C         --Find and store the index of the ID

            IX = LOCINT (ID, MAXSEL, IDSEL)

            IF (IX .LE. 0) THEN
               CALL INTSTR (1, 0, ID, WORD, LSTR)
               ERRMSG = SELMSG // ' ' //
     &            WORD(:LSTR) // ' does not exist, ignored'
               CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))

            ELSE IF (LOCINT (IX, NUMSEL, IXSEL) .LE. 0) THEN
               CALL FFADDI (ID, INLINE)
               IF (NUMSEL .GE. MAXSEL) THEN
                  ERRMSG = 'Too many ' // SELMSG // 's selected'
                  CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
                  GOTO 120
               END IF

               NUMSEL = NUMSEL + 1
               IXSEL(NUMSEL) = IX
            END IF

            GOTO 110
         END IF

  120    CONTINUE
         IF (NUMSEL .EQ. 0) THEN
            ERRMSG = 'No ' // SELMSG // 's are selected'
            CALL PRTERR ('CMDWARN', ERRMSG(:LENSTR(ERRMSG)))
         END IF
      END IF

      RETURN

  130 CONTINUE
      RETURN 1
      END
