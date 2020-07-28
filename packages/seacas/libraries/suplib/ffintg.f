C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFINTG (IFLD, INTYP, IFIELD, EXPECT, IDEFVL, IVAL, *)
C=======================================================================

C   --*** FFINTG *** (FFLIB) Parse free-field integer
C   --   Written by Amy Gilkey - revised 02/24/86
C   --
C   --FFINTG parses an integer field.  A default is supplied if the
C   --field is empty.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --   INTYP - IN - the input types from the free-field reader
C   --   IFIELD - IN - the integer fields
C   --   EXPECT - IN - the value to expect string, for error
C   --   IDEFVL - IN - the default value if field is empty
C   --   IVAL - OUT - the integer value, set only if no error
C   --   * - return statement if the field is invalid; message is printed

      INTEGER IFLD
      INTEGER INTYP(*)
      INTEGER IFIELD(*)
      CHARACTER*(*) EXPECT
      INTEGER IDEFVL, IVAL

      CHARACTER*80 ERRMSG

      IF (INTYP(IFLD) .GE. 2) THEN
         IVAL = IFIELD(IFLD)
      ELSE IF (INTYP(IFLD) .LE. -1) THEN
         IVAL = IDEFVL
      ELSE
         ERRMSG = 'Expected ' // EXPECT
         CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         GOTO 100
      END IF

      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN

  100 CONTINUE
      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN 1
      END
