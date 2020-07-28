C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FFCHAR (IFLD, INTYP, CFIELD, DEFVAL, CVAL)
C=======================================================================

C   --*** FFCHAR *** (FFLIB) Parse free-field character string
C   --   Written by Amy Gilkey - revised 02/24/86
C   --
C   --FFCHAR parses a character field.  A default is supplied if the
C   --field is empty.
C   --
C   --Parameters:
C   --   IFLD - IN/OUT - the index of the current field number, incremented
C   --   INTYP - IN - the input type from the free-field reader
C   --   CFIELD - IN - the character fields
C   --   DEFVAL - IN - the default value if field is empty
C   --   CVAL - OUT - the character value

      INTEGER IFLD
      INTEGER INTYP(*)
      CHARACTER*(*) CFIELD(*)
      CHARACTER*(*) DEFVAL, CVAL

      IF (INTYP(IFLD) .GE. 0) THEN
         CVAL = CFIELD(IFLD)
      ELSE IF (INTYP(IFLD) .LE. -1) THEN
         CVAL = DEFVAL
      END IF

      IF (INTYP(IFLD) .GE. -1) IFLD = IFLD + 1
      RETURN
      END
