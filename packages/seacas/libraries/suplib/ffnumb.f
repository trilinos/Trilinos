C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION FFNUMB (IFLD, INTYP)
C=======================================================================

C   --*** FFNUMB *** (FFLIB) Return number field status
C   --   Written by Amy Gilkey - revised 08/26/86
C   --
C   --FFNUMB returns true if and only if the field is a number (an
C   --integer or a real).
C   --
C   --Parameters:
C   --   IFLD - IN - the index of the current field number
C   --   INTYP - IN - the input type from the free-field reader

      INTEGER IFLD
      INTEGER INTYP(*)

      FFNUMB = (INTYP(IFLD) .EQ. 1) .OR. (INTYP(IFLD) .EQ. 2)

      RETURN
      END
