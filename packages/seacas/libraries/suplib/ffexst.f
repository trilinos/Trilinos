C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION FFEXST (IFLD, INTYP)
C=======================================================================

C   --*** FFEXST *** (FFLIB) Return end of fields status
C   --   Written by Amy Gilkey - revised 08/26/86
C   --
C   --FFEXST returns true if and only if it has not passed the end of the
C   --parsed fields (marked by a type less than -1).
C   --
C   --Parameters:
C   --   IFLD - IN - the index of the current field number
C   --   INTYP - IN - the input type from the free-field reader;
C   --      <-1 for end of fields

      INTEGER IFLD
      INTEGER INTYP(*)

      FFEXST = (INTYP(IFLD) .GE. -1)

      RETURN
      END
