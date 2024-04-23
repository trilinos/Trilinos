C Copyright(C) 1999-2020, 2023 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCREA (VALUE, NVALUES, VALUES)
C=======================================================================

C   --*** LOCREA *** (ETCLIB) Find closest real value
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --LOCREA returns the closest real value to the given value in a list of
C   --real values (which may not be ordered).
C   --
C   --Parameters:
C   --   VALUE - the value to be searched for
C   --   NVALUES - the number of values in the list
C   --   VALUES - the list of values

      REAL VALUES(*)

      DIF = ABS (VALUES(1) - VALUE)
      IX = 1
      DO 10 I = 2, NVALUES
         DIFI = ABS (VALUES(I) - VALUE)
         IF (DIF .GT. DIFI) THEN
            DIF = DIFI
            IX = I
         END IF
   10 CONTINUE

      LOCREA = IX

      RETURN
      END
