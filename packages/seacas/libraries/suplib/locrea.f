C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCREA (VALU, NVALUS, VALUS)
C=======================================================================

C   --*** LOCREA *** (ETCLIB) Find closest real value
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --LOCREA returns the closest real value to the given value in a list of
C   --real values (which may not be ordered).
C   --
C   --Parameters:
C   --   VALU - the value to be searched for
C   --   NVALUS - the number of values in the list
C   --   VALUS - the list of values

      REAL VALUS(*)

      DIF = ABS (VALUS(1) - VALU)
      IX = 1
      DO 10 I = 2, NVALUS
         DIFI = ABS (VALUS(I) - VALU)
         IF (DIF .GT. DIFI) THEN
            DIF = DIFI
            IX = I
         END IF
   10 CONTINUE

      LOCREA = IX

      RETURN
      END
