C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCRL (VALUE, NVALUS, VALUOK, VALUS)
C=======================================================================

C   --*** LOCRL *** (TIMSEL) Find closest real value of selected values
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --LOCRL returns the closest real value to the given value in a list of
C   --selected real values (which may not be ordered).
C   --
C   --Parameters:
C   --   VALUE - the value to be searched for
C   --   NVALUS - the number of values in the list
C   --   VALUOK - VALUS(i) is selected iff VALUOK(i)
C   --   VALUS - the list of values

      LOGICAL VALUOK(*)
      REAL VALUS(*)

      IX = 0
      DO 100 I = 1, NVALUS
         IF (VALUOK(I)) THEN
            DIFI = ABS (VALUS(I) - VALUE)
            IF ((IX .EQ. 0) .OR. (DIF .GT. DIFI)) THEN
               DIF = DIFI
               IX = I
            END IF
         END IF
  100 CONTINUE

      LOCRL = IX

      RETURN
      END
