C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LOCINT (INT, LENLST, INTLST)
C=======================================================================

C   --*** LOCINT *** (ETCLIB) Find integer in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --LOCINT returns the index of the given integer in a list of integers.
C   --If the integer is not in the list, returns 0.
C   --
C   --Parameters:
C   --   INT - IN - the integer to be searched for
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be searched

      INTEGER INT
      INTEGER LENLST
      INTEGER INTLST(LENLST)

      DO I = 1, LENLST
        IF (INT .EQ. INTLST(I)) GOTO 20
      END DO
      I = 0

 20   CONTINUE
      LOCINT = I
      RETURN
      END
