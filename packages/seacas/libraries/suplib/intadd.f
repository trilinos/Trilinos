C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION INTADD (LENLST, INTLST)
C=======================================================================

C   --*** INTADD *** (ETCLIB) Add all integers in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --INTADD returns the sum of all the integers in a list.
C   --
C   --Parameters:
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be added

      INTEGER LENLST
      INTEGER INTLST(*)

      INTADD = 0
      DO 100 I = 1, LENLST
         INTADD = INTADD + INTLST(I)
  100 CONTINUE

      RETURN
      END
