C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE COPINT (LEN, IFROM, ITO)
C=======================================================================

C   --*** COPINT *** (ETCLIB) Copy all integers in list
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --COPINT copies all the integers in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of integers in the list
C   --   IFROM - IN - the input list
C   --   ITO - IN - the copied list

      INTEGER LEN
      INTEGER IFROM(*), ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END
