C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPYINT (LEN, IFROM, ITO)
C=======================================================================
C   --*** CPYINT *** (ETCLIB) Copy all integers in list
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --CPYINT copies all the integers in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of integers in the list
C   --   IFROM - IN - the input list
C   --   ITO - OUT - the copied list

      INTEGER LEN
      INTEGER IFROM(*), ITO(*)

      DO 100 I = 1, LEN-7,8
         ITO(I+0) = IFROM(I+0)
         ITO(I+1) = IFROM(I+1)
         ITO(I+2) = IFROM(I+2)
         ITO(I+3) = IFROM(I+3)
         ITO(I+4) = IFROM(I+4)
         ITO(I+5) = IFROM(I+5)
         ITO(I+6) = IFROM(I+6)
         ITO(I+7) = IFROM(I+7)
  100 CONTINUE
      do 110 J = I, LEN
         ITO(J) = IFROM(J)
 110  continue

      RETURN
      END
