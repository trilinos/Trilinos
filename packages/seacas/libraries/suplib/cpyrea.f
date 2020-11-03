C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPYREA (LEN, RFROM, RTO)
C=======================================================================

C   --*** CPYREA *** (ETCLIB) Copy all real numbers in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --CPYREA copies all the real numbers in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of real numbers in the list
C   --   RFROM - IN - the input list
C   --   RTO - OUT - the copied list

      INTEGER LEN
      REAL RFROM(*), RTO(*)

      DO 100 I = 1, LEN
         RTO(I) = RFROM(I)
  100 CONTINUE

      RETURN
      END
