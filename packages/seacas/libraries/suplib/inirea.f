C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INIREA (LEN, RFROM, RTO)
C=======================================================================

C   --*** INIREA *** (ETCLIB) Initialize all real numbers in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --INIREA initializes all the real numbers in a list to a specified value.
C   --
C   --Parameters:
C   --   LEN - IN - the number of real numbers in the list
C   --   RFROM - IN - the initial value
C   --   RTO - OUT - the initialized list

      INTEGER LEN
      REAL RFROM
      REAL RTO(*)

      DO 100 I = 1, LEN
         RTO(I) = RFROM
  100 CONTINUE

      RETURN
      END
