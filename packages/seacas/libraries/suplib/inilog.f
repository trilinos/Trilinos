C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INILOG (LEN, LFROM, LTO)
C=======================================================================

C   --*** INILOG *** (ETCLIB) Initialize all logicals in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --INILOG initializes all the logicals in a list to a specified value.
C   --
C   --Parameters:
C   --   LEN - IN - the number of logicals in the list
C   --   LFROM - IN - the initial value
C   --   LTO - OUT - the initialized list

      INTEGER LEN
      LOGICAL LFROM
      LOGICAL LTO(*)

      DO 100 I = 1, LEN
         LTO(I) = LFROM
  100 CONTINUE

      RETURN
      END
