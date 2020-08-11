C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MOVREA (NMOV, RFROM, RTO)
C=======================================================================

C   --*** MOVREA *** (GJOIN) Moves real data
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --MOVREA moves real data.
C   --
C   --Parameters:
C   --   NMOV - IN - the number of reals to move
C   --   RFROM - IN - the reals to move
C   --   RTO - OUT - the array to move to

      REAL RFROM(*), RTO(*)

      DO 100 I = 1, NMOV
         RTO(I) = RFROM(I)
  100 CONTINUE

      RETURN
      END
