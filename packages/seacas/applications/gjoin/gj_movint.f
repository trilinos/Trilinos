C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MOVINT (NMOV, IFROM, ITO)
C=======================================================================

C   --*** MOVINT *** (GJOIN) Moves integer data
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --MOVINT moves integer data.
C   --
C   --Parameters:
C   --   NMOV - IN - the number of integers to move
C   --   IFROM - IN - the integers to move
C   --   ITO - OUT - the array to move to

      INTEGER IFROM(*), ITO(*)

      DO 100 I = 1, NMOV
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END
