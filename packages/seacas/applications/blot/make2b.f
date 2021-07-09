C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MAKE2B (NELBLK, LENE, IE2ELB)
C=======================================================================

C   --*** MAKE2B *** (BLOT) Create element to element block index
C   --   Written by Amy Gilkey - revised 09/08/87
C   --
C   --MAKE2B creates a list of the element block for each element.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks to read
C   --   LENE - IN - the cumulative element counts by element block
C   --   IE2ELB - OUT - the element block for each element
C   --      filled with the element block number, not ID

      INTEGER NELBLK
      INTEGER LENE(0:NELBLK)
      INTEGER IE2ELB(*)

      DO 110 IELB = 1, NELBLK
         DO 100 IEL = LENE(IELB-1)+1, LENE(IELB)
            IE2ELB(IEL) = IELB
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
