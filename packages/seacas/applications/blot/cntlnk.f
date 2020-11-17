C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNTLNK (NELBLK, LENE, NLNKE, LENLNK, NELEMS)
C=======================================================================

C   --*** CNTLNK *** (MESH) Return link array length and number elements
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --CNTLNK returns the length of the connectivity array and the number
C   --of elements defined.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LENLNK - OUT - the length of the connectivity array
C   --   NELEMS - OUT - the number of elements

      INTEGER LENE(0:*)
      INTEGER NLNKE(*)

      NELEMS = 0
      LENLNK = 0
      DO 100 IELB = 1, NELBLK
         NUME = LENE(IELB) - LENE(IELB-1)
         NELEMS = NELEMS + NUME
         LENLNK = LENLNK + NUME * NLNKE(IELB)
  100 CONTINUE

      RETURN
      END
