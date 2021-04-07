C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
C=======================================================================

C   --*** CNTELB *** (MESH) Counts selected element blocks
C   --   Written by Amy Gilkey - revised 01/23/87
C   --
C   --CNTELB counts the number of ON element blocks and the number of selected
C   --element blocks.
C   --
C   --Parameters:
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NELBLK - IN - the number of element blocks
C   --   NUMON - OUT - the number of ON element blocks
C   --   NUMSEL - OUT - the number of selected element blocks

      INTEGER IELBST(NELBLK)

C   --Count the number of ON element blocks

      NUMON = 0
      DO 100 IELB = 1, NELBLK
         IF (IELBST(IELB) .GE. 0) NUMON = NUMON + 1
  100 CONTINUE

C   --Count the number of selected element blocks

      NUMSEL = 0
      DO 110 IELB = 1, NELBLK
         IF (IELBST(IELB) .GT. 0) NUMSEL = NUMSEL + 1
  110 CONTINUE

      RETURN
      END
