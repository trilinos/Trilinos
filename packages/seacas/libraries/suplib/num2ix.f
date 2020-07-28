C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NUM2IX (NSETS, NINSET, IXSET)
C=======================================================================

C   --*** NUM2IX *** (ETCLIB) Change number-in-set to set indices
C   --   Written by Amy Gilkey - revised 10/23/87
C   --
C   --NUM2IX creates the set indices given the number in each set.
C   --
C   --Parameters:
C   --   NSETS - IN - the number of sets
C   --   NINSET - IN - the number in each set
C   --   IXSET - OUT - the set indices; set i = IXSET(i-1)+1 .. IXSET(i)

      INTEGER NINSET(*)
      INTEGER IXSET(0:*)

      IXSET(0) = 0
      DO 10 I = 1, NSETS
         IXSET(I) = IXSET(I-1) + NINSET(I)
   10 CONTINUE

      RETURN
      END
