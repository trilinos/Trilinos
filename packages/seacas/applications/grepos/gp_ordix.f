C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORDIX (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
C=======================================================================

C   --*** ORDIX *** (GJOIN) Order a list according to indices
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --ORDIX orders a list according to a list of indices.
C   --
C   --Parameters:
C   --   NORD - IN - the number of indices
C   --   IXORD - IN - the indices of the ordered items
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered list

      INTEGER IXORD(*)
      INTEGER IOLD(*)
      INTEGER ISCR(*)
      INTEGER INEW(*)

      DO 100 I = 1, LOLD
         ISCR(I) = IOLD(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END
