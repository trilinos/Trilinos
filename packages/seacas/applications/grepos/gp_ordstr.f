C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORDSTR (NORD, IXORD, LOLD, NAME, ISCR)
C=======================================================================
C   --*** ORDSTR *** (GJOIN) Order a list of strings according to indices
C   --   Written by Greg Sjaardema - revised 07/11/90
C   --   Modified from ORDIX Written by Amy Gilkey
C   --
C   --ORDSTR orders a list of strings according to a list of indices.
C   --
C   --Parameters:
C   --   NORD - IN - the number of indices
C   --   IXORD - IN - the indices of the ordered items
C   --   LOLD - IN - the length of NAME
C   --   NAME - IN - the unordered string list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered string list

      INTEGER IXORD(*)
      CHARACTER*(*) NAME(*)
      CHARACTER*(*) ISCR(*)

      DO 100 I = 1, LOLD
         ISCR(I) = NAME(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         NAME(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END
