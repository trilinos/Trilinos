C Copyright(C) 1999-2021 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORDSTR (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
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
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered string list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered string list

      include 'exodusII.inc'

      INTEGER IXORD(*)
      character*(MXSTLN) iold(*)
      character*(MXSTLN) iscr(*)
      character*(MXSTLN) inew(*)

      DO 100 I = 1, LOLD
         ISCR(I) = IOLD(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END

C=======================================================================
      SUBROUTINE ORDNAM (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
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
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered string list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered string list

      include 'gj_namlen.blk'

      INTEGER IXORD(*)
      character*(namlen) iold(*)
      character*(namlen) iscr(*)
      character*(namlen) inew(*)

      DO I = 1, LOLD
         ISCR(I) = IOLD(I)
      end do
      DO I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
      end do

      RETURN
      END
      
