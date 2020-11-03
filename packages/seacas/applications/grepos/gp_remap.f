C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE REMAP (LLIST, IX, LIST, ISCR)
C=======================================================================

C   --*** RENIX *** (GJOIN) Renumbers items according to an index list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENIX renumbers each item in the list to the index given in another
C   --list or adds an offset to the item.
C   --
C   --Parameters:
C   --   LLIST - IN - the length of LIST
C   --   IX - IN - the new index of the specified item
C   --   LIST - IN/OUT - the list of integers; renumbered

      INTEGER IX(*)
      INTEGER LIST(*)
      INTEGER ISCR(*)

      DO 100 I = 1, LLIST
        ISCR(I) = 0
 100  CONTINUE

      DO 110 I = LLIST, 1, -1
        INDEX = IX(I)
        if (index .ne. 0) ISCR(ABS(INDEX)) = LIST(I)
 110  CONTINUE

      DO 120 I = 1, LLIST
        LIST(I) = ISCR(I)
 120  CONTINUE

      RETURN
      END
