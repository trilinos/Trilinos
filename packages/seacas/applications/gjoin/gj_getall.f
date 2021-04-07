C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETALL (MATCH, LLIST, LIST, NINSET, INSET)
C=======================================================================

C   --*** GETALL *** (GJOIN) Get all items that match
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --GETALL checks a list of values and retrieves all items which
C   --match the given values.
C   --
C   --Parameters:
C   --   MATCH - IN - the value to match
C   --   LLIST - IN - the length of LIST
C   --   LIST - IN - the list of values
C   --   NINSET - OUT - the number of matching values
C   --   INSET - OUT - the indices of the matching items

      INTEGER LIST(*)
      INTEGER INSET(*)

      NINSET = 0
      DO 100 I = 1, LLIST
         IF (LIST(I) .EQ. MATCH) THEN
            NINSET = NINSET + 1
            INSET(NINSET) = I
         END IF
  100 CONTINUE

      RETURN
      END
