C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHGINT (IOLD, INEW, LIST, LLIST)
C=======================================================================

C   --*** CHGINT *** (GJOIN) Changes all occurrences in list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --CHGINT changes all occurrences of an integer in a list to a given
C   --integer.
C   --
C   --Parameters:
C   --   IOLD - IN - the integer to change
C   --   INEW - IN - the integer to replace IOLD
C   --   LIST - IN/OUT - the list of integers
C   --   LLIST - IN - the length of LIST

      INTEGER LIST(LLIST)

      DO 100 I = 1, LLIST
         IF (LIST(I) .EQ. IOLD) THEN
            LIST(I) = INEW
         END IF
  100 CONTINUE

      RETURN
      END
