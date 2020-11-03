C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION INTCNT (INT, LIST, LLIST)
C=======================================================================

C   --*** INTCNT *** (GJOIN) Returns the number of occurrences in list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --INTCNT returns the number of times an integer is found in a list.
C   --
C   --Parameters:
C   --   INT - IN - the integer to count
C   --   LIST - IN - the list of integers
C   --   LLIST - IN - the length of LIST

      INTEGER LIST(*)

      INTCNT = 0
      DO 100 I = 1, LLIST
         IF (LIST(I) .EQ. INT) THEN
            INTCNT = INTCNT + 1
         END IF
  100 CONTINUE

      RETURN
      END
