C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENIX (LLIST, IOFFIX, IX, LIST, ZEROIT)
C=======================================================================

C   --*** RENIX *** (GJOIN) Renumbers items according to an index list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENIX renumbers each item in the list to the index given in another
C   --list or adds an offset to the item.
C   --
C   --Parameters:
C   --   LLIST - IN - the length of LIST
C   --   IOFFIX - IN - the offset to be added; if <0, use IX
C   --   IX - IN - the new index of the specified item
C   --   LIST - IN/OUT - the list of integers; renumbered
C   --   ZEROIT - IN - logical, if true and IX(LIST(I)) <= 0,
C               set LIST(I) = 0, else set to ABS value.
      INTEGER IX(*)
      INTEGER LIST(*)
      LOGICAL ZEROIT

      IF (IOFFIX .GE. 0) THEN
         DO 100 I = 1, LLIST
            LIST(I) = LIST(I) + IOFFIX
  100    CONTINUE
      ELSE
         DO 110 I = 1, LLIST
           IF (IX(LIST(I)) .LE. 0 .AND. ZEROIT) THEN
              LIST(I) = 0
           ELSE
              LIST(I) = ABS(IX(LIST(I)))
           END IF
  110    CONTINUE
      END IF

      RETURN
      END
