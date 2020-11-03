C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INITIN (LIST, LLIST, IVAL)
C=======================================================================

C   --*** INITIN *** (GJOIN) Initializes a list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --INITIN initializes a list to a given value.
C   --
C   --Parameters:
C   --   LIST - OUT - the list of integers
C   --   LLIST - IN - the length of LIST
C   --   IVAL - IN - the initialization value

      INTEGER LIST(*)

      DO 100 I = 1, LLIST
         LIST(I) = IVAL
  100 CONTINUE

      RETURN
      END
