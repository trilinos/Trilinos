C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHGINT (IOLD, INEW, LIST, LLIST)
C=======================================================================
C $Id: chgint.f,v 1.1 1999/01/18 19:21:20 gdsjaar Exp $
C $Log: chgint.f,v $
C Revision 1.1  1999/01/18 19:21:20  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:25  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:09:27  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:09:26  gdsjaar
c Initial revision
c

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
