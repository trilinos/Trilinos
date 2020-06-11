C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RENIX (LLIST, IOFFIX, IX, LIST)
C=======================================================================
C $Id: renix.f,v 1.1 1999/01/18 19:21:26 gdsjaar Exp $
C $Log: renix.f,v $
C Revision 1.1  1999/01/18 19:21:26  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:02  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:36:01  gdsjaar
c Initial revision
c

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

      INTEGER IX(*)
      INTEGER LIST(*)

      IF (IOFFIX .GE. 0) THEN
         DO 100 I = 1, LLIST
            LIST(I) = LIST(I) + IOFFIX
  100    CONTINUE
      ELSE
         DO 110 I = 1, LLIST
            LIST(I) = IX(LIST(I))
  110    CONTINUE
      END IF

      RETURN
      END
