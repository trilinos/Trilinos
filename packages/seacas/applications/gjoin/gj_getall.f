C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETALL (MATCH, LLIST, LIST, NINSET, INSET)
C=======================================================================
C $Id: getall.f,v 1.1 1999/01/18 19:21:21 gdsjaar Exp $
C $Log: getall.f,v $
C Revision 1.1  1999/01/18 19:21:21  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:25  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:34:44  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:34:43  gdsjaar
c Initial revision
c

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
