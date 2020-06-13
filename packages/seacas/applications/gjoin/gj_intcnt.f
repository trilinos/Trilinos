C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION INTCNT (INT, LIST, LLIST)
C=======================================================================
C $Id: intcnt.f,v 1.1 1999/01/18 19:21:22 gdsjaar Exp $
C $Log: intcnt.f,v $
C Revision 1.1  1999/01/18 19:21:22  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:34:54  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:34:53  gdsjaar
c Initial revision
c

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
