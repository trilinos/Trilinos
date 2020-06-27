C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INITIN (LIST, LLIST, IVAL)
C=======================================================================
C $Id: initin.f,v 1.1 1999/01/18 19:21:22 gdsjaar Exp $
C $Log: initin.f,v $
C Revision 1.1  1999/01/18 19:21:22  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:34:52  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:34:51  gdsjaar
c Initial revision
c

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
