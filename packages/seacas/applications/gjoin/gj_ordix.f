C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORDIX (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
C=======================================================================
C $Id: ordix.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: ordix.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:22  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:21  gdsjaar
c Initial revision
c

C   --*** ORDIX *** (GJOIN) Order a list according to indices
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --ORDIX orders a list according to a list of indices.
C   --
C   --Parameters:
C   --   NORD - IN - the number of indices
C   --   IXORD - IN - the indices of the ordered items
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered list

      INTEGER IXORD(*)
      INTEGER IOLD(*)
      INTEGER ISCR(*)
      INTEGER INEW(*)

      DO 100 I = 1, LOLD
         ISCR(I) = IOLD(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END
