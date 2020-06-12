C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ORDSTR (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
C=======================================================================
C $Id: ordstr.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: ordstr.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:24  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:23  gdsjaar
c Initial revision
c

C   --*** ORDSTR *** (GJOIN) Order a list of strings according to indices
C   --   Written by Greg Sjaardema - revised 07/11/90
C   --   Modified from ORDIX Written by Amy Gilkey
C   --
C   --ORDSTR orders a list of strings according to a list of indices.
C   --
C   --Parameters:
C   --   NORD - IN - the number of indices
C   --   IXORD - IN - the indices of the ordered items
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered string list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered string list

      include 'exodusII.inc'

      INTEGER IXORD(*)
      character*(MXSTLN) iold(*)
      character*(MXSTLN) iscr(*)
      character*(MXSTLN) inew(*)

      DO 100 I = 1, LOLD
         ISCR(I) = IOLD(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END
