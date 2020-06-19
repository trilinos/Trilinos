C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION INTADD (LENLST, INTLST)
C=======================================================================
C$Id: intadd.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: intadd.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:10  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:09  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:32  gdsjaar
c Initial revision
c

C   --*** INTADD *** (ETCLIB) Add all integers in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --INTADD returns the sum of all the integers in a list.
C   --
C   --Parameters:
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be added

      INTEGER LENLST
      INTEGER INTLST(*)

      INTADD = 0
      DO 100 I = 1, LENLST
         INTADD = INTADD + INTLST(I)
  100 CONTINUE

      RETURN
      END
