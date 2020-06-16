C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION NUMEQI (INT, LENLST, INTLST)
C=======================================================================
C$Id: numeqi.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: numeqi.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:48  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:47  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:40  gdsjaar
c Initial revision
c

C   --*** NUMEQI *** (ETCLIB) Count number of occurrences of integer in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --NUMEQI returns the number of times the given integer occurs in a list
C   --of integers.
C   --
C   --Parameters:
C   --   INT - IN - the integer to be counted
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be searched

      INTEGER INT
      INTEGER LENLST
      INTEGER INTLST(*)

      NUMEQI = 0
      DO 10 I = 1, LENLST
         IF (INT .EQ. INTLST(I)) NUMEQI = NUMEQI + 1
   10 CONTINUE

      RETURN
      END
