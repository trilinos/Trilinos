C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION NUMEQL (TORF, LENLST, LOGLST)
C=======================================================================
C$Id: numeql.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: numeql.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:50  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:49  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:40  gdsjaar
c Initial revision
c

C   --*** NUMEQL *** (ETCLIB) Count number of occurrences of logical in list
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --NUMEQL returns the number of times the given logical occurs in a list
C   --of logicals.
C   --
C   --Parameters:
C   --   TORF - IN - the logical to be counted
C   --   LENLST - IN - the number of logicals in the list
C   --   LOGLST - IN - the list of logicals to be searched

      LOGICAL TORF
      INTEGER LENLST
      LOGICAL LOGLST(*)

      NUMEQL = 0
      DO 10 I = 1, LENLST
         IF (TORF .EQV. LOGLST(I)) NUMEQL = NUMEQL + 1
   10 CONTINUE

      RETURN
      END
