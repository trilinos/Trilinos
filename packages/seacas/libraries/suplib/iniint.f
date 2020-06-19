C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INIINT (LEN, IFROM, ITO)
C=======================================================================
C$Id: iniint.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: iniint.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:00  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:14:59  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:30  gdsjaar
c Initial revision
c

C   --*** INIINT *** (ETCLIB) Initialize all integers in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --INIINT initializes all the integers in a list to a specified value.
C   --
C   --Parameters:
C   --   LEN - IN - the number of integers in the list
C   --   IFROM - IN - the initial value
C   --   ITO - OUT - the initialized list

      INTEGER LEN
      INTEGER IFROM
      INTEGER ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM
  100 CONTINUE

      RETURN
      END
