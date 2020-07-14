C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INILOG (LEN, LFROM, LTO)
C=======================================================================
C$Id: inilog.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: inilog.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:03  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:02  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:31  gdsjaar
c Initial revision
c

C   --*** INILOG *** (ETCLIB) Initialize all logicals in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --INILOG initializes all the logicals in a list to a specified value.
C   --
C   --Parameters:
C   --   LEN - IN - the number of logicals in the list
C   --   LFROM - IN - the initial value
C   --   LTO - OUT - the initialized list

      INTEGER LEN
      LOGICAL LFROM
      LOGICAL LTO(*)

      DO 100 I = 1, LEN
         LTO(I) = LFROM
  100 CONTINUE

      RETURN
      END
