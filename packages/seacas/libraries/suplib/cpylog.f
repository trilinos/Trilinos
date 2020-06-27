C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPYLOG (LEN, LFROM, LTO)
C=======================================================================
C$Id: cpylog.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: cpylog.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:24  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:22  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:07  gdsjaar
c Initial revision
c

C   --*** CPYLOG *** (ETCLIB) Copy all logicals in list
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --CPYLOG copies all the logicals in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of logicals in the list
C   --   LFROM - IN - the input list
C   --   LTO - OUT - the copied list

      INTEGER LEN
      LOGICAL LFROM(*), LTO(*)

      DO 100 I = 1, LEN
         LTO(I) = LFROM(I)
  100 CONTINUE

      RETURN
      END
