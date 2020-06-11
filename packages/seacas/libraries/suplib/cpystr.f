C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPYSTR (LEN, IFROM, ITO)
C=======================================================================
C$Id: cpystr.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: cpystr.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:29  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:28  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:07  gdsjaar
c Initial revision
c

C   --*** CPYSTR *** (STRLIB) Copy all strings in list
C   --   Written by Amy Gilkey - revised 03/15/88
C   --
C   --CPYSTR copies all the strings in a list to another list.
C   --
C   --Parameters:
C   --   LEN - IN - the number of strings in the list
C   --   IFROM - IN - the input list
C   --   ITO - OUT - the copied list

      INTEGER LEN
      CHARACTER*(*) IFROM(*), ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END
