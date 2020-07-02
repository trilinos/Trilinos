C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MOVINT (NMOV, IFROM, ITO)
C=======================================================================
C $Id: movint.f,v 1.1 1999/01/18 19:21:23 gdsjaar Exp $
C $Log: movint.f,v $
C Revision 1.1  1999/01/18 19:21:23  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:04  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:03  gdsjaar
c Initial revision
c

C   --*** MOVINT *** (GJOIN) Moves integer data
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --MOVINT moves integer data.
C   --
C   --Parameters:
C   --   NMOV - IN - the number of integers to move
C   --   IFROM - IN - the integers to move
C   --   ITO - OUT - the array to move to

      INTEGER IFROM(*), ITO(*)

      DO 100 I = 1, NMOV
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END
