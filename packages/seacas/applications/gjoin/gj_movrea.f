C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MOVREA (NMOV, RFROM, RTO)
C=======================================================================
C $Id: movrea.f,v 1.1 1999/01/18 19:21:23 gdsjaar Exp $
C $Log: movrea.f,v $
C Revision 1.1  1999/01/18 19:21:23  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:07  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:06  gdsjaar
c Initial revision
c

C   --*** MOVREA *** (GJOIN) Moves real data
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --MOVREA moves real data.
C   --
C   --Parameters:
C   --   NMOV - IN - the number of reals to move
C   --   RFROM - IN - the reals to move
C   --   RTO - OUT - the array to move to

      REAL RFROM(*), RTO(*)

      DO 100 I = 1, NMOV
         RTO(I) = RFROM(I)
  100 CONTINUE

      RETURN
      END
