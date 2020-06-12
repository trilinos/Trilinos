C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMMAP (NUMEL, MAPEL)
C=======================================================================
C $Id: zmmap.f,v 1.1 1999/01/18 19:21:27 gdsjaar Exp $
C $Log: zmmap.f,v $
C Revision 1.1  1999/01/18 19:21:27  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:28  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:26  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:36:25  gdsjaar
c Initial revision
c

C   --*** ZMMAP *** (GJOIN) Compress element order map
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMMAP compresses the element order map by removing deleted elements.
C   --
C   --Parameters:
C   --   NUMEL - IN/OUT - the number of elements
C   --   MAPEL - IN/OUT - the element order map

      INTEGER MAPEL(*)

      JEL = 0
      DO 100 IEL = 1, NUMEL
         IF (MAPEL(IEL) .GT. 0) THEN
            JEL = JEL + 1
            MAPEL(JEL) = MAPEL(IEL)
         END IF
  100 CONTINUE

      NUMEL = JEL

      RETURN
      END
