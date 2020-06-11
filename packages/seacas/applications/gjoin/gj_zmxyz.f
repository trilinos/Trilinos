C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMXYZ (NDIM, NUMNP, IXNP, XN, YN, ZN)
C=======================================================================
C $Id: zmxyz.f,v 1.1 1999/01/18 19:21:28 gdsjaar Exp $
C $Log: zmxyz.f,v $
C Revision 1.1  1999/01/18 19:21:28  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:28  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:31  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:36:30  gdsjaar
c Initial revision
c

C   --*** ZMXYZ *** (GJOIN) Compress coordinates
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMXYZ compresses the coordinates by removing deleted nodes.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN/OUT - the number of nodes; returned compressed
C   --   IXNP - IN - the index of the compressed node; 0 if deleted
C   --   XN, YN, ZN - IN/OUT - the coordinates, returned compressed

      INTEGER IXNP(*)
      REAL XN(*), YN(*), ZN(*)

      JNP = 0
      DO 100 INP = 1, NUMNP
         IF (IXNP(INP) .GT. 0) THEN
            JNP = JNP + 1
            XN(JNP) = XN(INP)
            YN(JNP) = YN(INP)
            IF (NDIM .GE. 3) ZN(JNP) = ZN(INP)
         END IF
  100 CONTINUE

      NUMNP = JNP

      END
