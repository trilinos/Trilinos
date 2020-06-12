C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MUNXYZ (NDIM, NUMNP2, NUMNP1, IXNP2, XN2, YN2, ZN2)
C=======================================================================
C $Id: munxyz.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: munxyz.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:17  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:16  gdsjaar
c Initial revision
c

C   --*** MUNXYZ *** (GJOIN) Compress coordinates from the second database
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --MUNXYZ compresses the coordinates from the second database.  The
C   --coordinates that match the nodes in the first database are deleted
C   --from the list.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP2 - IN/OUT - the number of nodes in the second database;
C   --      returned compressed
C   --   NUMNP1 - IN - the number of nodes in the first database
C   --   IXNP2 - IN/OUT - the index of the compressed node; if input as <0,
C   --      there is no matching node and the new index is returned
C   --   XN2, YN2, ZN2 - IN/OUT - the coordinates, returned compressed

      INTEGER IXNP2(*)
      REAL XN2(*), YN2(*), ZN2(*)

      JNP = 0
      DO 100 INP = 1, NUMNP2
         IF (IXNP2(INP) .LT. 0) THEN
            JNP = JNP + 1
            IXNP2(INP) = JNP + NUMNP1
            XN2(JNP) = XN2(INP)
            YN2(JNP) = YN2(INP)
            IF (NDIM .GE. 3) ZN2(JNP) = ZN2(INP)
         END IF
  100 CONTINUE

      NUMNP2 = JNP

      END
