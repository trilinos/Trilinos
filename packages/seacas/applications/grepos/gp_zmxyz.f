C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZMXYZ (NDIM, NUMNP, IXNP, XN, YN, ZN)
C=======================================================================
C   --*** ZMXYZ *** (GJOIN) Compress coordinates
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMXYZ compresses the coordinates by removing deleted nodes.
C   --
C   --Parameters:
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN/OUT - the number of nodes; returned compressed
C   --   IXNP - IN/OUT - the index of the compressed node; 0 if deleted
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
           IXNP(INP) = JNP
        ELSE IF (IXNP(INP) .LT. 0) THEN
           IXNP(INP) = -IXNP(-IXNP(INP))
        END IF
 100  CONTINUE

      NUMNP = JNP

      END
