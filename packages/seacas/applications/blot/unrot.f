C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE UNROT (NUM, NPROT, ROTMAT, ROTCEN,
     &   HZ, VT, PD, XN, YN, ZN)
C=======================================================================

C   --*** UNROT *** (MESH) "Unrotate" 3D coordinates
C   --   Written by Amy Gilkey - revised 07/08/86
C   --
C   --UNROT returns the original coordinates from the rotated coordinates.
C   --It multiplies by the inverse of the rotation matrix and adds the
C   --rotation center.
C   --
C   --Parameters:
C   --   NUM - IN - the number of nodes to rotate
C   --   NPROT - IN - the node numbers of the nodes to rotate
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the center of the rotation
C   --   HZ, VT, PD - IN - the rotated nodal coordinates
C   --   XN, YN, ZN - OUT - the original nodal coordinates

      INTEGER NPROT(NUM)
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL HZ(*), VT(*), PD(*)
      REAL XN(*), YN(*), ZN(*)

      DO 100 IX = 1, NUM
         INP = NPROT(IX)
         X = HZ(INP)
         Y = VT(INP)
         Z = PD(INP)
         XN(INP) = X*ROTMAT(1,1) + Y*ROTMAT(1,2) + Z*ROTMAT(1,3)
         YN(INP) = X*ROTMAT(2,1) + Y*ROTMAT(2,2) + Z*ROTMAT(2,3)
         ZN(INP) = X*ROTMAT(3,1) + Y*ROTMAT(3,2) + Z*ROTMAT(3,3)
         XN(INP) = XN(INP) + ROTCEN(1)
         YN(INP) = YN(INP) + ROTCEN(2)
         ZN(INP) = ZN(INP) + ROTCEN(3)
  100 CONTINUE

      RETURN
      END
