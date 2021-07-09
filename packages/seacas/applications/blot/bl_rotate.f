C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE BL_ROTATE (NUM, NPROT, ROTMAT, ROTCEN,
     &   XN, YN, ZN, HZ, VT, PD)
C=======================================================================

C   --*** ROTATE *** (MESH) Rotate 3D coordinates
C   --   Written by Amy Gilkey - revised 09/09/87
C   --
C   --ROTATE rotates the 3D coordinates by subtracting the rotation center
C   --and multipling by the rotation matrix.
C   --
C   --Parameters:
C   --   NUM - IN - the number of nodes to rotate
C   --   NPROT - IN - the node numbers of the nodes to rotate
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the center of the rotation
C   --   XN, YN, ZN - IN - the original nodal coordinates
C   --   HZ, VT, PD - OUT - the rotated nodal coordinates

      INTEGER NPROT(NUM)
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL XN(*), YN(*), ZN(*)
      REAL HZ(*), VT(*), PD(*)

      DO 100 IX = 1, NUM
         INP = NPROT(IX)
         X = XN(INP) - ROTCEN(1)
         Y = YN(INP) - ROTCEN(2)
         Z = ZN(INP) - ROTCEN(3)
         HZ(INP) = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
         VT(INP) = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
         PD(INP) = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
  100 CONTINUE

      RETURN
      END
