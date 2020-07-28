C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE EXPLIM (NDIM, RDMESH, EXMESH)
C=======================================================================

C   --*** EXPLIM *** (BLOT) Expand 2D or 3D mesh limits by 5%
C   --   Written by Amy Gilkey - revised 06/30/86
C   --
C   --EXPLIM expands the mesh limits by 5% of the limits of the maximum
C   --dimension (2.5% on each side).
C   --
C   --Parameters:
C   --   NDIM - IN - the number of dimensions to be expanded
C   --   RDMESH - IN - the mesh limits
C   --      (left, right, bottom, top, near, far)
C   --   EXMESH - OUT - the expanded mesh limits (may be RDMESH)
C   --      (left, right, bottom, top, near, far)

      PARAMETER (PCT2 = 0.025)
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      REAL RDMESH(2*NDIM), EXMESH(2*NDIM)

      DIF = RDMESH(2) - RDMESH(1)
      DO 100 I = 2+2, 2*NDIM, 2
         DIF = MAX (DIF, (RDMESH(I) - RDMESH(I-1)))
  100 CONTINUE
      IF (DIF .EQ. 0.0) THEN
         DO 110 I = 2, 2*NDIM, 2
            DIF = MAX (DIF, ABS (RDMESH(I)))
  110    CONTINUE
         IF (DIF .EQ. 0.0) DIF = 1.0
      END IF
      DIF = PCT2 * DIF

      DO 120 I = 2, 2*NDIM, 2
         EXMESH(I-1) = RDMESH(I-1) - DIF
         EXMESH(I) = RDMESH(I) + DIF
  120 CONTINUE

      RETURN
      END
