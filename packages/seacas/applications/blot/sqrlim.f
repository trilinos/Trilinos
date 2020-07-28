C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SQRLIM (RDMESH, SQMESH)
C=======================================================================

C   --*** SQRLIM *** (MESH) Expand 2D mesh limits to square
C   --   Written by Amy Gilkey - revised 05/23/86
C   --
C   --SQRLIM expands the mesh limits to form a square.
C   --
C   --Parameters:
C   --   RDMESH - IN - the mesh limits
C   --      (left, right, bottom, top)
C   --   SQMESH - OUT - the square mesh limits (may be RDMESH)
C   --      (left, right, bottom, top)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      REAL RDMESH(KTOP), SQMESH(KTOP)

      RNG2 = .5 * MAX (RDMESH(KRGT) - RDMESH(KLFT)
     &   ,             RDMESH(KTOP) - RDMESH(KBOT))
      XCEN = .5 * (RDMESH(KRGT) + RDMESH(KLFT))
      YCEN = .5 * (RDMESH(KTOP) + RDMESH(KBOT))
      SQMESH(KLFT) = XCEN - RNG2
      SQMESH(KRGT) = XCEN + RNG2
      SQMESH(KBOT) = YCEN - RNG2
      SQMESH(KTOP) = YCEN + RNG2

      RETURN
      END
