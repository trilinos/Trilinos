C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCAL3D (MSCTYP, ROTMAT, ROTCEN, ALMESH, D2MESH)
C=======================================================================

C   --*** SCAL3D *** (MESH) Find the 2D limits of a 3D mesh
C   --   Written by Amy Gilkey - revised 06/30/86
C   --
C   --SCAL3D takes the mesh limits and scales them in one of two ways:
C   --   1) it finds the 2D limits enclosing the smallest sphere enclosing
C   --      a cube from the maximum 3D limits
C   --   2) it finds the 2D limits of the rotated mesh
C   --The mesh limits are expanded a little.
C   --
C   --Parameters:
C   --   MSCTYP - IN - mesh scaling flag (as in /MSHLIM/)
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the mesh center for the rotation
C   --   ALMESH - IN - the mesh limits
C   --   D2MESH - OUT - the 2D rotated mesh limits

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      CHARACTER*(*) MSCTYP
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL ALMESH(KFAR), D2MESH(KTOP)

      IF (MSCTYP .EQ. 'ALL') THEN
         RMAX = MAX (ALMESH(KRGT)-ROTCEN(1), ROTCEN(1)-ALMESH(KLFT)
     &      ,        ALMESH(KTOP)-ROTCEN(2), ROTCEN(2)-ALMESH(KBOT)
     &      ,        ALMESH(KFAR)-ROTCEN(3), ROTCEN(3)-ALMESH(KNEA))
         RAD = SQRT(3.0) * RMAX
         D2MESH(KLFT) = - RAD
         D2MESH(KRGT) = + RAD
         D2MESH(KBOT) = - RAD
         D2MESH(KTOP) = + RAD

      ELSE
         CALL BL_ROTATE (1, 1, ROTMAT, ROTCEN,
     &      ALMESH(KLFT), ALMESH(KBOT), ALMESH(KNEA), X1, Y1, RDUM)
         CALL BL_ROTATE (1, 1, ROTMAT, ROTCEN,
     &      ALMESH(KRGT), ALMESH(KTOP), ALMESH(KFAR), X2, Y2, RDUM)
         D2MESH(KLFT) = MIN (X1, X2)
         D2MESH(KRGT) = MAX (X1, X2)
         D2MESH(KBOT) = MIN (Y1, Y2)
         D2MESH(KTOP) = MAX (Y1, Y2)
      END IF

      CALL EXPLIM (2, D2MESH, D2MESH)

      RETURN
      END
