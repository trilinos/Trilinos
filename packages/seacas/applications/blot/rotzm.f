C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ROTZM (RDMESH, NNPSUR, NPSURF, XN, YN, ZN,
     &   ROTMSH, ROTMAT, ROTCEN, ZMMESH, ZMCEN, *)
C=======================================================================

C   --*** ROTZM *** (MESH) Find the object center within zoom limits
C   --   Written by Amy Gilkey - revised 09/09/87
C   --
C   --ROTZM averages the coordinates of all surface nodes within the zoom
C   --limits to find the center of the object visible within the zoom
C   --window.
C   --
C   --Parameters:
C   --   RDMESH - IN - the zoom mesh limits (may not be square)
C   --   NNPSUR - IN - the number of surface nodes
C   --   NPSURF - IN - the node numbers of the surface nodes
C   --   XN, YN, ZN - IN - the nodal coordinates (unrotated)
C   --   ROTMSH - IN - true iff the coordinates need to be rotated
C   --   ROTMAT - IN - the rotation matrix
C   --   ROTCEN - IN - the rotation center
C   --   ZMMESH - OUT - the zoom mesh limits (may not be square)
C   --   ZMCEN - OUT - the mesh center for the rotation
C   --   * - return statement if no nodes within window; message is printed

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      REAL RDMESH(KTOP)
      INTEGER NPSURF(*)
      REAL XN(*), YN(*), ZN(*)
      LOGICAL ROTMSH
      REAL ROTMAT(3,3), ROTCEN(3)
      REAL ZMMESH(KTOP)
      REAL ZMCEN(3)

      NIN = 0
      XTOT = 0.0
      YTOT = 0.0
      ZTOT = 0.0

      DO 100 IX = 1, NNPSUR
         INP = NPSURF(IX)
         IF (ROTMSH) THEN
            CALL BL_ROTATE (1, 1, ROTMAT, ROTCEN,
     &         XN(INP), YN(INP), ZN(INP), X, Y, Z)
         ELSE
            X = XN(INP)
            Y = YN(INP)
            Z = ZN(INP)
         END IF
         IF ((X .GE. RDMESH(KLFT)) .AND. (X .LE. RDMESH(KRGT)) .AND.
     &      (Y .GE. RDMESH(KBOT)) .AND. (Y .LE. RDMESH(KTOP))) THEN
            NIN = NIN + 1
            XTOT = XTOT + X
            YTOT = YTOT + Y
            ZTOT = ZTOT + Z
         END IF
  100 CONTINUE

      IF (NIN .LE. 0) THEN
         CALL PRTERR ('CMDERR',
     &      'No nodes within specified zoom window')
         GOTO 110
      END IF

      XTOT = XTOT / NIN
      YTOT = YTOT / NIN
      ZTOT = ZTOT / NIN

      CALL UNROT (1, 1, ROTMAT, ROTCEN,
     &   XTOT, YTOT, ZTOT, XCEN, YCEN, ZCEN)

      XDIF = ROTCEN(1) - XCEN
      YDIF = ROTCEN(2) - YCEN
      ZMMESH(KLFT) = RDMESH(KLFT) + XDIF
      ZMMESH(KRGT) = RDMESH(KRGT) + XDIF
      ZMMESH(KBOT) = RDMESH(KBOT) + YDIF
      ZMMESH(KTOP) = RDMESH(KTOP) + YDIF

      ZMCEN(1) = XCEN
      ZMCEN(2) = YCEN
      ZMCEN(3) = ZCEN

      RETURN

  110 CONTINUE
      RETURN 1
      END
