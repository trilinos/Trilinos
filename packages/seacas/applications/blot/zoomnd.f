C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ZOOMND( XN, YN, ZN, ZMMESH)
C=======================================================================

C   --*** NODEZM *** (MESH) ZOOM WINDOW TO TRACK A NODE
C   --   Written by RAY J. Meyers 20 June, 1990
C   --
C   --  NODEZM repositions the zoom window to be centered on the selected
C   --         node or x,y,z position for each rendering
C   --
C   --Parameters:

C   --   XN - IN - THE ROTATED, DEFORMED X COORDINATES
C   --   YN - IN - THE ROTATED, DEFORMED Y COORDINATES
C   --   ZN - IN - THE ROTATED, DEFORMED Z COORDINATES
C   --   ZMMESH - OUT - THE NEW ZOOM MESH WINDOW COORDINATES
C   --   NOTE: IF NODEZM =0, THE ROUTINE USES XZM,YZM,ZZM INSTEAD OF
C   --         XN,YN,ZN
C   --
C   --
C   --Common Variables:
C   --   Uses XZM,YZM,ZZM,RADZM,NODEZM OF /NODZOM /
C=======================================================================

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
      REAL ZMMESH(KTOP)

      REAL XN(*), YN(*), ZN(*)

      include 'nodzom.blk'
      include 'rotopt.blk'
      include 'd3nums.blk'

C FOR NODE TRACKING MODE, GET THE WINDOW CENTER COORDS FOR THE NODE

      IF( NODEZM .NE. 0) THEN
         XCEN = XN(NODEZM)
         YCEN = YN(NODEZM)
C -- FOR X,Y,Z MODE, ROTATE XZM,YZM,ZZM COORDS
      ELSE
         IF(IS3DIM) THEN
           CALL BL_ROTATE(1, 1, ROTMAT, ROTCEN, XZM, YZM, ZZM,
     &       XCEN, YCEN, ZCEN)
         ELSE
           XCEN = XZM
           YCEN = YZM
         END IF
       END IF

C -- CALCULATE WINDOW EXTENTS

      ZMMESH(KLFT) = XCEN - RADZM
      ZMMESH(KRGT) = XCEN + RADZM
      ZMMESH(KBOT) = YCEN - RADZM
      ZMMESH(KTOP) = YCEN + RADZM

      RETURN
      END
