C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
C
C FOR NODE TRACKING MODE, GET THE WINDOW CENTER COORDS FOR THE NODE
C
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
