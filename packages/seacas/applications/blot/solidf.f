C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SOLIDF (NLNKF, LINKF1, XN, YN, ZN)
C=======================================================================

C   --*** SOLIDF *** (DETOUR) Paint face
C   --   Written by Amy Gilkey - revised 09/24/85
C   --
C   --SOLIDF paints the a face of the mesh.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates

      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)
      COMMON /MSHLIM/ UNMESH(KFAR), ALMESH(KFAR),
     &   ZMMESH(KTOP), RDMESH(KTOP), TICMSH, SQMESH
      LOGICAL SQMESH

      REAL XPTS(20), YPTS(20)

      XMAX = -1.0E30
      XMIN =  1.0E30
      YMAX = -1.0E30
      YMIN =  1.0E30
      DO 100 ILINK = 1, NLNKF
         XPTS(ILINK) = XN(LINKF1(ILINK))
         YPTS(ILINK) = YN(LINKF1(ILINK))
         XMIN = MIN(XPTS(ILINK), XMIN)
         XMAX = MAX(XPTS(ILINK), XMAX)
         YMIN = MIN(YPTS(ILINK), YMIN)
         YMAX = MAX(YPTS(ILINK), YMAX)
  100 CONTINUE

      IF (XMAX .LT. ZMMESH(KLFT) .OR. XMIN .GT. ZMMESH(KRGT) .OR.
     *    YMAX .LT. ZMMESH(KBOT) .OR. YMIN .GT. ZMMESH(KTOP)) RETURN

      CALL MPD2PG (NLNKF, XPTS, YPTS, 'S')

      RETURN
      END
