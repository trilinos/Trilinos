C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GETDST( NODE1, NODE2, XN, YN, ZN, DIST)
C=======================================================================

C   --*** GETDST *** (MESH) GETS DISTANCE BETWEEN TWO NODES
C   --   Written by RAY J. Meyers 21 June, 1990
C   --
C   -- calculates the distance between two nodes
C   --    (written because of referencing problem using a(xn) directly )
C   --
C   --Parameters:
C   --
C   --   NODE1 - IN - THE INDEX OF FIRST NODE
C   --   NODE2 - IN - THE INDEX OF FIRST NODE
C   --   XN - IN - THE ROTATED, DEFORMED X COORDINATES
C   --   YN - IN - THE ROTATED, DEFORMED Y COORDINATES
C   --   ZN - IN - THE ROTATED, DEFORMED Z COORDINATES
C   --   DIST - OUT - THE DISTANCE BETWEEN NODE1 AND NODE2

C=======================================================================

      include 'd3nums.blk'
      REAL XN(*), YN(*), ZN(*)

      XDIST = XN(NODE2) - XN(NODE1)
      YDIST = YN(NODE2) - YN(NODE1)
      IF(IS3DIM) THEN
         ZDIST = ZN(NODE2) - ZN(NODE1)
      ELSE
         ZDIST = 0
      END IF

      DIST = SQRT( XDIST*XDIST + YDIST*YDIST + ZDIST*ZDIST)

      RETURN
      END
