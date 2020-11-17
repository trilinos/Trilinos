C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CALDIS (NNENUM, NENUM, XNE, YNE, ZNE, SDISTS)
C=======================================================================

C   --*** CALDIS *** (SPLOT) Calculate node/element distances
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --CALDIS calculates the distances between selected nodes or elements.
C   --
C   --Parameters:
C   --   NNENUM - IN  - the number of selected node/element numbers
C   --   NENUM  - IN  - the selected node/element numbers
C   --   XNE,
C   --   YNE,
C   --   ZNE    - IN  - the coordinates of the nodes (if nodes)
C   --                  or element centroids (if elements)
C   --   SDISTS - OUT - the calculated distances
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      INTEGER NENUM(NNENUM)
      REAL XNE(*), YNE(*), ZNE(*)
      REAL SDISTS(NNENUM)

      SDISTS(1) = 0.0
      INE = NENUM(1)
      XTHIS = XNE(INE)
      YTHIS = YNE(INE)
      IF (NDIM .EQ. 3) ZTHIS = ZNE(INE)

      DO 100 NUM = 2, NNENUM
         INE = NENUM(NUM)
         XLAST = XTHIS
         XTHIS = XNE(INE)
         YLAST = YTHIS
         YTHIS = YNE(INE)
         IF (NDIM .EQ. 2) THEN
            SDISTS(NUM) = SDISTS(NUM-1) +
     &         SQRT ((XTHIS-XLAST)**2 + (YTHIS-YLAST)**2)
         ELSE IF (NDIM .EQ. 3) THEN
            ZLAST = ZTHIS
            ZTHIS = ZNE(INE)
            SDISTS(NUM) = SDISTS(NUM-1) +
     &         SQRT ((XTHIS-XLAST)**2 + (YTHIS-YLAST)**2
     &         + (ZTHIS-ZLAST)**2)
         END IF
  100 CONTINUE

      RETURN
      END
