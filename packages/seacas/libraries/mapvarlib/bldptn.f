C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE BLDPTN(XB,YB,ZB,NDLSTB,XYZPTS)

C***********************************************************************

C BLDPTN CREATES ARRAYS XYZPTS FOR NODES

C Called by MAPVAR

C***********************************************************************

      include 'ebbyeb.blk'
      include 'amesh.blk'

      DIMENSION XB(*),YB(*),ZB(*),NDLSTB(*)
      DIMENSION XYZPTS(NUMNDB,3)

      DO 20 I = 1, NUMNDB
        XYZPTS(I,1) = XB(NDLSTB(I))
        XYZPTS(I,2) = YB(NDLSTB(I))
        IF (NDIMA .EQ. 2)THEN
          XYZPTS(I,3) = 0.
        ELSE
          XYZPTS(I,3) = ZB(NDLSTB(I))
        END IF
   20 CONTINUE

      RETURN
      END
