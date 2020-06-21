C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

      SUBROUTINE BLDSRF(XA,YA,ZA,XYZSRF)
C
C***********************************************************************
C
C BLDSRF CREATES ARRAYS XYZSRF
C
C Called by MAPVAR
C
C***********************************************************************
C
      include 'amesh.blk'
C
      DIMENSION XA(*),YA(*),ZA(*)
      DIMENSION XYZSRF(NODESA,3)
C
C
      DO 10 I = 1, NODESA
        XYZSRF(I,1) = XA(I)
        XYZSRF(I,2) = YA(I)
        IF (NDIMA .EQ. 2)THEN
          XYZSRF(I,3) = 0.
        ELSE
          XYZSRF(I,3) = ZA(I)
        END IF
   10 CONTINUE
C
      RETURN
      END
