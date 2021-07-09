C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTAV2(UMAP,N,X1,Y1,X2,Y2,TH,XL)
      REAL UMAP(*)
      INTEGER N
      REAL X1(*),Y1(*)
      REAL X2(*),Y2(*)
      REAL TH
      REAL XL
      REAL PX(32),PY(32),QX(32),QY(32)
      INTEGER MASK(1)
      include 'izbit.inc'

      MASK(1) = -1
      CALL PLTMV2(UMAP,N,MASK,X1,Y1,X2,Y2,PX,PY,QX,QY)
      DO 2000 I = 1,N
         IF (IAND(MASK(1),IZBIT(I)).NE.0) THEN
            CALL PLTARR(PX(I),PY(I),QX(I),QY(I),TH,XL)
         END IF

 2000 CONTINUE
      RETURN

      END
