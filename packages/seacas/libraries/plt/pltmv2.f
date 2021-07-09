C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMV2(UMAP,N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),
     *          QQX(*),QQY(*)

      AXX = UMAP(1)
      AYY = UMAP(4)
      AXY = UMAP(3)
      AYX = UMAP(2)
      BX = UMAP(5)
      BY = UMAP(6)
      DO 2120 I = 1,N
         PPX(I) = AXX*PX(I) + AXY*PY(I) + BX
         QQX(I) = AXX*QX(I) + AXY*QY(I) + BX
         PPY(I) = AYX*PX(I) + AYY*PY(I) + BY
         QQY(I) = AYX*QX(I) + AYY*QY(I) + BY
 2120 CONTINUE
      J = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J + 1
      KM = 1 + J/32
      J = J + JN
      CALL PLTVWV(UMAP(7),UMAP(11),JN,MASK(KM),PPX(J1),PPY(J1),QQX(J1),
     *            QQY(J1))
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
