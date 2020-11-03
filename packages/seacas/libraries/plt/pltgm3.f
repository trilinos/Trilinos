C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTGM3(PX,PY,PZ,S,UMAP)
      REAL UMAP(*)
      COMMON /CENBOD/XC,YC,ZC

      UMAP(17) = 1.
      UMAP(18) = PX
      UMAP(18+1) = PY
      UMAP(18+2) = PZ + S*2.5
      UMAP(21) = 1.
      UMAP(21+1) = 0.
      UMAP(21+2) = 0.
      UMAP(24) = 0.
      UMAP(24+1) = 1.
      UMAP(24+2) = 0.
      UMAP(27) = 0.
      UMAP(27+1) = 0.
      UMAP(27+2) = -1.
      UMAP(15) = S*.01
      UMAP(16) = S*100.
      UMAP(30) = 1.
      XC = PX
      YC = PY
      ZC = PZ
      RETURN

      END
