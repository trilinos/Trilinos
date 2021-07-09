C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMIX(UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      UMAP(21) = -UMAP(21)
      UMAP(24) = -UMAP(24)
      UMAP(27) = -UMAP(27)
      XC = -XC
      UMAP(18) = -UMAP(18)
      RETURN

      END
