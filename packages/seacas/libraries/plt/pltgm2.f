C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTGM2(XL,XU,YL,YU,PXL,PXU,PYL,PYU,UMAP)
      REAL UMAP(*)

      DU = XU - XL
      DP = PXU - PXL
      IF (DU.EQ.0.) THEN
         DU = 1.
      END IF

      UMAP(1) = DP/DU
      UMAP(1+1) = 0.
      UMAP(1+2) = 0.
      DV = YU - YL
      DQ = PYU - PYL
      IF (DV.EQ.0.) THEN
         DV = 1.
      END IF

      UMAP(1+3) = DQ/DV
      UMAP(5) = PXL - XL*UMAP(1)
      UMAP(5+1) = PYL - YL*UMAP(1+3)
      UMAP(7) = PXL
      UMAP(7+1) = PYL
      UMAP(9) = PXU
      UMAP(9+1) = PYL
      UMAP(11) = PXU
      UMAP(11+1) = PYU
      UMAP(13) = PXL
      UMAP(13+1) = PYU
      RETURN

      END
