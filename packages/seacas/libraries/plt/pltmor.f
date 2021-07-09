C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMOR(X,Y,Z,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL X,Y,Z
      REAL UMAP(*)

      DX = X - XC
      DY = Y - YC
      DZ = Z - ZC
      CALL PLTMMV(DX,DY,DZ,UMAP)
      RETURN

      END
