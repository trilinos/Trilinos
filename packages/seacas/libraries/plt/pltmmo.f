C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmmo.f,v 1.1 1993/07/16 16:48:51 gdsjaar Exp $
C $Log: pltmmo.f,v $
C Revision 1.1  1993/07/16 16:48:51  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTMMO(FACT,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      R = (UMAP(18)-XC)**2 + (UMAP(18+1)-YC)**2 + (UMAP(18+2)-ZC)**2
      R = SQRT(R)*FACT
      UMAP(18) = XC - R*UMAP(27)
      UMAP(18+1) = YC - R*UMAP(27+1)
      UMAP(18+2) = ZC - R*UMAP(27+2)
      RETURN

      END
