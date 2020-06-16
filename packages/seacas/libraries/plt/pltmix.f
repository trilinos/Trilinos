C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmix.f,v 1.1 1993/07/16 16:48:48 gdsjaar Exp $
C $Log: pltmix.f,v $
C Revision 1.1  1993/07/16 16:48:48  gdsjaar
C Changed plt to library rather than single source file.
C
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
