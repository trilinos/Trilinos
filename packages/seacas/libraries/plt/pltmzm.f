C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmzm.f,v 1.1 1993/07/16 16:48:58 gdsjaar Exp $
C $Log: pltmzm.f,v $
C Revision 1.1  1993/07/16 16:48:58  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTMZM(FACT,UMAP)
      REAL FACT
      REAL UMAP(*)

      UMAP(30) = FACT
      RETURN

      END
