C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltix2.f,v 1.1 1993/07/16 16:48:32 gdsjaar Exp $
C $Log: pltix2.f,v $
C Revision 1.1  1993/07/16 16:48:32  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTIX2(UMAP)
      REAL UMAP(*)

      UMAP(1+3) = -UMAP(1+3)
      RETURN

      END
