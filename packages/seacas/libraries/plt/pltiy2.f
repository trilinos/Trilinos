C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltiy2.f,v 1.1 1993/07/16 16:48:33 gdsjaar Exp $
C $Log: pltiy2.f,v $
C Revision 1.1  1993/07/16 16:48:33  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTIY2(UMAP)
      REAL UMAP(*)

      UMAP(1) = -UMAP(1)
      RETURN

      END
