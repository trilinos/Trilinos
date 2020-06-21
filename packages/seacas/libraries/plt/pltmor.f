C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltmor.f,v 1.1 1993/07/16 16:48:53 gdsjaar Exp $
C $Log: pltmor.f,v $
C Revision 1.1  1993/07/16 16:48:53  gdsjaar
C Changed plt to library rather than single source file.
C
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
