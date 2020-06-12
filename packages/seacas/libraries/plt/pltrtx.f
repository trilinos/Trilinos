C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltrtx.f,v 1.1 1993/07/16 16:49:24 gdsjaar Exp $
C $Log: pltrtx.f,v $
C Revision 1.1  1993/07/16 16:49:24  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTRTX(VAL,UMAP)
      REAL UMAP(*)
      REAL VAL
      REAL A(9),B(9)

      S = SIN(VAL)
      C = COS(VAL)
      DO 2280 I = 1,9
         A(I) = UMAP(21-1+I)
         B(I) = 0.
 2280 CONTINUE
      B(1) = 1.
      B(5) = C
      B(9) = C
      B(6) = S
      B(8) = -S
      CALL PLTROT(B,A,UMAP(21))
      RETURN

      END
