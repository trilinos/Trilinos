C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: ldscal.f,v 1.1 1993/07/16 16:46:35 gdsjaar Exp $
C $Log: ldscal.f,v $
C Revision 1.1  1993/07/16 16:46:35  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE LDSCAL(X,Y,Z,MAT)
      REAL MAT(4,4)

      CALL MXIDEN(4,MAT)
      MAT(1,1) = X
      MAT(2,2) = Y
      MAT(3,3) = Z
      RETURN

      END
