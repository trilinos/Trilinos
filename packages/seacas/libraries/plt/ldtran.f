C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LDTRAN(X,Y,Z,MAT)
      REAL MAT(4,4)

      CALL MXIDEN(4,MAT)
      MAT(4,1) = X
      MAT(4,2) = Y
      MAT(4,3) = Z
      RETURN

      END
