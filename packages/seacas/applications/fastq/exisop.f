C    Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EXISOP (MXND, XN, YN, LNODES, ANGLE, N1, XNEW, YNEW)
C***********************************************************************

C  SUBROUTINE EXISOP = CALCULATES A POSITION TO MAKE A PARALLELEPIPED

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)

      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)

      XNEW = XN (N0) + XN (N2) - XN (N1)
      YNEW = YN (N0) + YN (N2) - YN (N1)

      RETURN

      END
