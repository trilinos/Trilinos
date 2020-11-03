C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EXSIDE (MXND, XN, YN, LNODES, ANGLE, N1, XNEW, YNEW)
C***********************************************************************

C  SUBROUTINE EXSIDE = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A SIDE NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)

      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)

      ANG = ATAN2 (YN (N2) - YN (N1),  XN (N2) - XN (N1)) +
     &   (ANGLE (N1) * .5)

      DIST1 = SQRT ( (YN (N2) - YN (N1) ) **2 +
     &   ( XN (N2) - XN (N1) ) **2)
      DIST2 = SQRT ( (YN (N0) - YN (N1) ) **2 +
     &   ( XN (N0) - XN (N1) ) **2)

      DIST =  (DIST1 + DIST2) * .5
      XNEW =  (DIST * COS (ANG) ) + XN (N1)
      YNEW =  (DIST * SIN (ANG) ) + YN (N1)

      RETURN

      END
