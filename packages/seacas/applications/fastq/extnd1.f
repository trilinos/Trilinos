C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EXTND1 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
C***********************************************************************

C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), ANGLE (MXND)
      DIMENSION X(1), Y(1)

      CANG = (ANGLE (N2) * .5)
      ANG = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2)) - CANG
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N3) - YN (N2)) **2 +  (XN (N3) - XN (N2)) **2)
      DIST = (DIST1 + DIST2) * .5
      IF (CANG .EQ. 0.) THEN
         ADIST = DIST
      ELSE
         ADIST = DIST / SIN (CANG)
      ENDIF

      X(1) = ADIST * COS (ANG) + XN (N2)
      Y(1) = ADIST * SIN (ANG) + YN (N2)

      RETURN

      END
