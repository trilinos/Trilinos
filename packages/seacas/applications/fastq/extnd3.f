C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EXTND3 (MXND, XN, YN, ANGLE, N1, N2, N3, X, Y, DIST)
C***********************************************************************

C  SUBROUTINE EXTND3 = CALCULATES TWO POSITIONS AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE AND ONE AT 1/3 ANGLE INTERVALS

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND)
      DIMENSION X(3), Y(3)

      ANG = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2))
      CANG = ANGLE / 3.
      ANG1 = ANG - CANG
      ANG2 = ANG - (ANGLE * .5)
      ANG3 = ANG - 2. * CANG
      DIST1 = SQRT ((YN (N2) - YN (N1)) **2 +  (XN (N2) - XN (N1)) **2)
      DIST2 = SQRT ((YN (N3) - YN (N2)) **2 +  (XN (N3) - XN (N2)) **2)
      DIST = (DIST1 + DIST2) * .5
      IF (CANG. EQ. 0.) THEN
         ADIST = DIST
      ELSE
         ADIST = DIST / SIN (CANG)
      ENDIF

      X(1) = ADIST * COS (ANG1) + XN (N2)
      Y(1) = ADIST * SIN (ANG1) + YN (N2)
      X(2) = 1.4142 * ADIST * COS (ANG2) + XN (N2)
      Y(2) = 1.4142 * ADIST * SIN (ANG2) + YN (N2)
      X(3) = ADIST * COS (ANG3) + XN (N2)
      Y(3) = ADIST * SIN (ANG3) + YN (N2)

      RETURN

      END
