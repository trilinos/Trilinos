C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE EXCORN (MXND, XN, YN, LNODES, ANGLE, N0, N1, N2, XNEW,
     &   YNEW)
C***********************************************************************

C  SUBROUTINE EXCORN = CALCULATES A POSITION AN AVERAGE LENGTH AWAY
C                      FROM A CORNER NODE

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)

      LOGICAL SIDEP

C      XNEW = XN (N0) + XN (N2) - XN (N1)
C      YNEW = YN (N0) + YN (N2) - YN (N1)

      PID2 = 0.5 * ATAN2(0.0, -1.0)

C      ANG2 = ATAN2 (YN (N1)-YN (N0), XN (N1)-XN (N0))+PID2
      BANG1 = ATAN2 (YN (N1) - YN (N2), XN (N1) - XN (N2))
      BANG2 = ATAN2 (YN (N1) - YN (N0), XN (N1) - XN (N0))
      IF (SIDEP (ANGLE (N2)))THEN
         ANG1 = BANG1 - (ANGLE (N2) * .5)
      ELSE
         ANG1 = BANG1 - PID2
      ENDIF
      IF (SIDEP (ANGLE (N0)))THEN
         ANG2 = BANG2 + (ANGLE (N0) * .5)
      ELSE
         ANG2 = BANG2 + PID2
      ENDIF
      DIST1 = SQRT ((YN (N2)-YN (N1)) **2 +  (XN (N2)-XN (N1)) **2)
      DIST2 = SQRT ((YN (N0)-YN (N1)) **2 +  (XN (N0)-XN (N1)) **2)
      DIST = (DIST1 + DIST2) * .5
      XNEW = ( (DIST * COS (ANG1) + XN (N2)) +
     &   (DIST * COS (ANG2) + XN (N0)) ) * .5
      YNEW = ( (DIST * SIN (ANG1) + YN (N2)) +
     &   (DIST * SIN (ANG2) + YN (N0)) ) * .5
C      XNEW = DIST * COS (ANG1) + XN (N2)
C      YNEW = DIST * SIN (ANG1) + YN (N2)

      RETURN

      END
