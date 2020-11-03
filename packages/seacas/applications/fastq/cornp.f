C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      LOGICAL FUNCTION CORNP (ANGLE)
C***********************************************************************

C  FUNCTION CORNP = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ANGLE IS
C                   WITHIN THE CURRENT DEFINITION OF A CORNER

C***********************************************************************

      DATA EPS /.62/

      PI = ATAN2(0.0, -1.0)
      IF (ANGLE .LT. ( PI - EPS)) THEN
         CORNP=.TRUE.
      ELSE
         CORNP=.FALSE.
      ENDIF

      RETURN

      END
