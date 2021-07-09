C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      LOGICAL FUNCTION DISCTP (ANGLE)
C***********************************************************************

C  FUNCTION DISCTP = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ANGLE IS
C                    WITHIN THE CURRENT DEFINITION OF A
C                    DISSECTION VERTEX

C***********************************************************************

      DATA EPS /.31/

      PI = ATAN2(0.0, -1.0)
      IF (ANGLE .GT. (PI + EPS)) THEN
         DISCTP=.TRUE.
      ELSE
         DISCTP=.FALSE.
      ENDIF
      RETURN
      END
