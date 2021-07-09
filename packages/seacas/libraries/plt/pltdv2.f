C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTDV2(MAP,N,PX,PY,QX,QY)
      REAL MAP(*),PX(*),PY(*),QX(*),QY(*)
      DIMENSION PPX(32),PPY(32),QQX(32),QQY(32)
      INTEGER MASK(1)

      J = 0
 2470 IF (J.LT.N) THEN
         JN = MIN(N-J,32)
         J1 = J
         J = J + JN
         DO 2490 I = 1,JN
            PPX(I) = MAP(1)*PX(I+J1) + MAP(3)*PY(I+J1) + MAP(5)
            QQX(I) = MAP(1)*QX(I+J1) + MAP(3)*QY(I+J1) + MAP(5)
            PPY(I) = MAP(2)*PX(I+J1) + MAP(4)*PY(I+J1) + MAP(6)
            QQY(I) = MAP(2)*QX(I+J1) + MAP(4)*QY(I+J1) + MAP(6)
 2490    CONTINUE

         MASK(1) = -1
         CALL PLTVWV(MAP(7),MAP(11),JN,MASK,PPX,PPY,QQX,QQY)
         CALL PLTVCM(JN,MASK,PPX,PPY,QQX,QQY)
         GO TO 2470

      END IF
      RETURN

      END
