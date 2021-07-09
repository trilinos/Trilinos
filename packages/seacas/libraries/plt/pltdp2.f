C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTDP2(MAP,N,PX,PY)
      REAL MAP(*),PX(*),PY(*)
      REAL XWORK(32),YWORK(32)
      INTEGER MASK(1)

      J = 0
 2360 IF (.NOT. (J.LT.N)) GO TO 2370
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN

      do 2400 i=1, jn
        XWORK(I) = MAP(1)*PX(J1+I) + MAP(3)*PY(J1+I) + MAP(5)
        YWORK(I) = MAP(2)*PX(J1+I) + MAP(4)*PY(J1+I) + MAP(6)
 2400 CONTINUE

      MASK(1) = -1
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(7),MAP(9))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(9),MAP(11))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(11),MAP(13))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(13),MAP(7))
      CALL PLTPTM(JN,MASK,XWORK,YWORK)
      GO TO 2360

 2370 CONTINUE
      RETURN

      END
