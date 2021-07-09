C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCP3(N,MASK,PX,PY,PZ,V,Q)
      DIMENSION MASK(*),PX(*),PY(*),PZ(*),V(*),Q(*)
      include 'izbit.inc'

      J = 0
      KM = 0
 2060 IF (.NOT. (J.LT.N)) GO TO 2070
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2060

      END IF

      DO 2080 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).NE.0) THEN
            FP = (PX(J1+K)-V(1))*Q(1) + (PY(J1+K)-V(2))*Q(2) +
     *           (PZ(J1+K)-V(3))*Q(3)
            IF (FP.LT.0.) THEN
               M = IAND(M,NOT(JB))
            END IF

         END IF

 2080 CONTINUE
      MASK(KM) = M
      GO TO 2060

 2070 CONTINUE
      RETURN

      END
