C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTVWP(PLL,PUR,N,MASK,PX,PY)
      REAL PLL(2),PUR(2)
      INTEGER MASK(*)
      REAL PX(*),PY(*)
      include 'izbit.inc'

      PUR1 = PUR(1) + .0001
      PUR2 = PUR(2) + .0001
      PLL1 = PLL(1) - .0001
      PLL2 = PLL(2) - .0001
      J = 0
      KM = 0
 2300 IF (.NOT. (J.LT.N)) GO TO 2310
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      IF (MASK(KM).EQ.0) THEN
         GO TO 2300

      END IF

      DO 2320 K = 1,JN
         IF (IAND(MASK(KM),IZBIT(K)).EQ.0) THEN
            GO TO 2320

         END IF

         IF (PX(K+J1).LT.PLL1 .OR. PX(K+J1).GT.PUR1) THEN
            MASK(KM) = IAND(MASK(KM),NOT(IZBIT(K)))
            GO TO 2320

         END IF

         IF (PY(K+J1).LT.PLL2 .OR. PY(K+J1).GT.PUR2) THEN
            MASK(KM) = IAND(MASK(KM),NOT(IZBIT(K)))
            GO TO 2320

         END IF

 2320 CONTINUE
      GO TO 2300

 2310 CONTINUE
      RETURN

      END
