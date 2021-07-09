C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTZCP(ZNEAR,ZFAR,N,MASK,PZ)
      INTEGER N
      INTEGER MASK(*)
      REAL PZ(*)
      include 'izbit.inc'

      J = 0
      KM = 0
 2380 IF (.NOT. (J.LT.N)) GO TO 2390
      JN = MIN(N-J,32)
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2380

      END IF

      DO 2400 I = 1,JN
         JB = IZBIT(I)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2400

         END IF

         IF ((PZ(I).LT.ZNEAR) .OR. (PZ(I).GT.ZFAR)) THEN
            M = IAND(M,NOT(JB))
            GO TO 2400

         END IF

         MASK(KM) = M
 2400 CONTINUE
      GO TO 2380

 2390 CONTINUE
      RETURN

      END
