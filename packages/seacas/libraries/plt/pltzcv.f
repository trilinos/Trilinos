C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTZCV(ZNEAR,ZFAR,N,MASK,PX,PY,PZ,QX,QY,QZ)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL PZ(*)
      REAL QX(*)
      REAL QY(*)
      REAL QZ(*)
      include 'izbit.inc'

      J = 0
      KM = 0
      DZ = ZFAR - ZNEAR
 2420 IF (.NOT. (J.LT.N)) GO TO 2430
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2420

      END IF

      DO 2440 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2440

         END IF

         X1 = PX(K+J1)
         Y1 = PY(K+J1)
         Z1 = PZ(K+J1)
         X2 = QX(K+J1)
         Y2 = QY(K+J1)
         Z2 = QZ(K+J1)
         FP = Z1 - ZNEAR
         FQ = Z2 - ZNEAR
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         IF (FP.GT.DZ .AND. FQ.GT.DZ) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         DF = FQ - FP
         IF (DF.GT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FP.LT.0.) THEN
               Z1 = ZNEAR
               X1 = X1 - FP*TN
               Y1 = Y1 - FP*SN
            END IF

            IF (FQ.GT.DZ) THEN
               Z2 = ZFAR
               X2 = X2 + (DZ-FQ)*TN
               Y2 = Y2 + (DZ-FQ)*SN
            END IF

         ELSE IF (DF.LT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FQ.LT.0.) THEN
               Z2 = ZNEAR
               X2 = X2 - FQ*TN
               Y2 = Y2 - FQ*SN
            END IF

            IF (FP.GT.DZ) THEN
               Z1 = ZFAR
               X1 = X1 + (DZ-FP)*TN
               Y1 = Y1 + (DZ-FP)*SN
            END IF

         END IF

         PX(K+J1) = X1
         PY(K+J1) = Y1
         PZ(K+J1) = Z1
         QX(K+J1) = X2
         QY(K+J1) = Y2
         QZ(K+J1) = Z2
 2440 CONTINUE
      MASK(KM) = M
      GO TO 2420

 2430 CONTINUE
      RETURN

      END
