C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCV2(N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY,C1,C2)
      DIMENSION MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),QQX(*),
     *          QQY(*),C1(*),C2(*)
      include 'izbit.inc'

      CX = C1(1)
      CY = C1(2)
      DX = C2(1) - CX
      DY = C2(2) - CY
      J = 0
      KM = 0
 2100 IF (.NOT. (J.LT.N)) GO TO 2110
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2100

      END IF

      DO 2120 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(JB,M).EQ.0) THEN
            GO TO 2120

         END IF

         X1 = PX(J1+K)
         Y1 = PY(J1+K)
         X2 = QX(J1+K)
         Y2 = QY(J1+K)
         FP = (Y1-CY)*DX - (X1-CX)*DY
         FQ = (Y2-CY)*DX - (X2-CX)*DY
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))

         ELSE IF (FP.LT.0.) THEN
            XL = FQ/ (FQ-FP)
            X1 = X2 + XL* (X1-X2)
            Y1 = Y2 + XL* (Y1-Y2)

         ELSE IF (FQ.LT.0.) THEN
            XL = FP/ (FP-FQ)
            X2 = X1 + XL* (X2-X1)
            Y2 = Y1 + XL* (Y2-Y1)
         END IF

         PPX(K+J1) = X1
         PPY(K+J1) = Y1
         QQX(K+J1) = X2
         QQY(K+J1) = Y2
 2120 CONTINUE
      MASK(KM) = M
      GO TO 2100

 2110 CONTINUE
      RETURN

      END
