C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCV3(N,MASK,PX,PY,PZ,QX,QY,QZ,PPX,PPY,PPZ,QQX,QQY,
     *                  QQZ,V,Q)
      DIMENSION MASK(*),PX(*),PY(*),PZ(*),QX(*),QY(*),QZ(*),PPX(*),
     *          PPY(*),PPZ(*),QQX(*),QQY(*),QQZ(*),V(*),Q(*)
      include 'izbit.inc'

      CALL CPUMVU(PX,PPX,N)
      CALL CPUMVU(PY,PPY,N)
      CALL CPUMVU(PZ,PPZ,N)
      CALL CPUMVU(QX,QQX,N)
      CALL CPUMVU(QY,QQY,N)
      CALL CPUMVU(QZ,QQZ,N)
      J = 0
      KM = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J
      KM = 1 + KM
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2140

      END IF

      DO 2160 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(JB,M).EQ.0) THEN
            GO TO 2160

         END IF

         FP = (PPX(J1+K)-V(1))*Q(1) + (PPY(J1+K)-V(2))*Q(2) +
     *        (PPZ(J1+K)-V(3))*Q(3)
         FQ = (QQX(J1+K)-V(1))*Q(1) + (QQY(J1+K)-V(2))*Q(2) +
     *        (QQZ(J1+K)-V(3))*Q(3)
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))

         ELSE IF (FP.LT.0.) THEN
            XL = FP/ (FP-FQ)
            PPX(J1+K) = PPX(J1+K) + XL* (QQX(J1+K)-PPX(J1+K))
            PPY(J1+K) = PPY(J1+K) + XL* (QQY(J1+K)-PPY(J1+K))
            PPZ(J1+K) = PPZ(J1+K) + XL* (QQZ(J1+K)-PPZ(J1+K))

         ELSE IF (FQ.LT.0.) THEN
            XL = FQ/ (FQ-FP)
            QQX(J1+K) = QQX(J1+K) + XL* (PPX(J1+K)-QQX(J1+K))
            QQY(J1+K) = QQY(J1+K) + XL* (PPY(J1+K)-QQY(J1+K))
            QQZ(J1+K) = QQZ(J1+K) + XL* (PPZ(J1+K)-QQZ(J1+K))
         END IF

 2160 CONTINUE
      MASK(KM) = M
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
