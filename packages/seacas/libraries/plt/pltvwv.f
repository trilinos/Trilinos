C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTVWV(PLL,PUR,N,MASK,PX,PY,QX,QY)
      REAL PLL(2)
      REAL PUR(2)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL QX(*)
      REAL QY(*)
      include 'izbit.inc'

      PUR1 = PUR(1) + .0001
      PUR2 = PUR(2) + .0001
      PLL1 = PLL(1) - .0001
      PLL2 = PLL(2) - .0001
      DX = PUR1 - PLL1
      DY = PUR2 - PLL2
      J = 0
      KM = 0
   10 CONTINUE
      IF ((J.LT.N)) THEN
         JN = MIN(N-J,32)
         J1 = J
         J = J + JN
         KM = KM + 1
         M = MASK(KM)
         IF (M.NE.0) THEN

            DO 20 K = 1,JN
               JB = IZBIT(K)
               IF (IAND(M,JB).NE.0) THEN

                  X1 = PX(K+J1)
                  Y1 = PY(K+J1)
                  X2 = QX(K+J1)
                  Y2 = QY(K+J1)
                  FP = X1 - PLL1
                  FQ = X2 - PLL1
                  IF (FP.LT.0. .AND. FQ.LT.0.) THEN
                     M = IAND(M,NOT(JB))

                  ELSE IF (FP.GT.DX .AND. FQ.GT.DX) THEN
                     M = IAND(M,NOT(JB))

                  ELSE

                     DF = FQ - FP
                     IF (DF.GT.0.) THEN
                        TN = (Y2-Y1)/DF
                        IF (FP.LT.0.) THEN
                           X1 = PLL1
                           Y1 = Y1 - FP*TN
                        END IF

                        IF (FQ.GT.DX) THEN
                           X2 = PUR1
                           Y2 = Y2 + (DX-FQ)*TN
                        END IF

                     ELSE IF (DF.LT.0.) THEN
                        TN = (Y2-Y1)/DF
                        IF (FQ.LT.0.) THEN
                           X2 = PLL1
                           Y2 = Y2 - FQ*TN
                        END IF

                        IF (FP.GT.DX) THEN
                           X1 = PUR1
                           Y1 = Y1 + (DX-FP)*TN
                        END IF

                     END IF

                     FP = Y1 - PLL2
                     FQ = Y2 - PLL2
                     IF (FP.LT.0. .AND. FQ.LT.0.) THEN
                        M = IAND(M,NOT(JB))

                     ELSE IF (FP.GT.DY .AND. FQ.GT.DY) THEN
                        M = IAND(M,NOT(JB))

                     ELSE

                        DF = FQ - FP
                        IF (DF.GT.0.) THEN
                           TN = (X2-X1)/DF
                           IF (FP.LT.0.) THEN
                              Y1 = PLL2
                              X1 = X1 - FP*TN
                           END IF

                           IF (FQ.GT.DY) THEN
                              Y2 = PUR2
                              X2 = X2 + (DY-FQ)*TN
                           END IF

                        ELSE IF (DF.LT.0.) THEN
                           TN = (X2-X1)/DF
                           IF (FQ.LT.0.) THEN
                              Y2 = PLL2
                              X2 = X2 - FQ*TN
                           END IF

                           IF (FP.GT.DY) THEN
                              Y1 = PUR2
                              X1 = X1 + (DY-FP)*TN
                           END IF

                        END IF

                        PX(K+J1) = X1
                        PY(K+J1) = Y1
                        QX(K+J1) = X2
                        QY(K+J1) = Y2
                     END IF

                  END IF

               END IF

   20       CONTINUE
            MASK(KM) = M
         END IF

         GO TO 10

      END IF

      END
