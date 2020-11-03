C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMV3(UMAP,N,MASK,UX,UY,UZ,VX,VY,VZ,PX,PY,QX,QY)
      DIMENSION UMAP(*),MASK(*),UX(*),UY(*),UZ(*),VX(*),VY(*),VZ(*),
     *          PX(*),PY(*),QX(*),QY(*)
      DIMENSION TUX(32),TUY(32),TUZ(32),TVX(32),TVY(32),TVZ(32),
     *          TTUX(32),TTUY(32),TTUZ(32),TTVX(32),TTVY(32),TTVZ(32),
     *          V1(3),Q1(3),V2(3),Q2(3)
      include 'izbit.inc'

      DO 2160 L = 1,3
         V1(L) = UMAP(18+L-1) + UMAP(15)*UMAP(27+L-1)
         Q1(L) = UMAP(27+L-1)
         V2(L) = UMAP(18+L-1) + UMAP(16)*UMAP(27+L-1)
         Q2(L) = -UMAP(27+L-1)
 2160 CONTINUE
      J = 0
      KM = 0
 2180 IF (.NOT. (J.LT.N)) GO TO 2190
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      CALL PLTCV3(JN,MASK(KM),UX(J1),UY(J1),UZ(J1),VX(J1),VY(J1),VZ(J1),
     *            TUX,TUY,TUZ,TVX,TVY,TVZ,V1,Q1)
      CALL PLTCV3(JN,MASK(KM),TUX,TUY,TUZ,TVX,TVY,TVZ,TTUX,TTUY,TTUZ,
     *            TTVX,TTVY,TTVZ,V2,Q2)
      IF (UMAP(17).EQ.1.) THEN
         DO 2200 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK(KM)).NE.0) THEN
               PMS = (TTUX(K)-UMAP(18))*UMAP(27) +
     *               (TTUY(K)-UMAP(19))*UMAP(28) +
     *               (TTUZ(K)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TUX(K) = R* ((TTUX(K)-UMAP(18))*UMAP(21)+
     *                  (TTUY(K)-UMAP(19))*UMAP(22)+
     *                  (TTUZ(K)-UMAP(20))*UMAP(23))
               TUY(K) = R* ((TTUX(K)-UMAP(18))*UMAP(24)+
     *                  (TTUY(K)-UMAP(19))*UMAP(25)+
     *                  (TTUZ(K)-UMAP(20))*UMAP(26))
               PMS = (TTVX(K)-UMAP(18))*UMAP(27) +
     *               (TTVY(K)-UMAP(19))*UMAP(28) +
     *               (TTVZ(K)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TVX(K) = R* ((TTVX(K)-UMAP(18))*UMAP(21)+
     *                  (TTVY(K)-UMAP(19))*UMAP(22)+
     *                  (TTVZ(K)-UMAP(20))*UMAP(23))
               TVY(K) = R* ((TTVX(K)-UMAP(18))*UMAP(24)+
     *                  (TTVY(K)-UMAP(19))*UMAP(25)+
     *                  (TTVZ(K)-UMAP(20))*UMAP(26))
            END IF

 2200    CONTINUE

      ELSE IF (UMAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMV2(UMAP,JN,MASK(KM),TUX,TUY,TVX,TVY,PX(J1),PY(J1),QX(J1),
     *            QY(J1))
      GO TO 2180

 2190 CONTINUE
      RETURN

      END
