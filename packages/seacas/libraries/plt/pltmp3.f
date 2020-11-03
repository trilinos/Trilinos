C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMP3(UMAP,N,MASK,PX,PY,PZ,QX,QY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),PZ(*),QX(*),QY(*)
      DIMENSION Q1(3),V1(3),Q2(3),V2(3),TPX(32),TPY(32)
      include 'izbit.inc'

      DO 2060 L = 1,3
         V1(L) = UMAP(18+L-1) + UMAP(15)*UMAP(27+L-1)
         Q1(L) = UMAP(L+27-1)
         V2(L) = UMAP(18+L-1) + UMAP(16)*UMAP(27+L-1)
         Q2(L) = -UMAP(L+27-1)
 2060 CONTINUE
      J = 0
      KM = 0
 2080 IF (.NOT. (J.LT.N)) GO TO 2090
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V1,Q1)
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V2,Q2)
      IF (UMAP(17).EQ.1.) THEN
         M = MASK(KM)
         DO 2100 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,M).NE.0) THEN
               PMS = (PX(K+J1-1)-UMAP(18))*UMAP(27) +
     *               (PY(K+J1-1)-UMAP(19))*UMAP(28) +
     *               (PZ(K+J1-1)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TPX(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(21)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(22)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(23))
               TPY(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(24)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(25)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(26))
            END IF

 2100    CONTINUE
      END IF

      IF (UMAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMP2(UMAP(1),JN,MASK(KM),TPX,TPY,QX(J1),QY(J1))
      GO TO 2080

 2090 CONTINUE
      RETURN

      END
