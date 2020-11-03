C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTDV3(MAP,N,UX,UY,UZ,VX,VY,VZ)
      REAL MAP(*),UX(*),UY(*),UZ(*),VX(*),VY(*),VZ(*)
      DIMENSION TUX(32),TUY(32),TUZ(32),TVX(32),TVY(32),TVZ(32),
     *          TTUX(32),TTUY(32),TTUZ(32),TTVX(32),TTVY(32),TTVZ(32),
     *          V1(3),Q1(3),V2(3),Q2(3)
      INTEGER MASK(1)
      include 'izbit.inc'

      DO 2510 L = 1,3
         V1(L) = MAP(18+L-1) + MAP(15)*MAP(27+L-1)
         Q1(L) = MAP(27+L-1)
         V2(L) = MAP(18+L-1) + MAP(16)*MAP(27+L-1)
         Q2(L) = -MAP(27+L-1)
 2510 CONTINUE
      J = 0
 2530 IF (.NOT. (J.LT.N)) GO TO 2540
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      MASK(1) = -1
      CALL PLTCV3(JN,MASK,UX(J1),UY(J1),UZ(J1),VX(J1),VY(J1),VZ(J1),TUX,
     *            TUY,TUZ,TVX,TVY,TVZ,V1,Q1)
      CALL PLTCV3(JN,MASK,TUX,TUY,TUZ,TVX,TVY,TVZ,TTUX,TTUY,TTUZ,TTVX,
     *            TTVY,TTVZ,V2,Q2)
      IF (MAP(17).EQ.1.) THEN
         DO 2550 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK(1)).NE.0) THEN
               PMS = (TTUX(K)-MAP(18))*MAP(27) +
     *               (TTUY(K)-MAP(19))*MAP(28) +
     *               (TTUZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TUX(K) = R* ((TTUX(K)-MAP(18))*MAP(21)+
     *                  (TTUY(K)-MAP(19))*MAP(22)+
     *                  (TTUZ(K)-MAP(20))*MAP(23))
               TUY(K) = R* ((TTUX(K)-MAP(18))*MAP(24)+
     *                  (TTUY(K)-MAP(19))*MAP(25)+
     *                  (TTUZ(K)-MAP(20))*MAP(26))
               PMS = (TTVX(K)-MAP(18))*MAP(27) +
     *               (TTVY(K)-MAP(19))*MAP(28) +
     *               (TTVZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TVX(K) = R* ((TTVX(K)-MAP(18))*MAP(21)+
     *                  (TTVY(K)-MAP(19))*MAP(22)+
     *                  (TTVZ(K)-MAP(20))*MAP(23))
               TVY(K) = R* ((TTVX(K)-MAP(18))*MAP(24)+
     *                  (TTVY(K)-MAP(19))*MAP(25)+
     *                  (TTVZ(K)-MAP(20))*MAP(26))
            END IF

 2550    CONTINUE

      ELSE IF (MAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMV2(MAP,JN,MASK,TUX,TUY,TVX,TVY,TTUX,TTUY,TTVX,TTVY)
      CALL PLTVCM(JN,MASK,TTUX,TTUY,TTVX,TTVY)
      GO TO 2530

 2540 CONTINUE
      RETURN

      END
