C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTVCM(N,MASK,XX0,YY0,XX1,YY1)
      REAL DEVCAP(23)
      REAL DEFOUT(7)
      COMMON /STATUS/DEVCAP,DEFOUT
      REAL DEVP(5)
      COMMON /DEVICE/DEVP
      REAL COLP(3)
      REAL PALETT(3,16)
      COMMON /COLOR/COLP,PALETT
      REAL TEXTP(40)
      COMMON /TEXT/TEXTP
      REAL VECTP(5)
      REAL XCUR
      REAL YCUR
      COMMON /VECTRC/VECTP,XCUR,YCUR
      INTEGER IDEX(200,2)
      INTEGER NVECT(200,2)
      REAL XSIZE(200,2)
      REAL YSIZE(200,2)
      REAL X0(2300,2)
      REAL Y0(2300,2)
      REAL X1(2300,2)
      REAL Y1(2300,2)
      COMMON /FONT/IDEX,NVECT,XSIZE,YSIZE,X0,Y0,X1,Y1
      REAL GRAPHP(100)
      COMMON /GRAPH/GRAPHP
      COMMON /MAPPAR/MAPP(11)
      REAL MAPP
      COMMON /STORAG/MEMORY(1000)
      DIMENSION XX0(*),YY0(*),XX1(*),YY1(*),MASK(*)
      include 'izbit.inc'

      IF (VECTP(1).LE.0. .OR. VECTP(2).LE.0.) THEN
         XCUR = XX1(N)
         YCUR = YY1(N)
         RETURN

      END IF

      J = 0
      KM = 0
 2180 IF (.NOT. (J.LT.N)) GO TO 2190
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2180

      END IF

      DO 2200 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2200

         END IF

         IF (XX0(K+J1).NE.XCUR .OR. YY0(K+J1).NE.YCUR) THEN
            CALL PLTMOV(XX0(K+J1),YY0(K+J1))
         END IF

         CALL PLTDRW(XX1(K+J1),YY1(K+J1))
         XCUR = XX1(K+J1)
         YCUR = YY1(K+J1)
 2200 CONTINUE
      GO TO 2180

 2190 CONTINUE
      RETURN

      END
