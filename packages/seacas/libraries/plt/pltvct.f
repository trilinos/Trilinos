C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTVCT(N,XX0,YY0,XX1,YY1)
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
      DIMENSION XX0(*),YY0(*),XX1(*),YY1(*)

      IF (VECTP(1).LE.0. .OR. VECTP(2).LE.0.) THEN
         XCUR = XX1(N)
         YCUR = YY1(N)
         RETURN

      END IF

      DO 2160 I = 1,N
         IF (XX0(I).NE.XCUR .OR. YY0(I).NE.YCUR) THEN
            CALL PLTMOV(XX0(I),YY0(I))
         END IF

         CALL PLTDRW(XX1(I),YY1(I))
         XCUR = XX1(I)
         YCUR = YY1(I)
 2160 CONTINUE
      RETURN

      END
