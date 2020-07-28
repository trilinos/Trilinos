C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLIG(X,Y)
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
      INTEGER MASK(1)

      IF (MAPP(10).EQ.1.) THEN
         CALL PLTP2D(X,Y,XP,YP)
         CALL PLTP2D(XCUR,YCUR,XC,YC)

      ELSE
         XP = X
         YP = Y
         XC = XCUR
         YC = YCUR
      END IF

      MASK(1) = -1
      CALL PLTVWV(MAPP(6),MAPP(8),1,MASK,XC,YC,XP,YP)
      IF (MASK(1).EQ.-1) THEN
         CALL VDMOVA(XC,YC)
         CALL VDLINA(XP,YP)
      END IF

      RETURN

      END
