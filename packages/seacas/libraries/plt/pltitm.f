C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTITM
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
      DATA DPR/57.29577951/

      COSD(X) = COS(X/DPR)
      SIND(X) = SIN(X/DPR)
      TEXTP(12) = COSD(TEXTP(2)+TEXTP(3))
      TEXTP(13) = SIND(TEXTP(2)+TEXTP(3))
      TEXTP(14) = TEXTP(12)*TEXTP(1)
      TEXTP(15) = TEXTP(13)*TEXTP(1)
      TEXTP(16) = -TEXTP(13)*TEXTP(1)
      TEXTP(17) = TEXTP(12)*TEXTP(1)
      TEXTP(28) = COSD(TEXTP(2))
      TEXTP(29) = SIND(TEXTP(2))
      RETURN

      END
