C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTRST
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

      IF ((DEVCAP(15)+DEVCAP(16))/2..LE.800.) THEN
         CALL PLTSTT(12,2.)

      ELSE
         CALL PLTSTT(12,3.)
      END IF

      BUFF = 1./58.
      CALL VDSTCS(BUFF)
      TEXTP(35) = DEFOUT(6)
      TEXTP(36) = DEFOUT(7)
      TEXTP(1) = .015
      TEXTP(2) = 0.
      TEXTP(3) = 0.
      DO 2340 I = 4,11
         TEXTP(I) = 0.
 2340 CONTINUE
      TEXTP(20) = 0.
      TEXTP(21) = 0.
      TEXTP(22) = 1.
      TEXTP(23) = 0.
      TEXTP(24) = TEXTP(22)
      TEXTP(25) = .75
      TEXTP(26) = 0.
      TEXTP(27) = TEXTP(25)
      TEXTP(30) = 1.5
      TEXTP(31) = 1.
      TEXTP(32) = .7
      TEXTP(33) = .8
      TEXTP(34) = .5
      TEXTP(37) = 160.
      CALL PLTITM
      RETURN

      END
