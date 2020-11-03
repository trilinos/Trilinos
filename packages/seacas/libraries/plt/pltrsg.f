C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTRSG
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

      GRAPHP(7) = 1.
      GRAPHP(8) = 0.
      GRAPHP(9) = 0.
      GRAPHP(10) = 1.
      GRAPHP(11) = 0.
      GRAPHP(12) = 0.
      GRAPHP(13) = 0.
      GRAPHP(14) = 0.
      GRAPHP(15) = 1.
      GRAPHP(16) = 0.
      GRAPHP(17) = 1.
      GRAPHP(18) = .75
      GRAPHP(19) = 0.
      GRAPHP(20) = .75
      GRAPHP(1) = .15
      GRAPHP(2) = .075
      GRAPHP(3) = .75
      GRAPHP(4) = .6
      GRAPHP(5) = 1.
      GRAPHP(6) = 0.
      GRAPHP(23) = 1.
      GRAPHP(21) = 1.
      GRAPHP(22) = 1.
      GRAPHP(47) = 1. + 4
      GRAPHP(32) = 1.
      GRAPHP(35) = 0.
      GRAPHP(38) = DEFOUT(1)
      GRAPHP(36) = DEFOUT(1)
      GRAPHP(37) = DEFOUT(1)
      GRAPHP(39) = DEFOUT(1)
      GRAPHP(40) = 3.
      GRAPHP(41) = 1.
      GRAPHP(42) = 1.
      GRAPHP(44) = 5.
      GRAPHP(45) = 5.
      GRAPHP(88) = 5.
      GRAPHP(89) = 5.
      GRAPHP(46) = 5.
      GRAPHP(48) = 0.
      GRAPHP(49) = 0.
      GRAPHP(62) = 160.
      GRAPHP(63) = 160.
      GRAPHP(64) = 160.
      GRAPHP(65) = 160.
      GRAPHP(66) = 160.
      GRAPHP(67) = 160.
      GRAPHP(68) = 160.
      GRAPHP(69) = 160.
      GRAPHP(70) = 160.
      GRAPHP(71) = 2.
      GRAPHP(72) = DEFOUT(1)
      GRAPHP(73) = 0.
      GRAPHP(74) = DEFOUT(1)
      GRAPHP(75) = DEFOUT(1)
      GRAPHP(76) = DEFOUT(1)
      GRAPHP(77) = DEFOUT(1)
      GRAPHP(91) = -999999.
      GRAPHP(90) = -999999.
      GRAPHP(92) = 0.
      RETURN

      END
