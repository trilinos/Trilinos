C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTRED
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
      COMMON /PSAVE/TDEVP(5,10),TTEXTP(40,10),TVECTP(5,10),
     *       TGRAPH(100,10),TMAPP(11,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM
      EXTERNAL PLTBLK

      IF (DEVP(1).NE.TDEVP(1,IPOPD)) THEN
         CALL PLTSTD(1,TDEVP(1,IPOPD))
      END IF

      IF (DEVP(2).NE.TDEVP(2,IPOPD)) THEN
         CALL PLTSTD(2,TDEVP(2,IPOPD))
      END IF

      IF (DEVP(3).NE.TDEVP(3,IPOPD)) THEN
         CALL PLTSTD(3,TDEVP(3,IPOPD))
      END IF

      DO 2100 I = 1,5
         DEVP(I) = TDEVP(I,IPOPD)
 2100 CONTINUE
      IF (IPOPD.NE.1) THEN
         IPOPD = IPOPD - 1
      END IF

      RETURN

      END
