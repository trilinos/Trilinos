C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTREV
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

      IF (VECTP(1).NE.TVECTP(1,IPOPV)) THEN
         CALL PLTSTV(1,TVECTP(1,IPOPV))
      END IF

      IF (VECTP(2).NE.TVECTP(2,IPOPV)) THEN
         CALL PLTSTV(2,TVECTP(2,IPOPV))
      END IF

      DO 2200 I = 1,5
         VECTP(I) = TVECTP(I,IPOPV)
 2200 CONTINUE
      IF (IPOPV.NE.1) THEN
         IPOPV = IPOPV - 1
      END IF

      RETURN

      END
