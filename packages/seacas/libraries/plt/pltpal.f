C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTPAL(COL,R,G,B)
      CHARACTER*10 ECOLOR
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

      IF (COL.LT.8. .OR. COL.GT.15.) THEN
         CALL CHRIC(INT(COL),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Illegal palette color '//ECOLOR(1:L)//
     *               ' passed to PLTPAL; range is 8-15.',2)
         RETURN

      END IF

      IF (R.LT.0. .OR. R.GT.1.) THEN
         CALL CHRIC(INT(R),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Red value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      IF (G.LT.0. .OR. G.GT.1.) THEN
         CALL CHRIC(INT(G),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Green value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      IF (B.LT.0. .OR. B.GT.1.) THEN
         CALL CHRIC(INT(B),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTPAL','Blue value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
      END IF

      PALETT(1,INT(COL)) = R
      PALETT(2,INT(COL)) = G
      PALETT(3,INT(COL)) = B
      CALL VDSTCO(1,INT(COL),PALETT(1,INT(COL)),0)
      RETURN

      END
