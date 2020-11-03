C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTCOL(INDEX,R,G,B)
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
      DIMENSION COLARR(3)
      CHARACTER*10 ECOLOR

      FLAG = 0.
      IF (R.LT.0. .OR. R.GT.1.) THEN
         CALL CHRIC(INT(R),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTCOL','Red value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
         FLAG = 1.
      END IF

      IF (G.LT.0. .OR. G.GT.1.) THEN
         CALL CHRIC(INT(G),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTCOL','Green value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
         FLAG = 1.
      END IF

      IF (B.LT.0. .OR. B.GT.1.) THEN
         CALL CHRIC(INT(B),ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTCOL','Blue value of '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-1.',2)
         FLAG = 1.
      END IF

      IF (INDEX.LT.0 .OR. INDEX.GT.255) THEN
         CALL CHRIC(INDEX,ECOLOR,L)
         CALL PLTFLU
         CALL SIORPT('PLTCOL','Color index '//ECOLOR(1:L)//
     *               ' is out of range; range is 0-255.',2)
         FLAG = 1.
      END IF

      IF (FLAG.EQ.1.) THEN
         RETURN

      END IF

      COLARR(1) = R
      COLARR(2) = G
      COLARR(3) = B
      CALL VDSTCO(1,INDEX,COLARR,0)
      IF (INDEX.LE.15) THEN
         PALETT(1,INDEX) = R
         PALETT(2,INDEX) = G
         PALETT(3,INDEX) = B
      END IF

      RETURN

      END
