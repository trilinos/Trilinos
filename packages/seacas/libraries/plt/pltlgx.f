C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLGX(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
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
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      REAL INTERY,MINEXX,MAXEXX

      CALL VECRGS(IABS(NUM),X,XMAX,XMIN)
      CALL VECRGS(IABS(NUM),Y,YMAX,YMIN)
      IF (GRAPHP(22).EQ.1. .OR. GRAPHP(22).EQ.2.) THEN
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGX',
     *     'X <= 0 found on logarithmic X axis; ignoring X values <= 0.'
     *                  ,2)
         END IF

         TEMP = LOG10(XMIN)
         MINEXX = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXX.NE.TEMP) THEN
            MINEXX = MINEXX - 1
         END IF

         TEMP = LOG10(XMAX)
         MAXEXX = INT(TEMP)
         IF (TEMP.GT.0. .AND. TEMP.NE.MAXEXX) THEN
            MAXEXX = MAXEXX + 1
         END IF

         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         CALL PLTINO(YMIN,YMAX,FNLOWY,FNUPPY,INTERY,IEXPY,NMINY)
         YSTART = FNLOWY
         YEND = FNUPPY
         TNEXPY = 10.**IEXPY
         GRAPHP(24) = TENMNX
         GRAPHP(25) = TENMXX
         GRAPHP(28) = FNLOWY*TNEXPY
         GRAPHP(29) = FNUPPY*TNEXPY
         GRAPHP(30) = (FNUPPY-FNLOWY)/INTERY
         GRAPHP(31) = NMINY
         GRAPHP(78) = TENMNX
         GRAPHP(80) = TENMXX
         GRAPHP(83) = FNLOWY*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = FNUPPY*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.3.) THEN
         XMIN = GRAPHP(24)
         XMAX = GRAPHP(25)
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XTEMP,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGX',
     *'User scaling specified minimum X <= 0 on log X axis; using data t
     *o get min X',2)
         END IF

         IF (XMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLGX',
     *'User scaling specified maximum X <= 0 on log X axis; using data t
     *o get max X',2)
         END IF

         MINEXX = LOG10(XMIN)
         MAXEXX = LOG10(XMAX)
         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         TINT = (GRAPHP(29)-GRAPHP(28))/GRAPHP(30)
         IEXPY = NINT(LOG10(ABS(TINT)))
         TNEXPY = 10.**IEXPY
         FNLOWY = GRAPHP(28)/TNEXPY
         FNUPPY = GRAPHP(29)/TNEXPY
         INTERY = (FNUPPY-FNLOWY)/INT(GRAPHP(30))
         NMINY = INT(GRAPHP(31))
         YSTART = FNLOWY
         YEND = FNUPPY
         GRAPHP(78) = XMIN
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = XMAX
         GRAPHP(83) = YSTART*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = YEND*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.4.) THEN
         XMIN = GRAPHP(78)
         XMAX = GRAPHP(80)
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XTEMP,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGX',
     *'Exact scaling specified maximum X <= 0 on log X axis; using data
     *to get min X',2)
         END IF

         IF (XMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLGX',
     *'Exact scaling specified maximum X <= 0 on log X axis; using data
     *to get max X',2)
         END IF

         MINEXX = LOG10(XMIN)
         MAXEXX = LOG10(XMAX)
         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         IEXPY = NINT(LOG10(ABS(GRAPHP(86))))
         TNEXPY = 10.**IEXPY
         YSTART = GRAPHP(83)/TNEXPY
         YEND = GRAPHP(85)/TNEXPY
         FNLOWY = GRAPHP(84)/TNEXPY
         INTERY = GRAPHP(86)/TNEXPY
         NMINY = INT(GRAPHP(87))
         GRAPHP(24) = XMIN
         GRAPHP(25) = XMAX
         GRAPHP(28) = YSTART*TNEXPY
         GRAPHP(29) = YEND*TNEXPY
         GRAPHP(30) = (YSTART-YEND)/INTERY
         GRAPHP(31) = NMINY
      END IF

      IF (GRAPHP(90).NE.-999999.) THEN
         FAC = 10.** (IEXPY-GRAPHP(90))
         IEXPY = INT(GRAPHP(90))
         TNEXPY = 10.**IEXPY
         YSTART = YSTART*FAC
         YEND = YEND*FAC
         FNLOWY = FNLOWY*FAC
         INTERY = INTERY*FAC
      END IF

      IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.4.) THEN
         YSTART = YSTART*TNEXPY
         YEND = YEND*TNEXPY
         FNLOWY = FNLOWY*TNEXPY
         FNUPPY = FNUPPY*TNEXPY
         INTERY = INTERY*TNEXPY
         IEXPY = 0
         TNEXPY = 1.
      END IF

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'x',MINEXX,
     *            MAXEXX,XLAB,XUNIT)

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'y',YSTART,
     *            YEND,FNLOWY,INT(GRAPHP(42)),INTERY,NMINY,YLAB,YUNIT,
     *            IEXPY)

      CALL PLTGM2(MINEXX,MAXEXX,YSTART*TNEXPY,YEND*TNEXPY,GRAPHP(1),
     *            GRAPHP(1)+GRAPHP(3),GRAPHP(2),GRAPHP(2)+GRAPHP(4),
     *            GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)

      RETURN

      END
