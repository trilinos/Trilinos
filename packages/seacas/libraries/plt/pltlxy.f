C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLXY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
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
      REAL MINEXX,MAXEXX,MINEXY,MAXEXY
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      DIMENSION X(1),Y(1)

      XLENT = GRAPHP(3)
      YLENT = GRAPHP(4)
      CALL VECRGS(IABS(NUM),X,XMAX,XMIN)
      CALL VECRGS(IABS(NUM),Y,YMAX,YMIN)
      IF (GRAPHP(22).EQ.1. .OR. GRAPHP(22).EQ.2.) THEN
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
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

         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *     'Y <= 0 found on logarithmic Y axis; ignoring Y values <= 0.'
     *                  ,2)
         END IF

         TEMP = LOG10(YMIN)
         MINEXY = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXY.NE.TEMP) THEN
            MINEXY = MINEXY - 1
         END IF

         TEMP = LOG10(YMAX)
         MAXEXY = INT(TEMP)
         IF (TEMP.GT.0. .AND. MAXEXY.NE.TEMP) THEN
            MAXEXY = MAXEXY + 1
         END IF

         IF (GRAPHP(22).EQ.2.) THEN
            DELX = MAXEXX - MINEXX
            DELY = MAXEXY - MINEXY
            IF (DELX.GT.DELY) THEN
               MAXEXY = MINEXY + MAXEXX - MINEXX
            END IF

            IF (DELX.LT.DELY) THEN
               MAXEXX = MINEXX + MAXEXY - MINEXY
            END IF

            IF (GRAPHP(3).NE.GRAPHP(4)) THEN
               XLENT = AMIN1(GRAPHP(3),GRAPHP(4))
               YLENT = XLENT
            END IF

         END IF

         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
         GRAPHP(24) = TENMNX
         GRAPHP(25) = TENMXX
         GRAPHP(28) = TENMNY
         GRAPHP(29) = TENMXY
         GRAPHP(78) = TENMNX
         GRAPHP(80) = TENMXX
         GRAPHP(83) = TENMNY
         GRAPHP(85) = TENMXY

      ELSE IF (GRAPHP(22).EQ.3.) THEN
         XMIN = GRAPHP(24)
         XMAX = GRAPHP(25)
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XTEMP,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'User scaling specified minimum X <= 0 on log X axis; using data t
     *o get min X',2)
         END IF

         IF (XMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'User scaling specified maximum X <= 0 on log X axis; using data t
     *o get max X',2)
         END IF

         MINEXX = LOG10(XMIN)
         MAXEXX = LOG10(XMAX)
         YMIN = GRAPHP(28)
         YMAX = GRAPHP(29)
         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YTEMP,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'User scaling specified minimum Y <= 0 on log Y axis; using data t
     *o get min Y',2)
         END IF

         IF (YMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'User scaling specified maximum Y <= 0 on log Y axis; using data t
     *o get max Y',2)
         END IF

         MINEXY = LOG10(YMIN)
         MAXEXY = LOG10(YMAX)
         GRAPHP(78) = XMIN
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = XMAX
         GRAPHP(83) = YMIN
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = YMAX

      ELSE IF (GRAPHP(22).EQ.4.) THEN
         XMIN = GRAPHP(78)
         XMAX = GRAPHP(80)
         IF (XMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XTEMP,XMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'Exact scaling specified maximum X <= 0 on log X axis; using data
     *to get min X',2)
         END IF

         IF (XMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),X,XMAX,XTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'Exact scaling specified maximum X <= 0 on log X axis; using data
     *to get max X',2)
            GRAPHP(24) = XMIN
            GRAPHP(25) = XMAX
            GRAPHP(28) = YMIN
            GRAPHP(29) = YMAX
         END IF

         MINEXX = LOG10(XMIN)
         MAXEXX = LOG10(XMAX)
         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         YMIN = GRAPHP(83)
         YMAX = GRAPHP(85)
         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YTEMP,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'Exact scaling specified maximum Y <= 0 on log Y axis; using data
     *to get min Y',2)
         END IF

         IF (YMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLXY',
     *'Exact scaling specified maximum Y <= 0 on log Y axis; using data
     *to get max Y',2)
         END IF

         MINEXY = LOG10(YMIN)
         MAXEXY = LOG10(YMAX)
         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
      END IF

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'x',MINEXX,MAXEXX,
     *            XLAB,XUNIT)

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'y',MINEXY,MAXEXY,
     *            YLAB,YUNIT)

      CALL PLTGM2(MINEXX,MAXEXX,MINEXY,MAXEXY,GRAPHP(1),GRAPHP(1)+XLENT,
     *            GRAPHP(2),GRAPHP(2)+YLENT,GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)

      RETURN

      END
