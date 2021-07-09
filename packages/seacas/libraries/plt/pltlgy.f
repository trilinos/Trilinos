C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTLGY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
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
      REAL INTERX,MINEXY,MAXEXY
      DIMENSION X(1),Y(1)

      CALL VECRGS(IABS(NUM),X,XMAX,XMIN)
      CALL VECRGS(IABS(NUM),Y,YMAX,YMIN)
      IF (GRAPHP(22).EQ.1. .OR. GRAPHP(22).EQ.2.) THEN
         CALL PLTINO(XMIN,XMAX,FNLOWX,FNUPPX,INTERX,IEXPX,NMINX)
         TNEXPX = 10.**IEXPX
         XSTART = FNLOWX
         XEND = FNUPPX
         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGY',
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

         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
         GRAPHP(24) = FNLOWX*TNEXPX
         GRAPHP(25) = FNUPPX*TNEXPX
         GRAPHP(26) = (FNUPPX-FNLOWX)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = TENMNY
         GRAPHP(29) = TENMXY
         GRAPHP(78) = FNLOWX*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = FNUPPX*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = TENMNY
         GRAPHP(85) = TENMXY

      ELSE IF (GRAPHP(22).EQ.3.) THEN
         TINT = (GRAPHP(25)-GRAPHP(24))/GRAPHP(26)
         IEXPX = NINT(LOG10(ABS(TINT)))
         TNEXPX = 10.**IEXPX
         FNLOWX = GRAPHP(24)/TNEXPX
         FNUPPX = GRAPHP(25)/TNEXPX
         INTERX = (FNUPPX-FNLOWX)/INT(GRAPHP(26))
         NMINX = INT(GRAPHP(27))
         XSTART = FNLOWX
         XEND = FNUPPX
         YMIN = GRAPHP(28)
         YMAX = GRAPHP(29)
         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YTEMP,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGY',
     *'User scaling specified minimum Y <= 0 on log Y axis; using data t
     *o get min Y',2)
         END IF

         IF (YMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLGY',
     *'User scaling specified maximum Y <= 0 on log Y axis; using data t
     *o get max Y',2)
         END IF

         MINEXY = LOG10(YMIN)
         MAXEXY = LOG10(YMAX)
         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
         GRAPHP(78) = XSTART*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = XEND*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = YMIN
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = YMAX

      ELSE IF (GRAPHP(22).EQ.4.) THEN
         IEXPX = NINT(LOG10(ABS(GRAPHP(81))))
         TNEXPX = 10.**IEXPX
         XSTART = GRAPHP(78)/TNEXPX
         XEND = GRAPHP(80)/TNEXPX
         FNLOWX = GRAPHP(79)/TNEXPX
         INTERX = GRAPHP(81)/TNEXPX
         NMINX = INT(GRAPHP(82))
         YMIN = GRAPHP(83)
         YMAX = GRAPHP(85)
         IF (YMIN.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YTEMP,YMIN)
            CALL PLTFLU
            CALL SIORPT('PLTLGY',
     *'Exact scaling specified maximum Y <= 0 on log Y axis; using data
     *to get min Y',2)
         END IF

         IF (YMAX.LE.0) THEN
            CALL VECRGP(IABS(NUM),Y,YMAX,YTEMP)
            CALL PLTFLU
            CALL SIORPT('PLTLGY',
     *'Exact scaling specified maximum Y <= 0 on log Y axis; using data
     *to get max Y',2)
         END IF

         MINEXY = LOG10(YMIN)
         MAXEXY = LOG10(YMAX)
         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
         GRAPHP(24) = XSTART*TNEXPX
         GRAPHP(25) = XEND*TNEXPX
         GRAPHP(26) = (XSTART-XEND)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = YMIN
         GRAPHP(29) = YMAX
      END IF

      IF (GRAPHP(91).NE.-999999.) THEN
         FAC = 10.** (IEXPX-GRAPHP(91))
         IEXPX = INT(GRAPHP(91))
         TNEXPX = 10.**IEXPX
         XSTART = XSTART*FAC
         XEND = XEND*FAC
         FNLOWX = FNLOWX*FAC
         INTERX = INTERX*FAC
      END IF

      IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.4.) THEN
         XSTART = XSTART*TNEXPX
         XEND = XEND*TNEXPX
         FNLOWX = FNLOWX*TNEXPX
         FNUPPX = FNUPPX*TNEXPX
         INTERX = INTERX*TNEXPX
         IEXPX = 0
         TNEXPX = 1.
      END IF

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'x',XSTART,
     *            XEND,FNLOWX,INT(GRAPHP(41)),INTERX,NMINX,XLAB,XUNIT,
     *            IEXPX)

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'y',MINEXY,
     *            MAXEXY,YLAB,YUNIT)

      CALL PLTGM2(XSTART*TNEXPX,XEND*TNEXPX,MINEXY,MAXEXY,GRAPHP(1),
     *            GRAPHP(1)+GRAPHP(3),GRAPHP(2),GRAPHP(2)+GRAPHP(4),
     *            GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)

      RETURN

      END
