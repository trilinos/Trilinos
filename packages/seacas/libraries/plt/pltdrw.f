C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTDRW(X,Y)
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
      DIMENSION ISTYLE(0:7,5),VECTOR(7)
      REAL SAVLEN,LINLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV
      EXTERNAL PLTBLK

      SAVE VECTOR
      DATA ISTYLE/1,0,1,0,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,1,1,0,1,1,1,
     *     1,1,1,0,0,1,1,1,1,0,0,0,0/

      IF (VECTP(1).EQ.0.) THEN
         XCUR = X
         YCUR = Y
         RETURN

      END IF

      CALL VDIQOS(VECTOR)
      CALL VDSTLS(0)
      DX = X - XCUR
      DY = Y - YCUR
      LINLEN = SQRT(DX*DX+DY*DY)
      NSTYLE = INT(VECTOR(4))
      IF (NSTYLE.EQ.0 .OR. LINLEN.EQ.0) THEN
         CALL PLTLIG(X,Y)

      ELSE
         IDSHNO = IDSHSV
         DSHDRN = SAVLEN
         DSHLEN = MAX(0.02*VECTOR(5),.002)
         DSHLEN = MIN(DSHLEN,.005)
         DSHDX = 0.5*DX* (DSHLEN/LINLEN)
         DSHDY = 0.5*DY* (DSHLEN/LINLEN)
         DRWLEN = LINLEN - DSHLEN + DSHDRN
         IF (DRWLEN.LE.0) THEN
            XX = XCUR
            YY = YCUR

         ELSE
            IDSHNO = IDSHNO + 1
            XX = XCUR + DX* (DSHLEN-DSHDRN)/LINLEN
            YY = YCUR + DY* (DSHLEN-DSHDRN)/LINLEN
            IF (ISTYLE(MOD(IDSHNO,8),NSTYLE).NE.1) THEN
               IF (ISTYLE(MOD(IDSHNO+1,8),NSTYLE).EQ.1) THEN
                  CALL PLTMOV(XX,YY)
               END IF

            ELSE IF (ISTYLE(MOD(IDSHNO+1,8),NSTYLE).NE.1 .AND.
     *               DSHDRN.LT.DSHLEN/2) THEN
               CALL PLTLIG(0.5* (XCUR+XX),0.5* (YCUR+YY))
            END IF

            DSHDRN = 0.
            NUMDSH = MAX(NINT(DRWLEN/DSHLEN),1)
            DO 2020 I = 1,NUMDSH - 1
               IDSHNO = IDSHNO + 1
               X2 = XX + DSHDX
               Y2 = YY + DSHDY
               XX = X2 + DSHDX
               YY = Y2 + DSHDY
               IF (ISTYLE(MOD(IDSHNO,8),NSTYLE).NE.1) THEN
                  IF (ISTYLE(MOD(IDSHNO+1,8),NSTYLE).EQ.1) THEN
                     CALL PLTMOV(XX,YY)
                  END IF

               ELSE IF (ISTYLE(MOD(IDSHNO+1,8),NSTYLE).NE.1) THEN
                  CALL PLTLIG(X2,Y2)
               END IF

 2020       CONTINUE
         END IF

         DX = X - XX
         DY = Y - YY
         PRVDRN = DSHDRN
         DSHDRN = SQRT(DX*DX+DY*DY) + PRVDRN
         IF (ISTYLE(MOD(IDSHNO+1,8),NSTYLE).NE.1) THEN
            CALL PLTMOV(X,Y)

         ELSE IF (ISTYLE(MOD(IDSHNO+2,8),NSTYLE).NE.1 .AND.
     *            DSHDRN.GT.DSHLEN/2) THEN
            IF (PRVDRN.LT.DSHLEN/2) THEN
               CALL PLTLIG(XX+DSHDX,YY+DSHDY)
            END IF

            CALL PLTMOV(X,Y)

         ELSE
            CALL PLTLIG(X,Y)
         END IF

      END IF

      IDSHSV = IDSHNO
      SAVLEN = DSHDRN
      XCUR = X
      YCUR = Y
      CALL VDSTLS(INT(VECTOR(4)))
      RETURN

      END
