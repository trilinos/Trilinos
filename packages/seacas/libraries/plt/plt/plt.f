C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

      BLOCK DATA PLTBLK
      REAL SAVLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV
      COMMON /PSAVE/TDEVP(5,10),TTEXTP(40,10),TVECTP(5,10),
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM
      COMMON /MPSTCK/SVMAP(195,10),MAPDEP
      DATA IPOPD/0/,IPOPT/0/,IPOPV/0/,IPOPG/0/,IPOPM/0/
      DATA SAVLEN/0./,IDSHSV/0/
      DATA MAPDEP/0/
      end

      SUBROUTINE PLTARR(XTAIL,YTAIL,XHEAD,YHEAD,THETA,ARRLEN)
      REAL XTAIL,YTAIL
      REAL XHEAD,YHEAD
      REAL THETA
      REAL ARRLEN
      DATA SMALL/1.E-3/

      IF (ABS(XTAIL-XHEAD).LE.SMALL .AND.
     *    ABS(YTAIL-YHEAD).LE.SMALL) THEN
         RETURN

      END IF

      CALL PLTVCT(1,XTAIL,YTAIL,XHEAD,YHEAD)
      IF (ARRLEN.LE.0.) THEN
         RETURN

      END IF

      DX = XHEAD - XTAIL
      DY = YHEAD - YTAIL
      ALPHA = ATAN2(DY,DX)
      A = -ARRLEN*COS(ALPHA)*COS(THETA) + ARRLEN*SIN(ALPHA)*SIN(THETA) +
     *    XHEAD
      B = -ARRLEN*COS(ALPHA)*SIN(THETA) - ARRLEN*SIN(ALPHA)*COS(THETA) +
     *    YHEAD
      CALL PLTVCT(1,XHEAD,YHEAD,A,B)
      A = -ARRLEN*COS(ALPHA)*COS(-THETA) +
     *    ARRLEN*SIN(ALPHA)*SIN(-THETA) + XHEAD
      B = -ARRLEN*COS(ALPHA)*SIN(-THETA) -
     *    ARRLEN*SIN(ALPHA)*COS(-THETA) + YHEAD
      CALL PLTVCT(1,XHEAD,YHEAD,A,B)
      RETURN

      END
      SUBROUTINE PLTBGN
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

      CALL VDNWPG
      CALL PLTMOV(0.,0.)
      DO 2000 I = 4,11
         TEXTP(I) = 0.
 2000 CONTINUE
      RETURN

      END
      SUBROUTINE PLTBEL

      CALL VDBELL
      RETURN

      END
      LOGICAL FUNCTION PLTCRS(X,Y,KEY)
      CHARACTER KEY*1
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

      PLTCRS = .FALSE.
      CALL VDSTLA(X,Y)
      CALL VDAKGL(ICHAR,XT,YT)
      KEY = CHAR(ICHAR)
      CALL PLTD2P(XT,YT,X,Y)
      IF (X.LT.0. .OR. X.GT.DEVP(4) .OR. Y.LT.0. .OR. Y.GT.DEVP(5)) THEN
         CALL PLTFLU
         CALL SIORPT('PLTCRS',
     *              'The cursor is out of range of the current viewport'
     *               ,2)
         RETURN

      END IF

      PLTCRS = .TRUE.
      RETURN

      END
      SUBROUTINE PLTD2P(XD,YD,XP,YP)
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
      DATA DTR/.01745329/

      THETA = MAPP(5)*DTR
      SINT = SIN(THETA)
      COST = COS(THETA)
      A = (XD-MAPP(6))* (DEVP(4)/ (MAPP(8)-MAPP(6))) - DEVP(4)/2. -
     *    MAPP(3)
      B = (YD-MAPP(7))* (DEVP(5)/ (MAPP(9)-MAPP(7))) - DEVP(5)/2. -
     *    MAPP(4)
      XP = (A*COST+B*SINT)/MAPP(1) + DEVP(4)/2.
      YP = (B*COST-A*SINT)/MAPP(2) + DEVP(5)/2.
      RETURN

      END
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
      NSTYLE = VECTOR(4)
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
      SUBROUTINE PLTEND

      CALL VDFRAM(1)
      CALL VDTERM
      RETURN

      END
      SUBROUTINE PLTFLU

      CALL VDBUFL
      RETURN

      END
      SUBROUTINE PLTFRM(TYPE)
      INTEGER TYPE

      CALL VDFRAM(TYPE)
      RETURN

      END
      SUBROUTINE PLTCPY
      INTEGER SUPPO

      CALL PLTFLU
      CALL VDIQES(100,SUPPO)
      IF (SUPPO.NE.0) THEN
         CALL VDESCP(100,0)
      END IF

      RETURN

      END
      SUBROUTINE PLTINT
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
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      IF (FIRST) THEN
         ASP = 10./7.5
         CALL VDINIT(ASP,5)
         CALL VDFRAM(0)
         CALL PLTIQD(DEVCAP)
         CALL VDSTCS(1./58.)
         CALL VDIQOS(DEFOUT)
         DEFOUT(1) = 7.
         FIRST = .FALSE.
      END IF

      CALL VDIQND(DEVP(4),DEVP(5))
      CALL PLTRST
      CALL PLTRSD
      CALL PLTRSV
      CALL PLTRSG
      CALL PLTRSM
      CALL PLTRSC
      CALL MPINIT
      CALL MEMINI(MEMORY,1000)
      CALL PLTFLU
      RETURN

      END
      SUBROUTINE PLTIQD(ARRAY)
      DIMENSION ARRAY(*)

      DO 2040 I = 1,23
         CALL VDIQDC(I,ARRAY(I))
 2040 CONTINUE
      RETURN

      END
      SUBROUTINE PLTLIG(X,Y)
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

      IF (MAPP(10).EQ.1.) THEN
         CALL PLTP2D(X,Y,XP,YP)
         CALL PLTP2D(XCUR,YCUR,XC,YC)

      ELSE
         XP = X
         YP = Y
         XC = XCUR
         YC = YCUR
      END IF

      MASK = -1
      CALL PLTVWV(MAPP(6),MAPP(8),1,MASK,XC,YC,XP,YP)
      IF (MASK.EQ.-1) THEN
         CALL VDMOVA(XC,YC)
         CALL VDLINA(XP,YP)
      END IF

      RETURN

      END
      SUBROUTINE PLTP2D(X,Y,XN,YN)
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
      DATA DTR/.01745329/

      XT = X - DEVP(4)/2.
      YT = Y - DEVP(5)/2.
      THETA = MAPP(5)*DTR
      XT = XT*MAPP(1)
      YT = YT*MAPP(2)
      XN = (XT*COS(THETA)-YT*SIN(THETA)) + DEVP(4)/2. + MAPP(3)
      YN = (YT*COS(THETA)+XT*SIN(THETA)) + DEVP(5)/2. + MAPP(4)
      SX = (MAPP(8)-MAPP(6))/DEVP(4)
      SY = (MAPP(9)-MAPP(7))/DEVP(5)
      XN = XN*SX + MAPP(6)
      YN = YN*SY + MAPP(7)
      RETURN

      END
      SUBROUTINE PLTMOV(X,Y)
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
      REAL SAVLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV

      IF (MAPP(10).EQ.1.) THEN
         CALL PLTP2D(X,Y,XP,YP)

      ELSE
         XP = X
         YP = Y
      END IF

      MASK = -1
      CALL PLTVWP(MAPP(6),MAPP(8),1,MASK,XP,YP)
      IF (MASK.EQ.-1) THEN
         CALL VDMOVA(XP,YP)
      END IF

      XCUR = X
      YCUR = Y
      IDSHSV = 0
      SAVLEN = 0.
      RETURN

      END
      SUBROUTINE PLTPNT(N,X,Y)
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
      REAL SAVLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV
      DIMENSION X(*),Y(*)

      DO 2060 I = 1,N
         IF (MAPP(10).EQ.1.) THEN
            CALL PLTP2D(X(I),Y(I),XP,YP)

         ELSE
            XP = X(I)
            YP = Y(I)
         END IF

         MASK = -1
         CALL PLTVWP(MAPP(6),MAPP(8),1,MASK,XP,YP)
         IF (MASK.EQ.-1) THEN
            CALL VDPNTA(XP,YP)
         END IF

 2060 CONTINUE
      XCUR = X(N)
      YCUR = Y(N)
      IDSHSV = 0
      SAVLEN = 0.
      RETURN

      END
      SUBROUTINE PLTPLY(N,XA,YA)
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
      DIMENSION XA(*),YA(*)
      LOGICAL STATUS,MEMFRE
      REAL MEMORY

      IF (N.LE.0) THEN
         RETURN

      END IF

      IVX = MEMALL(N,MEMORY)
      IVY = MEMALL(N,MEMORY)
      IF (MAPP(10).EQ.1.) THEN
         DO 2080 I = 1,N
            IMX = IVX + I - 1
            IMY = IVY + I - 1
            CALL PLTP2D(XA(I),YA(I),MEMORY(IMX),MEMORY(IMY))
 2080    CONTINUE

      ELSE
         DO 2100 I = 1,N
            IMX = IVX + I - 1
            IMY = IVY + I - 1
            MEMORY(IMX) = XA(I)
            MEMORY(IMY) = YA(I)
 2100    CONTINUE
      END IF

      ICX = MEMALL(N+10,MEMORY)
      ICY = MEMALL(N+10,MEMORY)
      IZ1 = MEMALL(N,MEMORY)
      IZ2 = MEMALL(N+10,MEMORY)
      NO = N + 10
      CALL PLTVWG(MAPP(6),MAPP(8),N,MEMORY(IVX),MEMORY(IVY),MEMORY(IZ1),
     *            NO,MEMORY(ICX),MEMORY(ICY),MEMORY(IZ2))
      IF (NO.GT.0) THEN
         CALL VDPOLY(MEMORY(ICX),MEMORY(ICY),NO)
      END IF

      STATUS = MEMFRE(ICX,MEMORY)
      STATUS = MEMFRE(ICY,MEMORY)
      STATUS = MEMFRE(IZ1,MEMORY)
      STATUS = MEMFRE(IZ2,MEMORY)
      STATUS = MEMFRE(IVX,MEMORY)
      STATUS = MEMFRE(IVY,MEMORY)
      RETURN

      END
      SUBROUTINE PLTPTM(N,MASK,X,Y)
      DIMENSION X(*),Y(*),MASK(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      J = 0
      J1 = J
      KM = 0
 2120 IF (.NOT. (J.LT.N)) GO TO 2130
      JN = MIN(N-J,32)
      KM = KM + 1
      M = MASK(KM)
      J = J + JN
      IF (M.EQ.0) THEN
         GO TO 2120

      END IF

      DO 2140 K = 1,JN
         IF (IAND(M,IZBIT(K)).NE.0) THEN
            CALL PLTPNT(1,X(J1+K),Y(J1+K))
         END IF

 2140 CONTINUE
      GO TO 2120

 2130 CONTINUE
      RETURN

      END
      SUBROUTINE PLTRDC(XNDC,YNDC)

      CALL VDIQND(XNDC,YNDC)
      RETURN

      END
      SUBROUTINE PLTRXY(X,Y)
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

      X = XCUR
      Y = YCUR
      RETURN

      END
      SUBROUTINE PLTSTA
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

      CALL PLTIQD(DEVCAP)
      CALL VDSTCS(1./58.)
      CALL VDIQOS(DEFOUT)
      CALL VDIQND(DEVP(4),DEVP(5))
      CALL PLTRST
      CALL PLTRSD
      CALL PLTRSV
      CALL PLTRSG
      CALL PLTRSM
      CALL PLTRSC
      CALL MPINIT
      CALL MEMINI(MEMORY,1000)
      CALL PLTFLU
      RETURN

      END
      SUBROUTINE PLTVCT(N,XX0,YY0,XX1,YY1)
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
      DIMENSION XX0(*),YY0(*),XX1(*),YY1(*)

      IF (VECTP(1).LE.0. .OR. VECTP(2).LE.0.) THEN
         XCUR = XX1(N)
         YCUR = YY1(N)
         RETURN

      END IF

      DO 2160 I = 1,N
         IF (XX0(I).NE.XCUR .OR. YY0(I).NE.YCUR) THEN
            CALL PLTMOV(XX0(I),YY0(I))
         END IF

         CALL PLTDRW(XX1(I),YY1(I))
         XCUR = XX1(I)
         YCUR = YY1(I)
 2160 CONTINUE
      RETURN

      END
      SUBROUTINE PLTVCM(N,MASK,XX0,YY0,XX1,YY1)
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
      DIMENSION XX0(*),YY0(*),XX1(*),YY1(*),MASK(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      IF (VECTP(1).LE.0. .OR. VECTP(2).LE.0.) THEN
         XCUR = XX1(N)
         YCUR = YY1(N)
         RETURN

      END IF

      J = 0
      KM = 0
 2180 IF (.NOT. (J.LT.N)) GO TO 2190
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2180

      END IF

      DO 2200 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2200

         END IF

         IF (XX0(K+J1).NE.XCUR .OR. YY0(K+J1).NE.YCUR) THEN
            CALL PLTMOV(XX0(K+J1),YY0(K+J1))
         END IF

         CALL PLTDRW(XX1(K+J1),YY1(K+J1))
         XCUR = XX1(K+J1)
         YCUR = YY1(K+J1)
 2200 CONTINUE
      GO TO 2180

 2190 CONTINUE
      RETURN

      END
      SUBROUTINE PLTWAI

      CALL VDWAIT
      RETURN

      END
      SUBROUTINE PLTCG2(N,XV,YV,NO,XVO,YVO,C1,C2)
      INTEGER N
      REAL XV(*),YV(*)
      INTEGER NO
      REAL XVO(*),YVO(*)
      REAL C1(2),C2(2)
      REAL S(2),P(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTCG2')

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2000 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         CP = (C2(1)-C1(1))* (P(2)-C1(2)) - (C2(2)-C1(2))* (P(1)-C1(1))
         IF (CP.GE.0.) THEN
            CP = (C2(1)-C1(1))* (S(2)-C1(2)) -
     *           (C2(2)-C1(2))* (S(1)-C1(1))
            IF (CP.GE.0.) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               FP = (S(2)-C1(2))* (C2(1)-C1(1)) -
     *              (S(1)-C1(1))* (C2(2)-C1(2))
               FQ = (P(2)-C1(2))* (C2(1)-C1(1)) -
     *              (P(1)-C1(1))* (C2(2)-C1(2))
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            CP = (C2(1)-C1(1))* (S(2)-C1(2)) -
     *           (C2(2)-C1(2))* (S(1)-C1(1))
            IF (CP.GE.0.) THEN
               FP = (S(2)-C1(2))* (C2(1)-C1(1)) -
     *              (S(1)-C1(1))* (C2(2)-C1(2))
               FQ = (P(2)-C1(2))* (C2(1)-C1(1)) -
     *              (P(1)-C1(1))* (C2(2)-C1(2))
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2000 CONTINUE
      RETURN

      END
      SUBROUTINE PLTCP2(N,MASK,PX,PY,C1,C2)
      DIMENSION MASK(*),PX(*),PY(*),C1(*),C2(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      CX = C1(1)
      CY = C1(2)
      DX = C2(1) - CX
      DY = C2(2) - CY
      J = 0
      KM = 0
 2020 IF (.NOT. (J.LT.N)) GO TO 2030
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2020

      END IF

      DO 2040 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).NE.0) THEN
            FP = (PY(J1+K)-CY)*DX - (PX(J1+K)-CX)*DY
            IF (FP.LT.0.) THEN
               M = IAND(M,NOT(JB))
            END IF

         END IF

 2040 CONTINUE
      MASK(KM) = M
      GO TO 2020

 2030 CONTINUE
      RETURN

      END
      SUBROUTINE PLTCP3(N,MASK,PX,PY,PZ,V,Q)
      DIMENSION MASK(*),PX(*),PY(*),PZ(*),V(*),Q(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      J = 0
      KM = 0
 2060 IF (.NOT. (J.LT.N)) GO TO 2070
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2060

      END IF

      DO 2080 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).NE.0) THEN
            FP = (PX(J1+K)-V(1))*Q(1) + (PY(J1+K)-V(2))*Q(2) +
     *           (PZ(J1+K)-V(3))*Q(3)
            IF (FP.LT.0.) THEN
               M = IAND(M,NOT(JB))
            END IF

         END IF

 2080 CONTINUE
      MASK(KM) = M
      GO TO 2060

 2070 CONTINUE
      RETURN

      END
      SUBROUTINE PLTCV2(N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY,C1,C2)
      DIMENSION MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),QQX(*),
     *          QQY(*),C1(*),C2(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      CX = C1(1)
      CY = C1(2)
      DX = C2(1) - CX
      DY = C2(2) - CY
      J = 0
      KM = 0
 2100 IF (.NOT. (J.LT.N)) GO TO 2110
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2100

      END IF

      DO 2120 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(JB,M).EQ.0) THEN
            GO TO 2120

         END IF

         X1 = PX(J1+K)
         Y1 = PY(J1+K)
         X2 = QX(J1+K)
         Y2 = QY(J1+K)
         FP = (Y1-CY)*DX - (X1-CX)*DY
         FQ = (Y2-CY)*DX - (X2-CX)*DY
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))

         ELSE IF (FP.LT.0.) THEN
            XL = FQ/ (FQ-FP)
            X1 = X2 + XL* (X1-X2)
            Y1 = Y2 + XL* (Y1-Y2)

         ELSE IF (FQ.LT.0.) THEN
            XL = FP/ (FP-FQ)
            X2 = X1 + XL* (X2-X1)
            Y2 = Y1 + XL* (Y2-Y1)
         END IF

         PPX(K+J1) = X1
         PPY(K+J1) = Y1
         QQX(K+J1) = X2
         QQY(K+J1) = Y2
 2120 CONTINUE
      MASK(KM) = M
      GO TO 2100

 2110 CONTINUE
      RETURN

      END
      SUBROUTINE PLTCV3(N,MASK,PX,PY,PZ,QX,QY,QZ,PPX,PPY,PPZ,QQX,QQY,
     *                  QQZ,V,Q)
      DIMENSION MASK(*),PX(*),PY(*),PZ(*),QX(*),QY(*),QZ(*),PPX(*),
     *          PPY(*),PPZ(*),QQX(*),QQY(*),QQZ(*),V(*),Q(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      CALL CPUMVU(PX,PPX,N)
      CALL CPUMVU(PY,PPY,N)
      CALL CPUMVU(PZ,PPZ,N)
      CALL CPUMVU(QX,QQX,N)
      CALL CPUMVU(QY,QQY,N)
      CALL CPUMVU(QZ,QQZ,N)
      J = 0
      KM = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J
      KM = 1 + KM
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2140

      END IF

      DO 2160 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(JB,M).EQ.0) THEN
            GO TO 2160

         END IF

         FP = (PPX(J1+K)-V(1))*Q(1) + (PPY(J1+K)-V(2))*Q(2) +
     *        (PPZ(J1+K)-V(3))*Q(3)
         FQ = (QQX(J1+K)-V(1))*Q(1) + (QQY(J1+K)-V(2))*Q(2) +
     *        (QQZ(J1+K)-V(3))*Q(3)
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))

         ELSE IF (FP.LT.0.) THEN
            XL = FP/ (FP-FQ)
            PPX(J1+K) = PPX(J1+K) + XL* (QQX(J1+K)-PPX(J1+K))
            PPY(J1+K) = PPY(J1+K) + XL* (QQY(J1+K)-PPY(J1+K))
            PPZ(J1+K) = PPZ(J1+K) + XL* (QQZ(J1+K)-PPZ(J1+K))

         ELSE IF (FQ.LT.0.) THEN
            XL = FQ/ (FQ-FP)
            QQX(J1+K) = QQX(J1+K) + XL* (PPX(J1+K)-QQX(J1+K))
            QQY(J1+K) = QQY(J1+K) + XL* (PPY(J1+K)-QQY(J1+K))
            QQZ(J1+K) = QQZ(J1+K) + XL* (PPZ(J1+K)-QQZ(J1+K))
         END IF

 2160 CONTINUE
      MASK(KM) = M
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
      SUBROUTINE PLTVWG(PLL,PUR,N,XV,YV,ZV,NO,XVO,YVO,ZVO)
      DIMENSION XV(*),YV(*),ZV(*),XVO(*),YVO(*),ZVO(*),PLL(*),PUR(*)
      DIMENSION XVT(50),YVT(50),QUR(2),QLL(2)

      QUR(1) = PUR(1) + .0001
      QUR(2) = PUR(2) + .0001
      QLL(1) = PLL(1) - .0001
      QLL(2) = PLL(2) - .0001
      NOSAVE = NO
      NT = 50
      CALL PLTLI1(QLL,QUR,N,XV,YV,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI2(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      NT = 50
      CALL PLTLI3(QLL,QUR,NO,XVO,YVO,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI4(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      DO 2180 I = 1,NO
         ZVO(I) = PLTPGZ(NO,XV,YV,ZV,XVO(I),YVO(I))
 2180 CONTINUE
      RETURN

      END
      FUNCTION PLTPGZ(N,X,Y,Z,XQ,YQ)
      DIMENSION X(3),Y(3),Z(3),D(2)

      X21 = X(2) - X(1)
      X31 = X(3) - X(1)
      Y21 = Y(2) - Y(1)
      Y31 = Y(3) - Y(1)
      Z21 = Z(2) - Z(1)
      Z31 = Z(3) - Z(1)
      A = Y31*Z21 - Z31*Y21
      B = Z31*X21 - X31*Z21
      C = X31*Y21 - Y31*X21
      IF (C.EQ.0.) THEN
         DO 2200 I = 1,N - 2
            D(1) = X(I+1) - X(I)
            D(2) = Y(I+1) - Y(I)
            DDD = D(I)*D(I) + D(I+1)*D(I+1)
            IF (DDD.NE.0.) THEN
               GO TO 2210

            END IF

 2200    CONTINUE
 2210    CONTINUE
         IF (DDD.EQ.0.) THEN
            PLTPGZ = 0.0
            RETURN

         END IF

         DDP = D(1)*XQ + D(2)*YQ
         DDP1 = D(1)*X(1) + D(2)*Y(1)
         ALPHA = (DDP-DDP1)/DDD
         PLTPGZ = Z(1) + ALPHA* (Z(2)-Z(1))
         RETURN

      END IF

      PLTPGZ = ((A* (X(1)-XQ)+B* (Y(1)-YQ))/C) + Z(1)
      RETURN

      END
      SUBROUTINE PLTLI1(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI1')
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2220 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(2) .GE. PLL(2)
         IF (INSIDE) THEN
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            INSIDE = S(2) .GE. PLL(2)
            IF (INSIDE) THEN
               TEMP = PUR(1) - PLL(1)
               FP = (S(2)-PLL(2))*TEMP
               FQ = (P(2)-PLL(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2220 CONTINUE
      RETURN

      END
      SUBROUTINE PLTLI2(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI2')
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2240 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(1) .LE. PUR(1)
         IF (INSIDE) THEN
            INSIDE = S(1) .LE. PUR(1)
            IF (INSIDE) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               TEMP = PUR(2) - PLL(2)
               FP = - (S(1)-PUR(1))*TEMP
               FQ = - (P(1)-PUR(1))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            INSIDE = S(1) .LE. PUR(1)
            IF (INSIDE) THEN
               TEMP = PUR(2) - PLL(2)
               FP = - (S(1)-PUR(1))*TEMP
               FQ = - (P(1)-PUR(1))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2240 CONTINUE
      RETURN

      END
      SUBROUTINE PLTLI3(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI3')
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2260 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(2) .LE. PUR(2)
         IF (INSIDE) THEN
            INSIDE = S(2) .LE. PUR(2)
            IF (INSIDE) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               TEMP = PLL(1) - PUR(1)
               FP = (S(2)-PUR(2))*TEMP
               FQ = (P(2)-PUR(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            INSIDE = S(2) .LE. PUR(2)
            IF (INSIDE) THEN
               TEMP = PLL(1) - PUR(1)
               FP = (S(2)-PUR(2))*TEMP
               FQ = (P(2)-PUR(2))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2260 CONTINUE
      RETURN

      END
      SUBROUTINE PLTLI4(PLL,PUR,N,XV,YV,NO,XVO,YVO)
      DIMENSION PLL(2),PUR(2),XV(*),YV(*),XVO(*),YVO(*)
      DIMENSION P(2),S(2)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='PLTLI4')
      LOGICAL INSIDE

      NMAX = NO
      NO = 0
      S(1) = XV(N)
      S(2) = YV(N)
      DO 2280 I = 1,N
         P(1) = XV(I)
         P(2) = YV(I)
         INSIDE = P(1) .GE. PLL(1)
         IF (INSIDE) THEN
            INSIDE = S(1) .GE. PLL(1)
            IF (INSIDE) THEN
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)

            ELSE
               TEMP = PLL(2) - PUR(1)
               FP = - (S(1)-PLL(1))*TEMP
               FQ = - (P(1)-PLL(1))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1)
               YVO(NO) = P(2)
            END IF

         ELSE
            INSIDE = S(1) .GE. PLL(1)
            IF (INSIDE) THEN
               TEMP = PLL(2) - PUR(1)
               FP = - (S(1)-PLL(1))*TEMP
               FQ = - (P(1)-PLL(1))*TEMP
               XL = FQ/ (FQ-FP)
               NO = NO + 1
               IF (NO.GT.NMAX) THEN
                  CALL PLTFLU
                  CALL SIORPT(SUBNAM,
     *   'Not enough space for output vertices; polygon clip unfinished'
     *                        ,2)
                  RETURN

               END IF

               XVO(NO) = P(1) + XL* (S(1)-P(1))
               YVO(NO) = P(2) + XL* (S(2)-P(2))
            END IF

         END IF

         S(1) = P(1)
         S(2) = P(2)
 2280 CONTINUE
      RETURN

      END
      SUBROUTINE PLTVWP(PLL,PUR,N,MASK,PX,PY)
      REAL PLL(2),PUR(2)
      INTEGER MASK(*)
      REAL PX(*),PY(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      PUR1 = PUR(1) + .0001
      PUR2 = PUR(2) + .0001
      PLL1 = PLL(1) - .0001
      PLL2 = PLL(2) - .0001
      J = 0
      KM = 0
 2300 IF (.NOT. (J.LT.N)) GO TO 2310
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      IF (MASK(KM).EQ.0) THEN
         GO TO 2300

      END IF

      DO 2320 K = 1,JN
         IF (IAND(MASK(KM),IZBIT(K)).EQ.0) THEN
            GO TO 2320

         END IF

         IF (PX(K+J1).LT.PLL1 .OR. PX(K+J1).GT.PUR1) THEN
            MASK(KM) = IAND(MASK(KM),NOT(IZBIT(K)))
            GO TO 2320

         END IF

         IF (PY(K+J1).LT.PLL2 .OR. PY(K+J1).GT.PUR2) THEN
            MASK(KM) = IAND(MASK(KM),NOT(IZBIT(K)))
            GO TO 2320

         END IF

 2320 CONTINUE
      GO TO 2300

 2310 CONTINUE
      RETURN

      END
      SUBROUTINE PLTVWV(PLL,PUR,N,MASK,PX,PY,QX,QY)
      REAL PLL(2)
      REAL PUR(2)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL QX(*)
      REAL QY(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      PUR1 = PUR(1) + .0001
      PUR2 = PUR(2) + .0001
      PLL1 = PLL(1) - .0001
      PLL2 = PLL(2) - .0001
      DX = PUR1 - PLL1
      DY = PUR2 - PLL2
      J = 0
      KM = 0
 2340 IF (.NOT. (J.LT.N)) GO TO 2350
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2340

      END IF

      DO 2360 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2360

         END IF

         X1 = PX(K+J1)
         Y1 = PY(K+J1)
         X2 = QX(K+J1)
         Y2 = QY(K+J1)
         FP = X1 - PLL1
         FQ = X2 - PLL1
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))
            GO TO 2360

         END IF

         IF (FP.GT.DX .AND. FQ.GT.DX) THEN
            M = IAND(M,NOT(JB))
            GO TO 2360

         END IF

         DF = FQ - FP
         IF (DF.GT.0.) THEN
            TN = (Y2-Y1)/DF
            IF (FP.LT.0.) THEN
               X1 = PLL1
               Y1 = Y1 - FP*TN
            END IF

            IF (FQ.GT.DX) THEN
               X2 = PUR1
               Y2 = Y2 + (DX-FQ)*TN
            END IF

         ELSE IF (DF.LT.0.) THEN
            TN = (Y2-Y1)/DF
            IF (FQ.LT.0.) THEN
               X2 = PLL1
               Y2 = Y2 - FQ*TN
            END IF

            IF (FP.GT.DX) THEN
               X1 = PUR1
               Y1 = Y1 + (DX-FP)*TN
            END IF

         END IF

         FP = Y1 - PLL2
         FQ = Y2 - PLL2
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))
            GO TO 2360

         END IF

         IF (FP.GT.DY .AND. FQ.GT.DY) THEN
            M = IAND(M,NOT(JB))
            GO TO 2360

         END IF

         DF = FQ - FP
         IF (DF.GT.0.) THEN
            TN = (X2-X1)/DF
            IF (FP.LT.0.) THEN
               Y1 = PLL2
               X1 = X1 - FP*TN
            END IF

            IF (FQ.GT.DY) THEN
               Y2 = PUR2
               X2 = X2 + (DY-FQ)*TN
            END IF

         ELSE IF (DF.LT.0.) THEN
            TN = (X2-X1)/DF
            IF (FQ.LT.0.) THEN
               Y2 = PLL2
               X2 = X2 - FQ*TN
            END IF

            IF (FP.GT.DY) THEN
               Y1 = PUR2
               X1 = X1 + (DY-FP)*TN
            END IF

         END IF

         PX(K+J1) = X1
         PY(K+J1) = Y1
         QX(K+J1) = X2
         QY(K+J1) = Y2
 2360 CONTINUE
      MASK(KM) = M
      GO TO 2340

 2350 CONTINUE
      RETURN

      END
      SUBROUTINE PLTZCP(ZNEAR,ZFAR,N,MASK,PZ)
      INTEGER N
      INTEGER MASK(*)
      REAL PZ(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      J = 0
      KM = 0
 2380 IF (.NOT. (J.LT.N)) GO TO 2390
      JN = MIN(N-J,32)
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2380

      END IF

      DO 2400 I = 1,JN
         JB = IZBIT(I)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2400

         END IF

         IF ((PZ(I).LT.ZNEAR) .OR. (PZ(I).GT.ZFAR)) THEN
            M = IAND(M,NOT(JB))
            GO TO 2400

         END IF

         MASK(KM) = M
 2400 CONTINUE
      GO TO 2380

 2390 CONTINUE
      RETURN

      END
      SUBROUTINE PLTZCV(ZNEAR,ZFAR,N,MASK,PX,PY,PZ,QX,QY,QZ)
      INTEGER N
      INTEGER MASK(*)
      REAL PX(*)
      REAL PY(*)
      REAL PZ(*)
      REAL QX(*)
      REAL QY(*)
      REAL QZ(*)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      J = 0
      KM = 0
      DZ = ZFAR - ZNEAR
 2420 IF (.NOT. (J.LT.N)) GO TO 2430
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      KM = KM + 1
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2420

      END IF

      DO 2440 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).EQ.0) THEN
            GO TO 2440

         END IF

         X1 = PX(K+J1)
         Y1 = PY(K+J1)
         Z1 = PZ(K+J1)
         X2 = QX(K+J1)
         Y2 = QY(K+J1)
         Z2 = QZ(K+J1)
         FP = Z1 - ZNEAR
         FQ = Z2 - ZNEAR
         IF (FP.LT.0. .AND. FQ.LT.0.) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         IF (FP.GT.DZ .AND. FQ.GT.DZ) THEN
            M = IAND(M,NOT(JB))
            GO TO 2440

         END IF

         DF = FQ - FP
         IF (DF.GT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FP.LT.0.) THEN
               Z1 = ZNEAR
               X1 = X1 - FP*TN
               Y1 = Y1 - FP*SN
            END IF

            IF (FQ.GT.DZ) THEN
               Z2 = ZFAR
               X2 = X2 + (DZ-FQ)*TN
               Y2 = Y2 + (DZ-FQ)*SN
            END IF

         ELSE IF (DF.LT.0.) THEN
            TN = (X2-X1)/DF
            SN = (Y2-Y1)/DF
            IF (FQ.LT.0.) THEN
               Z2 = ZNEAR
               X2 = X2 - FQ*TN
               Y2 = Y2 - FQ*SN
            END IF

            IF (FP.GT.DZ) THEN
               Z1 = ZFAR
               X1 = X1 + (DZ-FP)*TN
               Y1 = Y1 + (DZ-FP)*SN
            END IF

         END IF

         PX(K+J1) = X1
         PY(K+J1) = Y1
         PZ(K+J1) = Z1
         QX(K+J1) = X2
         QY(K+J1) = Y2
         QZ(K+J1) = Z2
 2440 CONTINUE
      MASK(KM) = M
      GO TO 2420

 2430 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION PLTCNM(VALUE,CLR)
      CHARACTER*(*) CLR

      PLTCNM = .TRUE.
      IF (VALUE.EQ.0.) THEN
         CLR = 'BLACK'

      ELSE IF (VALUE.EQ.1.) THEN
         CLR = 'RED'

      ELSE IF (VALUE.EQ.2.) THEN
         CLR = 'GREEN'

      ELSE IF (VALUE.EQ.3.) THEN
         CLR = 'YELLOW'

      ELSE IF (VALUE.EQ.4.) THEN
         CLR = 'BLUE'

      ELSE IF (VALUE.EQ.6.) THEN
         CLR = 'CYAN'

      ELSE IF (VALUE.EQ.5.) THEN
         CLR = 'MAGENTA'

      ELSE IF (VALUE.EQ.7.) THEN
         CLR = 'WHITE'

      ELSE IF (VALUE.EQ.8.) THEN
         CLR = 'GRAY'

      ELSE IF (VALUE.EQ.10.) THEN
         CLR = 'DKGRAY'

      ELSE IF (VALUE.EQ.9.) THEN
         CLR = 'LTGRAY'

      ELSE IF (VALUE.EQ.12.) THEN
         CLR = 'LIME'

      ELSE IF (VALUE.EQ.11.) THEN
         CLR = 'PINK'

      ELSE IF (VALUE.EQ.15.) THEN
         CLR = 'ORANGE'

      ELSE IF (VALUE.EQ.14.) THEN
         CLR = 'VIOLET'

      ELSE IF (VALUE.EQ.13.) THEN
         CLR = 'LTBLUE'

      ELSE
         PLTCNM = .FALSE.
         CLR = 'UNKNOWN'
      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTICL(CLR,VAL)
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
      CHARACTER*(*) CLR
      CHARACTER*16 ICOL
      REAL VAL
      LOGICAL CHRCMP

      PLTICL = .TRUE.
      icol = '                '
      CALL CHRUP(CLR,ICOL)
      IF (CHRCMP(ICOL,'BLA','CK')) THEN
         VAL = 0.

      ELSE IF (CHRCMP(ICOL,'RED',' ')) THEN
         VAL = 1.

      ELSE IF (CHRCMP(ICOL,'GRE','EN')) THEN
         VAL = 2.

      ELSE IF (CHRCMP(ICOL,'YEL','LOW')) THEN
         VAL = 3.

      ELSE IF (CHRCMP(ICOL,'BLU','E')) THEN
         VAL = 4.

      ELSE IF (CHRCMP(ICOL,'MAG','ENTA')) THEN
         VAL = 5.

      ELSE IF (CHRCMP(ICOL,'CYA','N')) THEN
         VAL = 6.

      ELSE IF (CHRCMP(ICOL,'WHI','TE')) THEN
         VAL = 7.

      ELSE IF (CHRCMP(ICOL,'P1',' ')) THEN
         VAL = 8.

      ELSE IF (CHRCMP(ICOL,'P2',' ')) THEN
         VAL = 9.

      ELSE IF (CHRCMP(ICOL,'P3',' ')) THEN
         VAL = 10.

      ELSE IF (CHRCMP(ICOL,'P4',' ')) THEN
         VAL = 11.

      ELSE IF (CHRCMP(ICOL,'P5',' ')) THEN
         VAL = 12.

      ELSE IF (CHRCMP(ICOL,'P6',' ')) THEN
         VAL = 13.

      ELSE IF (CHRCMP(ICOL,'P7',' ')) THEN
         VAL = 14.

      ELSE IF (CHRCMP(ICOL,'P8',' ')) THEN
         VAL = 15.

      ELSE IF (CHRCMP(ICOL,'GR','AY')) THEN
         VAL = 8.

      ELSE IF (CHRCMP(ICOL,'LTG','RAY')) THEN
         VAL = 9.

      ELSE IF (CHRCMP(ICOL,'DKG','RAY')) THEN
         VAL = 10.

      ELSE IF (CHRCMP(ICOL,'PI','NK')) THEN
         VAL = 11.

      ELSE IF (CHRCMP(ICOL,'LI','ME')) THEN
         VAL = 12.

      ELSE IF (CHRCMP(ICOL,'LTB','LUE')) THEN
         VAL = 13.

      ELSE IF (CHRCMP(ICOL,'VIO','LET')) THEN
         VAL = 14.

      ELSE IF (CHRCMP(ICOL,'OR','ANGE')) THEN
         VAL = 15.

      ELSE
         PLTICL = .FALSE.
         VAL = 7
      END IF

      RETURN

      END
      SUBROUTINE PLTIQC(ICOLOR,R,G,B)
      DIMENSION CARRAY(3)

      CALL VDIQCO(1,ICOLOR,CARRAY,0)
      R = CARRAY(1)
      G = CARRAY(2)
      B = CARRAY(3)
      RETURN

      END
      SUBROUTINE PLTISP(SV,CV)
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

      IF (SV.LT.0.) THEN
         SV = 0.
      END IF

      IF (SV.GT.1.) THEN
         SV = 1.
      END IF

      CV = COLP(2) + REAL(NINT(SV*COLP(3)))
      IF (CV.GT.DEVCAP(4)-1.) THEN
         CV = DEVCAP(4) - 1.
      END IF

      RETURN

      END
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
      SUBROUTINE PLTSPC(S1,RED1,GREEN1,BLUE1,S2,RED2,GREEN2,BLUE2)
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
      DIMENSION DCOLOR(3,256),ICOL(256)

      IF (S1.LT.0. .OR. S1.GE.1.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *               'Starting value of color spectrum is out of range',
     *               2)
         RETURN

      END IF

      IF (S2.LE.0. .OR. S2.GT.1.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *               'Ending value of color spectrum is out of range',2)
         RETURN

      END IF

      IF (S2.LE.S1) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *'Starting value of color spectrum must be less than the ending val
     *ue',2)
         RETURN

      END IF

      IF (COLP(3).EQ.0.) THEN
         RETURN

      END IF

      I1 = INT(COLP(2)) + NINT(S1*COLP(3))
      I2 = INT(COLP(2)) + NINT(S2*COLP(3)-1.)
      DS = I2 - I1 + 1
      DR = (RED2-RED1)/DS
      DG = (GREEN2-GREEN1)/DS
      DB = (BLUE2-BLUE1)/DS
      ICOL(1) = I1
      DCOLOR(1,1) = RED1
      DCOLOR(2,1) = GREEN1
      DCOLOR(3,1) = BLUE1
      DO 2000 I = 2,DS
         ICOL(I) = ICOL(I-1) + 1
         DCOLOR(1,I) = DCOLOR(1,I-1) + DR
         DCOLOR(1,I) = MAX(0.,DCOLOR(1,I))
         DCOLOR(1,I) = MIN(1.,DCOLOR(1,I))
         DCOLOR(2,I) = DCOLOR(2,I-1) + DG
         DCOLOR(2,I) = MAX(0.,DCOLOR(2,I))
         DCOLOR(2,I) = MIN(1.,DCOLOR(2,I))
         DCOLOR(3,I) = DCOLOR(3,I-1) + DB
         DCOLOR(3,I) = MAX(0.,DCOLOR(3,I))
         DCOLOR(3,I) = MIN(1.,DCOLOR(3,I))
 2000 CONTINUE
      CALL VDSTCO(INT(DS),ICOL,DCOLOR,0)
      RETURN

      END
      SUBROUTINE PLTXTH(X,Y,TEXT)
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
      DIMENSION INCHAR(133)
      CHARACTER*(*) TEXT

      CALL PLTRIM(TEXT,L)
      IF (L.LE.0) THEN
         RETURN

      END IF

      NC = 0
      DO 2000 I = 1,L
         IC = ICHAR(TEXT(I:I))
         IF (IC.NE.0) THEN
            NC = NC + 1
            INCHAR(NC) = IC
         END IF

 2000 CONTINUE
      IF (NC.EQ.0) THEN
         RETURN

      END IF

      CALL PLTMOV(X,Y)
      CALL VDTEXT(NC,INCHAR)
      CALL VDIQCP(TEXTP(4),TEXTP(5))
      CALL PLTMOV(X,Y)
      TEXTP(6) = X
      TEXTP(7) = Y
      RETURN

      END
      SUBROUTINE PLTXTS(X,Y,TEXT)
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
      LOGICAL CPUIFC
      CHARACTER*(*) TEXT
      CHARACTER*1 LASCHR,ESCCHR
      CHARACTER*20 ESC
      INTEGER ASCII
      LOGICAL STATUS,CHRCI
      DATA ESCCHR/'\\'/

      IFONT = 1
      CALL PLTSVV
      CALL PLTSTV(1,1.)
      CALL PLTSTV(2,TEXTP(37))
      CALL PLTRIM(TEXT,NCHAR)
      TEXTP(18) = X
      TEXTP(19) = Y
      IFLAG = 0
      YSAVE = Y
      YCHRSZ = TEXTP(1)
      LASCHR = 'M'
      I = 1
 2020 IF (.NOT. (I.LE.NCHAR)) GO TO 2040
      ASCII = ICHAR(TEXT(I:I))
      IF (ASCII.LT.1 .OR. ASCII.GT.126) THEN
         CALL CHRIC(ASCII,LASCHR,LI)
         CALL PLTFLU
         CALL SIORPT('PLTXTS','Invalid character "'//LASCHR(1:LI)//
     *               '" in text string; rest of string ignored',2)
         RETURN

      END IF

      IF (ASCII.EQ.ICHAR(ESCCHR) .AND. TEXT(I+1:I+1).EQ.ESCCHR) THEN
         I = I + 1

      ELSE IF (ASCII.EQ.ICHAR(ESCCHR)) THEN
         CALL PLTESC(TEXT,I,ESC)
         CALL CHRUP(ESC,ESC)
         IF (ESC.EQ.'^') THEN
            IF (IFLAG.NE.1) THEN
               IF (IFLAG.EQ.-1) THEN
                  CALL PLTNOR(YBUMP,YCHRSZ)
               END IF

               CALL PLTSUP(YBUMP,YCHRSZ)
               IFLAG = 1
               GO TO 2030

            END IF

         ELSE IF (ESC.EQ.'_') THEN
            IF (IFLAG.NE.-1) THEN
               IF (IFLAG.EQ.1) THEN
                  CALL PLTNOR(YBUMP,YCHRSZ)
               END IF

               CALL PLTSUB(YBUMP,YCHRSZ)
            END IF

            IFLAG = -1
            GO TO 2030

         ELSE IF (ESC.EQ.'-') THEN
            IF (IFLAG.NE.0) THEN
               CALL PLTNOR(YBUMP,YCHRSZ)
            END IF

            IFLAG = 0
            GO TO 2030

         ELSE IF (ESC.EQ.'CLO') THEN
            ASCII = 4

         ELSE IF (ESC.EQ.'CSQ') THEN
            ASCII = 5

         ELSE IF (ESC.EQ.'CDI') THEN
            ASCII = 6

         ELSE IF (ESC.EQ.'CCS') THEN
            ASCII = 7

         ELSE IF (ESC.EQ.'CX') THEN
            ASCII = 8

         ELSE IF (ESC.EQ.'CTR') THEN
            ASCII = 9

         ELSE IF (ESC.EQ.'CCI') THEN
            ASCII = 10

         ELSE IF (ESC.EQ.'CDO') THEN
            ASCII = 11

         ELSE IF (ESC.EQ.'LO') THEN
            ASCII = 12

         ELSE IF (ESC.EQ.'SQ') THEN
            ASCII = 13

         ELSE IF (ESC.EQ.'DI') THEN
            ASCII = 14

         ELSE IF (ESC.EQ.'CS') THEN
            ASCII = 15

         ELSE IF (ESC.EQ.'X') THEN
            ASCII = 16

         ELSE IF (ESC.EQ.'TR') THEN
            ASCII = 17

         ELSE IF (ESC.EQ.'CI') THEN
            ASCII = 18

         ELSE IF (ESC.EQ.'DO') THEN
            ASCII = 19

         ELSE IF (ESC.EQ.'PLUSMIN') THEN
            ASCII = 20

         ELSE IF (ESC.EQ.'LEQ') THEN
            ASCII = 21

         ELSE IF (ESC.EQ.'GEQ') THEN
            ASCII = 22

         ELSE IF (ESC.EQ.'NEQ') THEN
            ASCII = 23

         ELSE IF (ESC.EQ.'PRIME') THEN
            ASCII = 24

         ELSE IF (ESC.EQ.'NLEQ') THEN
            ASCII = 25

         ELSE IF (ESC.EQ.'NGEQ') THEN
            ASCII = 26

         ELSE IF (ESC.EQ.'LL') THEN
            ASCII = 27

         ELSE IF (ESC.EQ.'GG') THEN
            ASCII = 28

         ELSE IF (ESC.EQ.'SUM') THEN
            ASCII = 29

         ELSE IF (ESC.EQ.'NLT') THEN
            ASCII = 30

         ELSE IF (ESC.EQ.'NGT') THEN
            ASCII = 31

         ELSE IF (ESC.EQ.'APPROX') THEN
            ASCII = 127

         ELSE IF (ESC.EQ.'CR') THEN
            TEXTP(18) = X + TEXTP(30)*TEXTP(1)*TEXTP(29)
            GO TO 2030

         ELSE IF (ESC.EQ.'LF') THEN
            TEXTP(19) = YSAVE - TEXTP(30)*TEXTP(1)*TEXTP(28)
            YSAVE = TEXTP(19)
            GO TO 2030

         ELSE IF (ESC.EQ.'CL') THEN
            TEXTP(18) = X + TEXTP(30)*TEXTP(1)*TEXTP(29)
            TEXTP(19) = YSAVE - TEXTP(30)*TEXTP(1)*TEXTP(28)
            YSAVE = TEXTP(19)
            GO TO 2030

         ELSE IF (ESC.EQ.'ENG') THEN
            IFONT = 1
            GO TO 2030

         ELSE IF (ESC.EQ.'GR') THEN
            IFONT = 2
            GO TO 2030

         ELSE IF (ESC.EQ.'DDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,3.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'DLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,2.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'LDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,5.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'MDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,6.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'SDLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,4.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE IF (ESC.EQ.'SLINE') THEN
            CALL PLTSVV
            CALL PLTSTV(1,1.)
            CALL PLTDV2(TEXTP(14),1,0.,.5,4.,.5)
            TEXTP(18) = TEXTP(18) + 4.*TEXTP(28)*YCHRSZ*TEXTP(31)
            TEXTP(19) = TEXTP(19) + 4.*TEXTP(29)*YCHRSZ*TEXTP(31)
            LASCHR = 'M'
            CALL PLTREV
            GO TO 2030

         ELSE
            STATUS = CHRCI(ESC,IESC)
            IF (STATUS) THEN
               ASCII = IESC

            ELSE
               CALL PLTRIM(ESC,L)
               CALL PLTFLU
               CALL SIORPT('PLTXTS','Invalid escape sequence "'//
     *                     ESC(1:L)//'"; escape sequence ignored.',2)
               GO TO 2030

            END IF

         END IF

      END IF

      NOVECT = NVECT(ASCII,IFONT)
      J = 0
 2050 IF (.NOT. (J.LT.NOVECT)) GO TO 2060
      JN = MIN(32,NOVECT-J)
      CALL PLTDV2(TEXTP(14),JN,X0(IDEX(ASCII,IFONT)+J,IFONT),
     *            Y0(IDEX(ASCII,IFONT)+J,IFONT),
     *            X1(IDEX(ASCII,IFONT)+J,IFONT),
     *            Y1(IDEX(ASCII,IFONT)+J,IFONT))
      J = J + JN
      GO TO 2050

 2060 CONTINUE
      TEXTP(18) = TEXTP(18) + XSIZE(ASCII,IFONT)*TEXTP(28)*YCHRSZ*
     *            TEXTP(31)
      TEXTP(19) = TEXTP(19) + XSIZE(ASCII,IFONT)*TEXTP(29)*YCHRSZ*
     *            TEXTP(31)
      IF (I.LE.LEN(TEXT)) THEN
         LASCHR = TEXT(I:I)
      END IF

      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2040

      END IF

 2030 I = I + 1
      GO TO 2020

 2040 CONTINUE
      CALL PLTMOV(X,Y)
      TEXTP(8) = TEXTP(18)
      TEXTP(9) = TEXTP(19)
      TEXTP(10) = X
      TEXTP(11) = YSAVE
      IF (IFLAG.NE.0) THEN
         DO 2070 I = 14,17
            TEXTP(I) = TEXTP(I)/TEXTP(32)
 2070    CONTINUE
      END IF

      CALL PLTREV
      RETURN

      END
      SUBROUTINE PLTESC(TEXT,I,ESC)
      CHARACTER*(*) TEXT,ESC

      IP = I + 1
      LT = LEN(TEXT)
      IF (TEXT(IP:IP).EQ.'^') THEN
         ESC = '^'
         I = I + 1
         RETURN

      ELSE IF (TEXT(IP:IP).EQ.'-') THEN
         ESC = '-'
         I = I + 1
         RETURN

      ELSE IF (TEXT(IP:IP).EQ.'_') THEN
         ESC = '_'
         I = I + 1
         RETURN

      END IF

      DO 2090 J = IP,LT
         IF (TEXT(J:J).EQ.' ' .OR. TEXT(J:J).EQ.'\\') THEN
            GO TO 2100

         END IF

 2090 CONTINUE
 2100 CONTINUE
      ESC = TEXT(IP:J-1)
      I = J
      IF (I.GE.LT) THEN
         RETURN

      END IF

      IF (TEXT(J:J).EQ.'\\') THEN
         I = J - 1
      END IF

      RETURN

      END
      SUBROUTINE PLTSUP(YBUMP,YCHRSZ)
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

      YBUMP = TEXTP(1)*TEXTP(33)
      TEXTP(18) = TEXTP(18) - YBUMP*TEXTP(29)
      TEXTP(19) = TEXTP(19) + YBUMP*TEXTP(28)
      DO 2110 J = 14,17
         TEXTP(J) = TEXTP(J)*TEXTP(32)
 2110 CONTINUE
      YCHRSZ = TEXTP(1)*TEXTP(32)
      RETURN

      END
      SUBROUTINE PLTSUB(YBUMP,YCHRSZ)
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

      YBUMP = -TEXTP(1)*TEXTP(34)
      TEXTP(18) = TEXTP(18) - YBUMP*TEXTP(29)
      TEXTP(19) = TEXTP(19) + YBUMP*TEXTP(28)
      DO 2130 J = 14,17
         TEXTP(J) = TEXTP(J)*TEXTP(32)
 2130 CONTINUE
      YCHRSZ = TEXTP(1)*TEXTP(32)
      RETURN

      END
      SUBROUTINE PLTNOR(YBUMP,YCHRSZ)
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

      DO 2150 J = 14,17
         TEXTP(J) = TEXTP(J)/TEXTP(32)
 2150 CONTINUE
      YCHRSZ = TEXTP(1)
      TEXTP(18) = TEXTP(18) + YBUMP*TEXTP(29)
      TEXTP(19) = TEXTP(19) - YBUMP*TEXTP(28)
      YBUMP = 0.
      RETURN

      END
      SUBROUTINE PLTFNT(FILENM)
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
      CHARACTER*(*) FILENM
      CHARACTER*64 LOCFIL
      CHARACTER*132 TLINE
      CHARACTER*10 IERROR
      INTEGER CHRLEN

      CALL CPUNAL(IUNIT)
      LOCFIL = FILENM
      L = CHRLEN(LOCFIL)
      IF (L.EQ.0) THEN
         RETURN

      END IF

      OPEN (UNIT=IUNIT,FILE=LOCFIL(1:L),FORM='unformatted',STATUS='old',
     *     ERR=10,IOSTAT=IOS)
      READ (IUNIT) NOVECT,TEXTP(39),TEXTP(40)
      IF (NOVECT.GT.2300) THEN
         CALL PLTFLU
         TLINE = 'Too many vectors in font file "'//LOCFIL(1:L)//
     *           '"; no font read in.'
         CALL SIORPT('PLTFNT',TLINE,2)
         RETURN

      END IF

      READ (IUNIT,ERR=30,IOSTAT=IOS) (IDEX(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (NVECT(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (XSIZE(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (YSIZE(I,1),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X0(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y0(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X1(I,1),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y1(I,1),I=1,NOVECT)
      READ (IUNIT) NOVECT,TEXTP(39),TEXTP(40)
      IF (NOVECT.GT.2300) THEN
         CALL PLTFLU
         TLINE = 'Too many vectors in font file "'//LOCFIL(1:L)//
     *           '"; no font read in.'
         CALL SIORPT('PLTFNT',TLINE,2)
         RETURN

      END IF

      READ (IUNIT,ERR=30,IOSTAT=IOS) (IDEX(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (NVECT(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (XSIZE(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (YSIZE(I,2),I=1,200)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X0(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y0(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (X1(I,2),I=1,NOVECT)
      READ (IUNIT,ERR=30,IOSTAT=IOS) (Y1(I,2),I=1,NOVECT)
      CLOSE (IUNIT)
      CALL CPUNDE(IUNIT)
      RETURN

   30 CONTINUE
      CALL CHRIC(IOS,IERROR,LI)
      CALL PLTFLU
      TLINE = 'I/O error '//IERROR(1:LI)//' in reading font file "'//
     *        LOCFIL(1:L)//'"; no font read in.'
      CALL SIORPT('PLTFNT',TLINE,2)
      RETURN

   10 CONTINUE
      CALL CHRIC(IOS,IERROR,LI)
      CALL PLTFLU
      TLINE = 'Open error '//IERROR(1:LI)//' in reading font file "'//
     *        LOCFIL(1:L)//'"; no font read in.'
      CALL SIORPT('PLTFNT',TLINE,2)
      RETURN

      END
      SUBROUTINE PLTSBM(N,MASK,X,Y,SYMB)
      DIMENSION X(*),Y(*),MASK(*)
      CHARACTER*(*) SYMB
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      J = 0
      KM = 0
 2170 IF (.NOT. (J.LT.N)) GO TO 2180
      JN = MIN(N-J,32)
      KM = KM + 1
      M = MASK(KM)
      J1 = J
      J = J + JN
      IF (M.EQ.0) THEN
         GO TO 2170

      END IF

      DO 2190 K = 1,JN
         IF (IAND(M,IZBIT(K)).NE.0) THEN
            CALL PLTXTS(X(J1+K),Y(J1+K),SYMB)
         END IF

 2190 CONTINUE
      GO TO 2170

 2180 CONTINUE
      RETURN

      END
      SUBROUTINE PLTSYM(X,Y,S)
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
      INTEGER S
      REAL X,Y
      CHARACTER CH

      IF (S.LT.32 .OR. S.GT.126) THEN
         RETURN

      END IF

      CH = CHAR(S)
      CALL PLTXTH(X-TEXTP(36)* (5./7.)/2.,Y-TEXTP(35)/2.,CH)
      RETURN

      END
      SUBROUTINE PLTXTC(X,Y,LINE)
      CHARACTER*(*) LINE

      CALL PLTRIM(LINE,L)
      IF (L.LE.0) THEN
         RETURN

      END IF

      CALL PLTXSL(LINE(1:L),XL)
      CALL PLTXTS(X-XL/2.,Y,LINE(1:L))
      RETURN

      END
      SUBROUTINE PLTXHE(X,Y)
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

      X = TEXTP(4)
      Y = TEXTP(5)
      RETURN

      END
      SUBROUTINE PLTXHL(CHARST,LENGTH)
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
      CHARACTER*(*) CHARST
      REAL LENGTH

      LENGTH = REAL(LEN(CHARST))*TEXTP(36)
      RETURN

      END
      SUBROUTINE PLTXHN(X,Y)
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

      X = TEXTP(6)
      Y = TEXTP(7) - TEXTP(35)
      RETURN

      END
      SUBROUTINE PLTXSE(X,Y)
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

      X = TEXTP(8)
      Y = TEXTP(9)
      RETURN

      END
      SUBROUTINE PLTXSL(CHARST,LENGTH)
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
      CHARACTER*(*) CHARST
      CHARACTER ESC*20,ESCCHR*1
      REAL LENGTH,LENGT1
      INTEGER ASCII
      LOGICAL STATUS,CHRCI
      DATA ESCCHR/'\\'/

      LENGTH = 0.
      LENGT1 = 0.
      YSZE = TEXTP(1)
      CALL PLTRIM(CHARST,NUM)
      MODE = 0
      IFONT = 1
      I = 1
 2210 IF (.NOT. (I.LE.NUM)) GO TO 2230
      ASCII = ICHAR(CHARST(I:I))
      IF (ASCII.LT.1 .OR. ASCII.GT.126) THEN
         CALL CHRIC(ASCII,ESCCHR,LI)
         CALL PLTFLU
         CALL SIORPT('PLTXSL','Invalid character "'//ESCCHR(1:LI)//
     *               '" in text string; rest of string ignored',2)
         RETURN

      END IF

      RLINE = 0.
      IF (ASCII.EQ.ICHAR(ESCCHR) .AND. CHARST(I+1:I+1).EQ.ESCCHR) THEN
         I = I + 1

      ELSE IF (ASCII.EQ.ICHAR(ESCCHR)) THEN
         CALL PLTESC(CHARST,I,ESC)
         CALL CHRUP(ESC,ESC)
         IF (ESC.EQ.'^' .OR. ESC.EQ.'_') THEN
            IF (MODE.EQ.0) THEN
               YSZE = TEXTP(1)*TEXTP(32)
               MODE = 1
            END IF

            GO TO 2220

         ELSE IF (ESC.EQ.'-') THEN
            YSZE = TEXTP(1)
            MODE = 0
            GO TO 2220

         ELSE IF (ESC.EQ.'CLO') THEN
            ASCII = 4

         ELSE IF (ESC.EQ.'CSQ') THEN
            ASCII = 5

         ELSE IF (ESC.EQ.'CDI') THEN
            ASCII = 6

         ELSE IF (ESC.EQ.'CCS') THEN
            ASCII = 7

         ELSE IF (ESC.EQ.'CX') THEN
            ASCII = 8

         ELSE IF (ESC.EQ.'CTR') THEN
            ASCII = 9

         ELSE IF (ESC.EQ.'CCI') THEN
            ASCII = 10

         ELSE IF (ESC.EQ.'CDO') THEN
            ASCII = 11

         ELSE IF (ESC.EQ.'LO') THEN
            ASCII = 12

         ELSE IF (ESC.EQ.'SQ') THEN
            ASCII = 13

         ELSE IF (ESC.EQ.'DI') THEN
            ASCII = 14

         ELSE IF (ESC.EQ.'CS') THEN
            ASCII = 15

         ELSE IF (ESC.EQ.'X') THEN
            ASCII = 16

         ELSE IF (ESC.EQ.'TR') THEN
            ASCII = 17

         ELSE IF (ESC.EQ.'CI') THEN
            ASCII = 18

         ELSE IF (ESC.EQ.'DO') THEN
            ASCII = 19

         ELSE IF (ESC.EQ.'PLUSMIN') THEN
            ASCII = 20

         ELSE IF (ESC.EQ.'LEQ') THEN
            ASCII = 21

         ELSE IF (ESC.EQ.'GEQ') THEN
            ASCII = 22

         ELSE IF (ESC.EQ.'NEQ') THEN
            ASCII = 23

         ELSE IF (ESC.EQ.'PRIME') THEN
            ASCII = 24

         ELSE IF (ESC.EQ.'NLEQ') THEN
            ASCII = 25

         ELSE IF (ESC.EQ.'NGEQ') THEN
            ASCII = 26

         ELSE IF (ESC.EQ.'LL') THEN
            ASCII = 27

         ELSE IF (ESC.EQ.'GG') THEN
            ASCII = 28

         ELSE IF (ESC.EQ.'SUM') THEN
            ASCII = 29

         ELSE IF (ESC.EQ.'NLT') THEN
            ASCII = 30

         ELSE IF (ESC.EQ.'NGT') THEN
            ASCII = 31

         ELSE IF (ESC.EQ.'APPROX') THEN
            ASCII = 127

         ELSE IF (ESC.EQ.'CR') THEN
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'LF') THEN
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'CL') THEN
            IF (LENGTH.GT.LENGT1) THEN
               LENGT1 = LENGTH
            END IF

            LENGTH = 0.
            GO TO 2220

         ELSE IF (ESC.EQ.'ENG') THEN
            IFONT = 1
            GO TO 2220

         ELSE IF (ESC.EQ.'GR') THEN
            IFONT = 2
            GO TO 2220

         ELSE IF (ESC.EQ.'DDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'DLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'LDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'MDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'SDLINE') THEN
            RLINE = 4.

         ELSE IF (ESC.EQ.'SLINE') THEN
            RLINE = 4.

         ELSE
            STATUS = CHRCI(ESC,IESC)
            IF (STATUS) THEN
               ASCII = IESC

            ELSE
               CALL PLTRIM(ESC,L)
               CALL PLTFLU
               CALL SIORPT('PLTXTS','Invalid escape sequence "'//
     *                     ESC(1:L)//'"; escape sequence ignored.',2)
               GO TO 2220

            END IF

         END IF

      END IF

      IF (RLINE.EQ.4.) THEN
         LENGTH = LENGTH + 4.*YSZE**TEXTP(31)

      ELSE
         LENGTH = LENGTH + XSIZE(ASCII,IFONT)*YSZE*TEXTP(31)
      END IF

 2220 I = I + 1
      GO TO 2210

 2230 CONTINUE
      IF (LENGT1.GT.LENGTH) THEN
         LENGTH = LENGT1
      END IF

      IF (NUM.EQ.1 .AND. RLINE.NE.4. .AND. ASCII.GT.32) THEN
         TEMP = TEXTP(39)/TEXTP(40)
         LENGTH = LENGTH - TEMP*YSZE*TEXTP(31)
      END IF

      RETURN

      END
      SUBROUTINE PLTXSN(X,Y)
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

      X = TEXTP(10) + TEXTP(30)*TEXTP(1)*TEXTP(29)
      Y = TEXTP(11) - TEXTP(30)*TEXTP(1)*TEXTP(28)
      RETURN

      END
      SUBROUTINE PLTRIM(LINE,L)
      CHARACTER*(*) LINE
      CHARACTER CH

      L1 = LEN(LINE)
 2240 IF (.NOT. (L1.GT.0)) GO TO 2260
      CH = LINE(L1:L1)
      IF (CH.EQ.' ') THEN
         GO TO 2250

      END IF

      GO TO 2260

 2250 L1 = L1 - 1
      GO TO 2240

 2260 CONTINUE
      L = L1
      RETURN

      END
      SUBROUTINE PLTCUR(X,Y,NUM)
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
      LOGICAL CPUIFC
      DIMENSION X(*),Y(*),X0T(32),Y0T(32),X1T(32),Y1T(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/
      DATA SYMRAT/120./

      IF (NUM.LE.0) THEN
         RETURN

      END IF

      CALL PLTSVV
      CALL PLTSVD
      CALL PLTSVT
      CALL PLTSTV(1,GRAPHP(5))
      CALL PLTSTD(1,GRAPHP(38))
      CALL PLTSTV(2,GRAPHP(63))
      IF (GRAPHP(6).EQ.1.) THEN
         CALL PLTSTT(3,0.)
         CALL PLTSTT(4,0.)
         SYMSIZ = (GRAPHP(3)+GRAPHP(4))*GRAPHP(46)/ (10.*SYMRAT)
         CALL PLTSTT(2,SYMSIZ)
         CALL PLTSTT(11,GRAPHP(66))
      END IF

      J = 0
      NP = 0
      IF (NUM.EQ.1) THEN
         IF (GRAPHP(6).EQ.1. .AND. SYMSIZ.GT.0.0) THEN
            IF (GRAPHP(21).EQ.1.) THEN
               X0T(1) = X(1)
               Y0T(1) = Y(1)

            ELSE IF (GRAPHP(21).EQ.2.) THEN
               IF (X(1).LE.0.) THEN
                  RETURN

               END IF

               X0T(1) = ALOG10(X(1))
               Y0T(1) = Y(1)

            ELSE IF (GRAPHP(21).EQ.3.) THEN
               X0T(1) = X(1)
               IF (Y(1).LE.0.) THEN
                  RETURN

               END IF

               Y0T(1) = ALOG10(Y(1))

            ELSE IF (GRAPHP(21).EQ.4.) THEN
               IF (X(1).LE.0.) THEN
                  RETURN

               END IF

               IF (Y(1).LE.0.) THEN
                  RETURN

               END IF

               X0T(1) = ALOG10(X(1))
               Y0T(1) = ALOG10(Y(1))
            END IF

            MASKS = -1
            CALL PLTMP2(GRAPHP(7),1,MASKS,X0T,Y0T,X1T,Y1T)
            JB = IZBIT(1)
            IF (IAND(MASKS,JB).NE.0) THEN
               CALL PLTSTD(1,GRAPHP(75))
               CALL PLTXTS(X1T(1),Y1T(1),CHAR(INT(GRAPHP(47))))
            END IF

         END IF

      ELSE
         IF (GRAPHP(21).EQ.1.) THEN
            XSAV = X(1)
            YSAV = Y(1)

         ELSE IF (GRAPHP(21).EQ.2.) THEN
            DO 2000 I = 1,NUM
               IF (X(I).LE.0.) THEN
                  GO TO 2000

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2010

 2000       CONTINUE
 2010       CONTINUE

         ELSE IF (GRAPHP(21).EQ.3.) THEN
            DO 2020 I = 1,NUM
               IF (Y(I).LE.0.) THEN
                  GO TO 2020

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2030

 2020       CONTINUE
 2030       CONTINUE

         ELSE IF (GRAPHP(21).EQ.4.) THEN
            DO 2040 I = 1,NUM
               IF (X(I).LE.0. .OR. Y(I).LE.0.) THEN
                  GO TO 2040

               END IF

               XSAV = X(I)
               YSAV = Y(I)
               GO TO 2050

 2040       CONTINUE
 2050       CONTINUE
         END IF

 2060    IF (.NOT. (J.LT.NUM-1)) GO TO 2070
         K = MIN(NUM-J,32)
         IF (GRAPHP(21).EQ.1.) THEN
            NV = 0
            DO 2080 I = 1,K - 1
               NV = NV + 1
               X0T(I) = XSAV
               Y0T(I) = YSAV
               X1T(I) = X(I+J+1)
               Y1T(I) = Y(I+J+1)
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2080       CONTINUE

         ELSE IF (GRAPHP(21).EQ.2.) THEN
            NV = 0
            DO 2100 I = 1,K - 1
               IF (X(J+I+1).LE.0.) THEN
                  GO TO 2100

               END IF

               NV = NV + 1
               X0T(NV) = ALOG10(XSAV)
               Y0T(NV) = YSAV
               X1T(NV) = ALOG10(X(J+I+1))
               Y1T(NV) = Y(J+I+1)
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2100       CONTINUE

         ELSE IF (GRAPHP(21).EQ.3.) THEN
            NV = 0
            DO 2120 I = 1,K - 1
               IF (Y(J+I+1).LE.0.) THEN
                  GO TO 2120

               END IF

               NV = NV + 1
               X0T(NV) = XSAV
               Y0T(NV) = ALOG10(YSAV)
               X1T(NV) = X(J+I+1)
               Y1T(NV) = ALOG10(Y(J+I+1))
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2120       CONTINUE

         ELSE IF (GRAPHP(21).EQ.4.) THEN
            NV = 0
            DO 2140 I = 1,K - 1
               IF (X(J+I+1).LE.0.) THEN
                  GO TO 2140

               END IF

               IF (Y(J+I+1).LE.0.) THEN
                  GO TO 2140

               END IF

               NV = NV + 1
               X0T(NV) = ALOG10(XSAV)
               Y0T(NV) = ALOG10(YSAV)
               X1T(NV) = ALOG10(X(J+I+1))
               Y1T(NV) = ALOG10(Y(J+I+1))
               XSAV = X(I+J+1)
               YSAV = Y(I+J+1)
 2140       CONTINUE
         END IF

         CALL PLTDV2(GRAPHP(7),NV,X0T,Y0T,X1T,Y1T)
         IF (GRAPHP(6).EQ.1. .AND. SYMSIZ.GT.0.0) THEN
            CALL PLTSTD(1,GRAPHP(75))
            MASKS = -1
            XT = X1T(NV)
            YT = Y1T(NV)
            CALL PLTMP2(GRAPHP(7),NV,MASKS,X0T,Y0T,X1T,Y1T)
            DO 2160 L = 1,NV
               JB = IZBIT(L)
               IF (IAND(MASKS,JB).NE.0 .AND.
     *             MOD(L+NP+INT(GRAPHP(23))-1,INT(GRAPHP(23))).EQ.
     *             0) THEN
                  CALL PLTXTS(X1T(L),Y1T(L),CHAR(INT(GRAPHP(47))))
               END IF

 2160       CONTINUE
            CALL PLTSTD(1,GRAPHP(38))
         END IF

         NP = NP + NV
         J = J + K - 1
         IF (J+1.GE.NUM .AND. (GRAPHP(6).EQ.1..AND.SYMSIZ.GT.0.0)) THEN
            X0T(1) = XT
            Y0T(1) = YT
            CALL PLTSTD(1,GRAPHP(75))
            MASKS = -1
            CALL PLTMP2(GRAPHP(7),1,MASKS,X0T,Y0T,X1T,Y1T)
            JB = IZBIT(1)
            IF (IAND(MASKS,JB).NE.0 .AND.
     *          MOD(NP+INT(GRAPHP(23)),INT(GRAPHP(23))).EQ.0) THEN
               CALL PLTXTS(X1T(1),Y1T(1),CHAR(INT(GRAPHP(47))))
            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2070

         END IF

         GO TO 2060

 2070    CONTINUE
      END IF

      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
      LOGICAL FUNCTION PLTD2G(XD,YD,XG,YG)
      DIMENSION UMAP(14)

      CALL PLTGTG(27,UMAP)
      CALL PLTGTG(9,TYPE)
      PLTD2G = .FALSE.
      TOP = XD*UMAP(4) - YD*UMAP(3) - UMAP(5)*UMAP(4) + UMAP(6)*UMAP(3)
      BOTTOM = UMAP(1)*UMAP(4) - UMAP(2)*UMAP(3)
      XG = TOP/BOTTOM
      TOP = XD*UMAP(2) - YD*UMAP(1) - UMAP(5)*UMAP(2) + UMAP(6)*UMAP(1)
      BOTTOM = UMAP(3)*UMAP(2) - UMAP(4)*UMAP(1)
      YG = TOP/BOTTOM
      IF (TYPE.EQ.2.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG

      ELSE IF (TYPE.EQ.3.) THEN
         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         YG = 10.**YG

      ELSE IF (TYPE.EQ.4.) THEN
         IF (XD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         IF (YD.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTD2G',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XG = 10.**XG
         YG = 10.**YG
      END IF

      PLTD2G = .TRUE.
      RETURN

      END
      LOGICAL FUNCTION PLTG2D(XG,YG,XD,YD)
      DIMENSION UMAP(14)

      PLTG2D = .FALSE.
      CALL PLTGTG(27,UMAP)
      CALL PLTGTG(9,TYPE)
      IF (TYPE.EQ.1.) THEN
         XD = XG*UMAP(1) + YG*UMAP(3) + UMAP(5)
         YD = XG*UMAP(2) + YG*UMAP(4) + UMAP(6)

      ELSE IF (TYPE.EQ.2.) THEN
         IF (XG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         XD = ALOG10(XG)*UMAP(1) + YG*UMAP(3) + UMAP(5)
         YD = ALOG10(XG)*UMAP(2) + YG*UMAP(4) + UMAP(6)

      ELSE IF (TYPE.EQ.3.) THEN
         IF (YG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XD = XG*UMAP(1) + ALOG10(YG)*UMAP(3) + UMAP(5)
         YD = XG*UMAP(2) + ALOG10(YG)*UMAP(4) + UMAP(6)

      ELSE IF (TYPE.EQ.4.) THEN
         IF (XG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log x axis.'
     *                  ,2)
            RETURN

         END IF

         IF (YG.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTG2D',
     *                 'Cannot convert a negative number on log y axis.'
     *                  ,2)
            RETURN

         END IF

         XD = ALOG10(XG)*UMAP(1) + ALOG10(YG)*UMAP(3) + UMAP(5)
         YD = ALOG10(XG)*UMAP(2) + ALOG10(YG)*UMAP(4) + UMAP(6)
      END IF

      PLTG2D = .TRUE.
      RETURN

      END
      SUBROUTINE PLTGPH(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
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
      DIMENSION X(1),Y(1)

      IF (GRAPHP(21).EQ.1.) THEN
         CALL PLTNXY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)

      ELSE IF (GRAPHP(21).EQ.2.) THEN
         CALL PLTLGX(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)

      ELSE IF (GRAPHP(21).EQ.3.) THEN
         CALL PLTLGY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)

      ELSE IF (GRAPHP(21).EQ.4.) THEN
         CALL PLTLXY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
      END IF

      RETURN

      END
      SUBROUTINE PLTAXS(X,Y,XLENG,YLENG,TYPE,XMIN,XMAX,XSTART,NDEC,
     *                  INTER,MININT,LABEL,UNITS,EXP)
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
      LOGICAL CPUIFC
      INTEGER EXP
      REAL LENMAJ,LENMIN,MAJRAT,MINRAT,NUMRAT,NUMSIZ,LABSIZ,MINTIC,
     *     LABRAT,MAJTIC,LONNUM,LABSCA,INTER
      CHARACTER*(*) LABEL,UNITS,TYPE
      CHARACTER TTYPE*1,FORMA*32,LINE*15,TLABEL*132,EXPL*4,CINPUT*132,
     *          CNDEC*1
      LOGICAL PLTNER
      DATA MAJRAT/50./
      DATA MINRAT/2./
      DATA NUMRAT/50./
      DATA LABRAT/40./
      DATA FUDGE/.0025/

      IF (INTER.EQ.0.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTAXS',
     *               'Interval between major ticks must be nonzero.',2)
         RETURN

      END IF

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      CALL PLTSVT
      CALL PLTSVD
      CALL PLTSVV
      CALL PLTSTT(3,0.)
      CALL PLTSTT(4,0.)
      CALL PLTSTV(1,1.)
      CALL CHRIC(NDEC,CNDEC,L)
      FORMA = '(f15.'//CNDEC//')'
 2180 CONTINUE
      IF (TTYPE.EQ.'X') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(44)/5.))/NUMRAT
         LABSIZ = (LABSCA* (GRAPHP(45)/5.))/LABRAT
         YOFF = Y - NUMSIZ*1.8
         LENMAJ = XLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         MAJTIC = XLENG/ ((XMAX-XMIN)/INTER)
         IF (MININT.NE.0) THEN
            MINTIC = MAJTIC/MININT

         ELSE
            MINTIC = 0.
         END IF

         XLOW = XSTART
         IF (XSTART.NE.XMIN) THEN
            XLOW = XSTART - INTER
         END IF

         FMJTIC = X - (XMIN-XLOW)*MAJTIC/INTER
         XMAJ = FMJTIC
         XNUM = XSTART
         IF (PLTNER(FMJTIC,X) .AND. XSTART.NE.XMIN) THEN
            XNUM = XSTART - INTER
         END IF

         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTVCT(1,X,Y,X+XLENG,Y)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X,Y+YLENG,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,X) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=120) XNUM
  120       DO 2210 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2220

               END IF

 2210       CONTINUE
 2220       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(X-TLEN/2.,YOFF,LINE(J:ILAST))
            XNUM = XNUM + INTER
         END IF

 2230    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2250

         END IF

         DO 2260 I = 1,MININT - 1
            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2270

            END IF

            XNEW = XMAJ + FLOAT(I)*MINTIC
            IF (XNEW-FUDGE.LE.X) THEN
               GO TO 2260

            END IF

            IF (XNEW.GT.X+XLENG-FUDGE) THEN
               GO TO 2270

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,XNEW,Y,XNEW,Y+LENMIN)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XNEW,Y+YLENG,XNEW,Y+YLENG-LENMIN)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG-LENMIN)

               ELSE
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG)
               END IF

               CALL PLTSTV(1,1.)
            END IF

 2260    CONTINUE
 2270    CONTINUE
         XMAJ = XMAJ + MAJTIC
         IF (XMAJ.GT.X+XLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2250

         ELSE IF (XMAJ.GT.X+XLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2250

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,XMAJ,Y,XMAJ,Y+LENMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,XMAJ,Y+YLENG,XMAJ,Y+YLENG-LENMAJ)
         END IF

         IF (GRAPHP(73).NE.0. .OR. GRAPHP(35).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG-LENMAJ)

            ELSE
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=121) XNUM
  121       DO 2280 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2290

               END IF

 2280       CONTINUE
 2290       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE(J:ILAST))
         END IF

         XNUM = XNUM + INTER
 2240    GO TO 2230

 2250    CONTINUE
         IF (PLTNER(XMAJ,X+XLENG) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=122) XNUM
  122       DO 2300 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2310

               END IF

 2300       CONTINUE
 2310       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE(J:ILAST))
         END IF

         IF ((XMIN.LT.0..AND.XMAX.GT.0.) .AND.
     *       (GRAPHP(35).EQ.0..AND.GRAPHP(73).EQ.0.) .AND.
     *       GRAPHP(71).NE.0.) THEN
            CALL PLTSTV(1,GRAPHP(71))
            CALL PLTSTD(1,GRAPHP(72))
            CALL PLTSTV(2,GRAPHP(70))
            X0LINE = X - (XMIN*XLENG)/ (XMAX-XMIN)
            CALL PLTVCT(1,X0LINE,Y+LENMAJ,X0LINE,Y+YLENG-LENMAJ)
            CALL PLTSTV(1,1.)
         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (EXP.NE.0. .AND. NUMSIZ.GT.0.) THEN
               IF (LL.GT.0) THEN
                  LL = LL + 1
               END IF

               EXPL = ' '
               WRITE (EXPL,'(i4)',ERR=110) EXP
  110          DO 2320 I = 1,4
                  IF (EXPL(I:I).NE.' ') THEN
                     GO TO 2330

                  END IF

 2320          CONTINUE
 2330          CONTINUE
               TLABEL(LL+1:) = '( *10\^'//EXPL(I:)//'\- )'
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X + (XLENG-TLEN)/2.
               YLAB = YOFF - LABSIZ*2.
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2200

            END IF

         END IF

      ELSE IF (TTYPE.EQ.'Y') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(88)/5.))/NUMRAT
         LABSIZ = (LABSCA* (GRAPHP(89)/5.))/LABRAT
         LONNUM = 0.
         XOFF = NUMSIZ*.8
         YOFF = NUMSIZ/2.
         LENMAJ = YLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         MAJTIC = YLENG/ ((XMAX-XMIN)/INTER)
         IF (MININT.NE.0) THEN
            MINTIC = MAJTIC/MININT

         ELSE
            MINTIC = 0.
         END IF

         XLOW = XSTART
         IF (XSTART.NE.XMIN) THEN
            XLOW = XSTART - INTER
         END IF

         FMJTIC = Y - (XMIN-XLOW)*MAJTIC/INTER
         YMAJ = FMJTIC
         XNUM = XSTART
         IF (PLTNER(FMJTIC,Y) .AND. XSTART.NE.XMIN) THEN
            XNUM = XSTART - INTER
         END IF

         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTVCT(1,X,Y,X,Y+YLENG)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,Y,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,Y) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=130) XNUM
  130       DO 2340 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2350

               END IF

 2340       CONTINUE
 2350       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),Y-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

            XNUM = XNUM + INTER
         END IF

 2360    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2380

         END IF

         DO 2390 I = 1,MININT - 1
            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2400

            END IF

            YNEW = YMAJ + FLOAT(I)*MINTIC
            IF (YNEW-FUDGE.LE.Y) THEN
               GO TO 2390

            END IF

            IF (YNEW.GT.Y+YLENG-FUDGE) THEN
               GO TO 2400

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,X,YNEW,X+LENMIN,YNEW)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+XLENG,YNEW,X+XLENG-LENMIN,YNEW)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG-LENMIN,YNEW)

               ELSE
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG,YNEW)
               END IF

               CALL PLTSTV(1,1.)
            END IF

 2390    CONTINUE
 2400    CONTINUE
         YMAJ = YMAJ + MAJTIC
         IF (YMAJ.GT.Y+YLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2380

         ELSE IF (YMAJ.GT.Y+YLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2380

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,X,YMAJ,X+LENMAJ,YMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,YMAJ,X+XLENG-LENMAJ,YMAJ)
         END IF

         IF (GRAPHP(35).NE.0. .OR. GRAPHP(73).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG-LENMAJ,YMAJ)

            ELSE
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG,YMAJ)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=131) XNUM
  131       DO 2410 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2420

               END IF

 2410       CONTINUE
 2420       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),YMAJ-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

         END IF

         XNUM = XNUM + INTER
 2370    GO TO 2360

 2380    CONTINUE
         IF (PLTNER(YMAJ,Y+YLENG) .AND. NUMSIZ.GT.0.0) THEN
            LINE = ' '
            WRITE (LINE,FORMA,ERR=132) XNUM
  132       DO 2430 J = 1,15
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2440

               END IF

 2430       CONTINUE
 2440       CONTINUE
            ILAST = 15
            IF (NDEC.EQ.0) THEN
               ILAST = 14
            END IF

            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            CALL PLTXSL(LINE(J:ILAST),TLEN)
            IF (GRAPHP(92).EQ.1.) THEN
               LONNUM = NUMSIZ
               CALL PLTSTT(3,90.)
               CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE(J:ILAST))
               CALL PLTSTT(3,0.)

            ELSE
               CALL PLTXTS(X- (TLEN+XOFF),YMAJ-NUMSIZ/2.,LINE(J:ILAST))
               IF (TLEN.GT.LONNUM) THEN
                  LONNUM = TLEN
               END IF

            END IF

         END IF

         IF ((XMIN.LT.0..AND.XMAX.GT.0.) .AND.
     *       (GRAPHP(35).EQ.0..AND.GRAPHP(73).EQ.0.) .AND.
     *       GRAPHP(71).NE.0.) THEN
            CALL PLTSTV(1,GRAPHP(71))
            CALL PLTSTD(1,GRAPHP(72))
            CALL PLTSTV(2,GRAPHP(70))
            Y0LINE = Y - (XMIN*YLENG)/ (XMAX-XMIN)
            CALL PLTVCT(1,X+LENMAJ,Y0LINE,X+XLENG-LENMAJ,Y0LINE)
            CALL PLTSTV(1,1.)
         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (EXP.NE.0. .AND. NUMSIZ.GT.0.) THEN
               IF (LL.GT.0) THEN
                  LL = LL + 1
               END IF

               EXPL = ' '
               WRITE (EXPL,'(i4)',ERR=134) EXP
  134          DO 2450 I = 1,4
                  IF (EXPL(I:I).NE.' ') THEN
                     GO TO 2460

                  END IF

 2450          CONTINUE
 2460          CONTINUE
               TLABEL(LL+1:) = '( *10\^'//EXPL(I:)//'\- )'
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTSTT(3,90.)
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X - LONNUM - LABSIZ*1.4
               YLAB = Y + (YLENG-TLEN)/2.
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
               CALL PLTSTT(3,0.)
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2200

            END IF

         END IF

      ELSE
         CALL PLTFLU
         CALL SIORPT('PLTAXS','Invalid axis type - '//TTYPE,2)
      END IF

 2190 IF (.NOT. (.TRUE.)) GO TO 2180
 2200 CONTINUE
      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
      SUBROUTINE PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER,NINTER
      INTEGER EXP
      REAL NI
      DATA SMALL/1.E-4/,DESINT/5./

      DELTA = MAX - MIN
      TMAX = MAX
      TMIN = MIN
      IF (DELTA.NE.0) THEN
         EXP = NINT(ALOG10(ABS(DELTA))) - 1

      ELSE
         IF (MIN.NE.0) THEN
            EXP = NINT(ALOG10(ABS(MIN))) - 1
            EPS = .01*ABS(MIN)

         ELSE
            EXP = 0
            EPS = .1
         END IF

         TMAX = MAX + EPS
         TMIN = TMIN - EPS
      END IF

      TENEXP = 10.**EXP
      RMIN = 1.E20
      J = 1
      DO 2470 I = 1,5
         NINTER = DELTA/ (FLOAT(I)*TENEXP)
         TEMP = ABS(DESINT-NINTER)
         IF (TEMP.LT.RMIN) THEN
            J = I
            RMIN = TEMP
         END IF

 2470 CONTINUE
      INTER = FLOAT(J)
      IF (DELTA.EQ.0.) THEN
         INTER = 1.
      END IF

      IF (INTER.EQ.1.) THEN
         NMIN = 5

      ELSE IF (INTER.EQ.2.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.3.) THEN
         NMIN = 6

      ELSE IF (INTER.EQ.4.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.5.) THEN
         NMIN = 5
      END IF

      TENI = INTER*TENEXP
      IF (TMIN.GE.0.) THEN
         ADD = 0.

      ELSE
         ADD = -1.
      END IF

      RNI = TMIN/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      START = (NI+ADD)*TENI
      IF (TMAX.GE.0.) THEN
         ADD = 1.

      ELSE
         ADD = 0.
      END IF

      RNI = TMAX/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      REND = (NI+ADD)*TENI
      START = START/TENEXP
      REND = REND/TENEXP
      IF (REND.NE.0.) THEN
 2490    IF (.NOT. (ABS(REND).GT.10.)) GO TO 2500
         REND = REND/10.
         START = START/10.
         INTER = INTER/10.
         EXP = EXP + 1
         GO TO 2490

 2500    CONTINUE
      END IF

      IF (START.NE.0.) THEN
 2510    IF (.NOT. (ABS(START).LT.1.)) GO TO 2520
         REND = REND*10.
         START = START*10.
         INTER = INTER*10.
         EXP = EXP - 1
         GO TO 2510

 2520    CONTINUE
      END IF

      IF (START.EQ.0 .OR. TMIN.EQ.0) THEN
         RETURN

      END IF

      IF (START.EQ.REND) THEN
         REND = START + INTER
      END IF

      IF (ABS(START-REND).EQ.INTER) THEN
         NMIN = 10
      END IF

      RETURN

      END
      SUBROUTINE PLTINI(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER
      INTEGER EXP

      DELTA = MAX - MIN
      IF (DELTA.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTINI',
     *               'Maximum value must be greater than minimum value.'
     *               ,2)
         RETURN

      END IF

      CALL PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      TENEXP = 10.**EXP
      IEXP = NINT(ALOG10(ABS(INTER))) - 2
      SMALL = 10.**IEXP
      IF (ABS(START-MIN/TENEXP).GT.SMALL) THEN
         START = START + INTER
      END IF

      IF (ABS(REND-MAX/TENEXP).GT.SMALL) THEN
         REND = REND - INTER
      END IF

      RETURN

      END
      SUBROUTINE PLTLAX(X,Y,XLENG,YLENG,TYPE,MINEXP,MAXEXP,LABEL,UNITS)
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
      LOGICAL CPUIFC
      REAL LENMAJ,LENMIN,MAJRAT,MINRAT,NUMRAT,NUMSIZ,LABSIZ,NUMSZM,
     *     LABRAT,MAJTIC,LONNUM,LABSCA,LOGTAB(8),MINEXP,MAXEXP
      CHARACTER*(*) LABEL,UNITS,TYPE
      CHARACTER TTYPE*1,LINE1*10,LINE*10,TLABEL*132,CINPUT*132
      INTEGER PLTITL
      LOGICAL FIRST,PLTNER
      DATA MAJRAT/50./
      DATA MINRAT/2./
      DATA NUMRAT/50./
      DATA LABRAT/40./
      DATA FIRST/.TRUE./
      DATA FUDGE/.0025/

      IF (MAXEXP.LE.MINEXP) THEN
         CALL PLTFLU
         CALL SIORPT('PLTLAX',
     *         'Maximum exponent must be greater than minimum exponent.'
     *               ,2)
         RETURN

      END IF

      IF (FIRST) THEN
         DO 2530 I = 1,8
            LOGTAB(I) = ALOG10(FLOAT(I+1))
 2530    CONTINUE
         FIRST = .FALSE.
      END IF

      LOWEXP = PLTITL(MINEXP)
      NUM = LOWEXP
      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      CALL PLTSVT
      CALL PLTSVD
      CALL PLTSVV
      CALL PLTSTT(3,0.)
      CALL PLTSTT(4,0.)
      CALL PLTSTV(1,1.)
      LINE(1:4) = '10\^'
 2550 CONTINUE
      IF (TTYPE.EQ.'X') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(44)/5.))/NUMRAT
         NUMSZM = NUMSIZ*.8
         LABSIZ = (LABSCA* (GRAPHP(45)/5.))/LABRAT
         LENMAJ = XLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         LFLAG = 0
         MAJTIC = XLENG/ (MAXEXP-MINEXP)
         FMJTIC = X - (MINEXP-LOWEXP)*MAJTIC
         XMAJ = FMJTIC
         YOFF = Y - NUMSIZ*1.8
         YOFFM = YOFF + (NUMSIZ-NUMSZM)/2.
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTVCT(1,X,Y,X+XLENG,Y)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X,Y+YLENG,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,X) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(X-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(X-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

 2580    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2600

         END IF

         DO 2610 I = 1,8
            XNEW = XMAJ + LOGTAB(I)*MAJTIC
            IF (PLTNER(XNEW,X) .AND. NUMSZM.GT.0.0 .AND.
     *          GRAPHP(48).EQ.1.) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
               END IF

            END IF

            IF (XNEW-FUDGE.LE.X) THEN
               GO TO 2610

            END IF

            IF (XNEW.GT.X+XLENG-FUDGE) THEN
               GO TO 2620

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,XNEW,Y,XNEW,Y+LENMIN)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XNEW,Y+YLENG,XNEW,Y+YLENG-LENMIN)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG-LENMIN)

               ELSE
                  CALL PLTVCT(1,XNEW,Y+LENMIN,XNEW,Y+YLENG)
               END IF

               CALL PLTSTV(1,1.)
            END IF

            IF (GRAPHP(48).EQ.1. .AND. NUMSIZ.GT.0.0) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
               END IF

            END IF

 2610    CONTINUE
 2620    CONTINUE
         XMAJ = XMAJ + MAJTIC
         NUM = NUM + 1
         IF (XMAJ.GT.X+XLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2600

         ELSE IF (XMAJ.GT.X+XLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2600

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,XMAJ,Y,XMAJ,Y+LENMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,XMAJ,Y+YLENG,XMAJ,Y+YLENG-LENMAJ)
         END IF

         IF (GRAPHP(35).NE.0. .OR. GRAPHP(73).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG-LENMAJ)

            ELSE
               CALL PLTVCT(1,XMAJ,Y+LENMAJ,XMAJ,Y+YLENG)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

 2590    GO TO 2580

 2600    CONTINUE
         IF (PLTNER(XNEW,X+XLENG) .AND. NUMSZM.GT.0.0 .AND.
     *       GRAPHP(48).EQ.1.) THEN
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTSTT(2,NUMSZM)
               CALL PLTLOA(LINE1,I+1,0)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XNEW-TLEN/2.,YOFFM,LINE1)

            ELSE
               CALL PLTSTT(2,NUMSIZ)
               CALL PLTLOD(LINE1,I+1,NUM-1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XNEW-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

         IF (PLTNER(XMAJ,X+XLENG) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.1. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE)
               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               CALL PLTXTS(XMAJ-TLEN/2.,YOFF,LINE1)
            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2570

         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (LFLAG.EQ.0 .AND. (GRAPHP(40).EQ.1..OR.
     *          GRAPHP(40).EQ.3.) .AND. LOWEXP.NE.0 .AND.
     *          NUMSIZ.GT.0.) THEN
               CALL PLTLOA(LINE1,LOWEXP,0)
               CALL CHRTRM(LINE1,L)
               TLABEL(LL+1:) = ' (*10\^'//LINE1(:L)//'\-)'
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X + (XLENG-TLEN)/2.
               YLAB = YOFF - LABSIZ*2.0
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2570

            END IF

         END IF

      ELSE IF (TTYPE.EQ.'Y') THEN
         LABSCA = (XLENG+YLENG)/2.
         NUMSIZ = (LABSCA* (GRAPHP(88)/5.))/NUMRAT
         NUMSZM = NUMSIZ*.8
         LABSIZ = (LABSCA* (GRAPHP(89)/5.))/LABRAT
         LENMAJ = YLENG/MAJRAT
         LENMIN = LENMAJ/MINRAT
         LFLAG = 0
         MAJTIC = YLENG/ (MAXEXP-MINEXP)
         FMJTIC = Y - (MINEXP-LOWEXP)*MAJTIC
         YMAJ = FMJTIC
         XOFF = NUMSIZ*.8
         YOFF = NUMSIZ/2.
         YOFFM = NUMSZM/2.
         LONNUM = 0.
         CALL PLTSTV(2,GRAPHP(62))
         CALL PLTSTD(1,GRAPHP(37))
         CALL PLTVCT(1,X,Y,X,Y+YLENG)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,Y,X+XLENG,Y+YLENG)
         END IF

         IF (PLTNER(FMJTIC,Y) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),Y-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,Y-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),Y-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

 2630    CONTINUE
         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2650

         END IF

         DO 2660 I = 1,8
            YNEW = YMAJ + LOGTAB(I)*MAJTIC
            IF (PLTNER(YNEW,Y) .AND. NUMSZM.GT.0.0 .AND.
     *          GRAPHP(49).EQ.1.) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFFM,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

                  CALL PLTSTT(2,NUMSIZ)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

               END IF

            END IF

            IF (YNEW-FUDGE.LE.Y) THEN
               GO TO 2660

            END IF

            IF (YNEW.GT.Y+YLENG-FUDGE) THEN
               GO TO 2670

            END IF

            CALL PLTSTD(1,GRAPHP(77))
            CALL PLTSTV(2,GRAPHP(67))
            CALL PLTVCT(1,X,YNEW,X+LENMIN,YNEW)
            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+XLENG,YNEW,X+XLENG-LENMIN,YNEW)
            END IF

            IF (GRAPHP(73).NE.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))
               CALL PLTSTD(1,GRAPHP(74))
               CALL PLTSTV(2,GRAPHP(69))
               IF (GRAPHP(32).EQ.1.) THEN
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG-LENMIN,YNEW)

               ELSE
                  CALL PLTVCT(1,X+LENMIN,YNEW,X+XLENG,YNEW)
               END IF

               CALL PLTSTV(1,1.)
            END IF

            IF (GRAPHP(49).EQ.1. .AND. NUMSIZ.GT.0.0) THEN
               CALL PLTSTT(11,GRAPHP(64))
               CALL PLTSTD(1,GRAPHP(76))
               IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
                  CALL PLTSTT(2,NUMSZM)
                  CALL PLTLOA(LINE1,I+1,0)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

                  CALL PLTSTT(2,NUMSIZ)

               ELSE
                  CALL PLTSTT(2,NUMSIZ)
                  CALL PLTLOD(LINE1,I+1,NUM)
                  CALL PLTXSL(LINE1,TLEN)
                  IF (GRAPHP(92).EQ.1.) THEN
                     LONNUM = NUMSIZ
                     CALL PLTSTT(3,90.)
                     CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                     CALL PLTSTT(3,0.)

                  ELSE
                     CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                     IF (TLEN.GT.LONNUM) THEN
                        LONNUM = TLEN
                     END IF

                  END IF

               END IF

            END IF

 2660    CONTINUE
 2670    CONTINUE
         YMAJ = YMAJ + MAJTIC
         NUM = NUM + 1
         IF (YMAJ.GT.Y+YLENG+FUDGE .AND. GRAPHP(32).EQ.0.) THEN
            GO TO 2650

         END IF

         IF (YMAJ.GT.Y+YLENG-FUDGE .AND. GRAPHP(32).EQ.1.) THEN
            GO TO 2650

         END IF

         CALL PLTSTD(1,GRAPHP(77))
         CALL PLTSTV(2,GRAPHP(67))
         CALL PLTVCT(1,X,YMAJ,X+LENMAJ,YMAJ)
         IF (GRAPHP(32).EQ.1.) THEN
            CALL PLTVCT(1,X+XLENG,YMAJ,X+XLENG-LENMAJ,YMAJ)
         END IF

         IF (GRAPHP(35).NE.0. .OR. GRAPHP(73).NE.0.) THEN
            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(1,GRAPHP(73))

            ELSE
               CALL PLTSTV(1,GRAPHP(35))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTD(1,GRAPHP(74))

            ELSE
               CALL PLTSTD(1,GRAPHP(36))
            END IF

            IF (GRAPHP(35).EQ.0.) THEN
               CALL PLTSTV(2,GRAPHP(69))

            ELSE
               CALL PLTSTV(2,GRAPHP(68))
            END IF

            IF (GRAPHP(32).EQ.1.) THEN
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG-LENMAJ,YMAJ)

            ELSE
               CALL PLTVCT(1,X+LENMAJ,YMAJ,X+XLENG,YMAJ)
            END IF

            CALL PLTSTV(1,1.)
         END IF

         IF (NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

 2640    GO TO 2630

 2650    CONTINUE
         IF (PLTNER(YNEW,Y+YLENG) .AND. NUMSIZ.GT.0.0 .AND.
     *       GRAPHP(49).EQ.1.) THEN
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTSTT(2,NUMSZM)
               CALL PLTLOA(LINE1,I+1,0)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFFM,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            ELSE
               CALL PLTSTT(2,NUMSIZ)
               CALL PLTLOD(LINE1,I+1,NUM-1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YNEW-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YNEW-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

         IF (PLTNER(YMAJ,Y+YLENG) .AND. NUMSIZ.GT.0.0) THEN
            CALL PLTSTT(2,NUMSIZ)
            CALL PLTSTT(11,GRAPHP(64))
            CALL PLTSTD(1,GRAPHP(76))
            IF (GRAPHP(40).EQ.2. .OR. GRAPHP(40).EQ.3.) THEN
               CALL PLTLOA(LINE1,NUM,0)
               LINE(5:) = LINE1
               CALL PLTXSL(LINE,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

               LFLAG = 1

            ELSE
               CALL PLTLOA(LINE1,NUM,1)
               CALL PLTXSL(LINE1,TLEN)
               IF (GRAPHP(92).EQ.1.) THEN
                  LONNUM = NUMSIZ
                  CALL PLTSTT(3,90.)
                  CALL PLTXTS(X-XOFF,YMAJ-TLEN/2.,LINE1)
                  CALL PLTSTT(3,0.)

               ELSE
                  CALL PLTXTS(X- (TLEN+XOFF),YMAJ-YOFF,LINE1)
                  IF (TLEN.GT.LONNUM) THEN
                     LONNUM = TLEN
                  END IF

               END IF

            END IF

         END IF

         IF (CPUIFC(.FALSE.)) THEN
            GO TO 2570

         END IF

         IF (LABSIZ.GT.0.0) THEN
            TLABEL = ' '
            LL = 0
            IF (LABEL.NE.' ') THEN
               TLABEL = LABEL
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (UNITS.NE.' ') THEN
               TLABEL(LL+6:) = UNITS
               CALL CHRTRM(TLABEL,LL)
               TLABEL(LL+1:LL+2) = '\-'
               LL = LL + 2
            END IF

            IF (LFLAG.EQ.0 .AND. (GRAPHP(40).EQ.2..OR.
     *          GRAPHP(40).EQ.3.) .AND. LOWEXP.NE.0 .AND.
     *          NUMSIZ.GT.0.) THEN
               CALL PLTLOA(LINE1,LOWEXP,0)
               CALL CHRTRM(LINE1,L)
               TLABEL(LL+1:) = ' (*10\^'//LINE1(:L)//'\-)'
            END IF

            CINPUT = TLABEL
            CALL CHRSTR(CINPUT,TLABEL,LL)
            IF (LL.GT.0) THEN
               CALL PLTSTT(2,LABSIZ)
               CALL PLTSTT(11,GRAPHP(65))
               CALL PLTSTD(1,GRAPHP(39))
               CALL PLTSTT(3,90.)
               CALL PLTXSL(TLABEL(1:LL),TLEN)
               XLAB = X - LONNUM - LABSIZ*1.4
               YLAB = Y + (YLENG-TLEN)/2.
               CALL PLTXTS(XLAB,YLAB,TLABEL(1:LL))
               CALL PLTSTT(3,0.)
            END IF

            IF (CPUIFC(.FALSE.)) THEN
               GO TO 2570

            END IF

         END IF

      ELSE
         CALL PLTFLU
         CALL SIORPT('PLTLAX','Invalid axis type - '//TTYPE,2)
      END IF

 2560 IF (.NOT. (.TRUE.)) GO TO 2550
 2570 CONTINUE
      CALL PLTRET
      CALL PLTRED
      CALL PLTREV
      RETURN

      END
      LOGICAL FUNCTION PLTNER(X,Y)

      IF (ABS(X-Y).LT..001) THEN
         PLTNER = .TRUE.
         RETURN

      ELSE
         PLTNER = .FALSE.
         RETURN

      END IF

      END
      INTEGER FUNCTION PLTITL(REALN)

      IF (PLTFRC(REALN).EQ.0.) THEN
         PLTITL = INT(REALN)
         RETURN

      END IF

      AREAL = ABS(REALN)
      PLTITL = INT(AREAL)
      IF (REALN.LT.0.) THEN
         PLTITL = -PLTITL - 1
      END IF

      RETURN

      END
      FUNCTION PLTFRC(REALN)

      PLTFRC = REALN - INT(REALN)
      RETURN

      END
      SUBROUTINE PLTLOA(LINE1,NUM,TYPE)
      CHARACTER*10 LINE
      CHARACTER*(*) LINE1
      INTEGER TYPE

      LINE1 = ' '
      LINE = ' '
      IF (TYPE.EQ.1) THEN
         IF (NUM.GE.0) THEN
            WRITE (LINE,10,ERR=20) INT(10.**NUM)

   10       FORMAT (I10)

   20       DO 2680 J = 1,10
               IF (LINE(J:J).NE.' ') THEN
                  GO TO 2690

               END IF

 2680       CONTINUE
 2690       CONTINUE
            LINE1 = LINE(J:)
            RETURN

         END IF

         LINE1(1:1) = '.'
         DO 2700 I = 1,ABS(NUM) - 1
            LINE1(I+1:I+1) = '0'
 2700    CONTINUE
         LINE1(I+1:I+1) = '1'

      ELSE
         WRITE (LINE,10,ERR=40) NUM
   40    DO 2720 J = 1,10
            IF (LINE(J:J).NE.' ') THEN
               GO TO 2730

            END IF

 2720    CONTINUE
 2730    CONTINUE
         LINE1 = LINE(J:)
      END IF

      RETURN

      END
      SUBROUTINE PLTLOD(LINE1,J,NUM)
      CHARACTER*(*) LINE1
      CHARACTER*10 LINE

      LINE1 = ' '
      LINE = ' '
      IF (NUM.GE.0) THEN
         WRITE (LINE,10,ERR=20) INT(FLOAT(J)*10.**NUM)

   10    FORMAT (I10)

   20    DO 2740 I = 1,10
            IF (LINE(I:I).NE.' ') THEN
               GO TO 2750

            END IF

 2740    CONTINUE
 2750    CONTINUE
         LINE1 = LINE(I:)
         RETURN

      END IF

      LINE1(1:1) = '.'
      DO 2760 I = 1,ABS(NUM) - 1
         LINE1(I+1:I+1) = '0'
 2760 CONTINUE
      LINE1(I+1:I+1) = CHAR(J+48)
      RETURN

      END
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
      LOGICAL CPUIFC
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      REAL INTERY,MINEXX,MAXEXX

 2780 CONTINUE
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

         TEMP = ALOG10(XMIN)
         MINEXX = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXX.NE.TEMP) THEN
            MINEXX = MINEXX - 1
         END IF

         TEMP = ALOG10(XMAX)
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

         MINEXX = ALOG10(XMIN)
         MAXEXX = ALOG10(XMAX)
         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         TINT = (GRAPHP(29)-GRAPHP(28))/GRAPHP(30)
         IEXPY = NINT(ALOG10(ABS(TINT)))
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

         MINEXX = ALOG10(XMIN)
         MAXEXX = ALOG10(XMAX)
         TENMNX = 10.**MINEXX
         TENMXX = 10.**MAXEXX
         IEXPY = NINT(ALOG10(ABS(GRAPHP(86))))
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
         IEXPY = GRAPHP(90)
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
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2800

      END IF

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'y',YSTART,
     *            YEND,FNLOWY,INT(GRAPHP(42)),INTERY,NMINY,YLAB,YUNIT,
     *            IEXPY)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2800

      END IF

      CALL PLTGM2(MINEXX,MAXEXX,YSTART*TNEXPY,YEND*TNEXPY,GRAPHP(1),
     *            GRAPHP(1)+GRAPHP(3),GRAPHP(2),GRAPHP(2)+GRAPHP(4),
     *            GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2800

      END IF

 2790 IF (.NOT. (.TRUE.)) GO TO 2780
 2800 CONTINUE
      RETURN

      END
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
      LOGICAL CPUIFC
      REAL MINEXX,MAXEXX,MINEXY,MAXEXY
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      DIMENSION X(1),Y(1)

 2810 CONTINUE
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

         TEMP = ALOG10(XMIN)
         MINEXX = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXX.NE.TEMP) THEN
            MINEXX = MINEXX - 1
         END IF

         TEMP = ALOG10(XMAX)
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

         TEMP = ALOG10(YMIN)
         MINEXY = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXY.NE.TEMP) THEN
            MINEXY = MINEXY - 1
         END IF

         TEMP = ALOG10(YMAX)
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

         MINEXX = ALOG10(XMIN)
         MAXEXX = ALOG10(XMAX)
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

         MINEXY = ALOG10(YMIN)
         MAXEXY = ALOG10(YMAX)
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

         MINEXX = ALOG10(XMIN)
         MAXEXX = ALOG10(XMAX)
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

         MINEXY = ALOG10(YMIN)
         MAXEXY = ALOG10(YMAX)
         TENMNY = 10.**MINEXY
         TENMXY = 10.**MAXEXY
      END IF

      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2830

      END IF

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'x',MINEXX,MAXEXX,
     *            XLAB,XUNIT)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2830

      END IF

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'y',MINEXY,MAXEXY,
     *            YLAB,YUNIT)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2830

      END IF

      CALL PLTGM2(MINEXX,MAXEXX,MINEXY,MAXEXY,GRAPHP(1),GRAPHP(1)+XLENT,
     *            GRAPHP(2),GRAPHP(2)+YLENT,GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2830

      END IF

 2820 IF (.NOT. (.TRUE.)) GO TO 2810
 2830 CONTINUE
      RETURN

      END
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
      LOGICAL CPUIFC
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      REAL INTERX,MINEXY,MAXEXY
      DIMENSION X(1),Y(1)

 2840 CONTINUE
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

         TEMP = ALOG10(YMIN)
         MINEXY = INT(TEMP)
         IF (TEMP.LT.0. .AND. MINEXY.NE.TEMP) THEN
            MINEXY = MINEXY - 1
         END IF

         TEMP = ALOG10(YMAX)
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
         IEXPX = NINT(ALOG10(ABS(TINT)))
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

         MINEXY = ALOG10(YMIN)
         MAXEXY = ALOG10(YMAX)
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
         IEXPX = NINT(ALOG10(ABS(GRAPHP(81))))
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

         MINEXY = ALOG10(YMIN)
         MAXEXY = ALOG10(YMAX)
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
         IEXPX = GRAPHP(91)
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

      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2860

      END IF

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'x',XSTART,
     *            XEND,FNLOWX,INT(GRAPHP(41)),INTERX,NMINX,XLAB,XUNIT,
     *            IEXPX)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2860

      END IF

      CALL PLTLAX(GRAPHP(1),GRAPHP(2),GRAPHP(3),GRAPHP(4),'y',MINEXY,
     *            MAXEXY,YLAB,YUNIT)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2860

      END IF

      CALL PLTGM2(XSTART*TNEXPX,XEND*TNEXPX,MINEXY,MAXEXY,GRAPHP(1),
     *            GRAPHP(1)+GRAPHP(3),GRAPHP(2),GRAPHP(2)+GRAPHP(4),
     *            GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTCUR(X,Y,NUM)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2860

      END IF

 2850 IF (.NOT. (.TRUE.)) GO TO 2840
 2860 CONTINUE
      RETURN

      END
      SUBROUTINE PLTNCF(X,TYPE,FN,NE)
      CHARACTER*(*) TYPE
      CHARACTER*1 TTYPE
      REAL FNICE(17)
      DATA FNICE/-10.,-8.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,
     *     8.,10./

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      IF (X.EQ.0.) THEN
         FN = 0.
         RETURN

      END IF

      F1 = X/10.**NE
      IF (TTYPE.EQ.'O') THEN
         DO 2870 I = 1,17
            IF (F1.LE.FNICE(I)) THEN
               GO TO 2880

            END IF

 2870    CONTINUE
 2880    CONTINUE
      END IF

      IF (TTYPE.EQ.'U') THEN
         DO 2890 I = 17,1,-1
            IF (F1.GE.FNICE(I)) THEN
               GO TO 2900

            END IF

 2890    CONTINUE
 2900    CONTINUE
      END IF

      IF (I.GT.17) THEN
         I = 17
      END IF

      IF (I.LT.1) THEN
         I = 1
      END IF

      FN = FNICE(I)
      RETURN

      END
      SUBROUTINE PLTNIC(X,TYPE,FN,NE,INTER,NMIN)
      CHARACTER*(*) TYPE
      CHARACTER*1 TTYPE
      REAL FNICE(11),INTERA(11),INTER
      INTEGER NMINA(11)
      DATA FNICE/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
      DATA INTERA/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./
      DATA NMINA/10,10,10,10,5,5,5,5,5,5,5/

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      IF (X.EQ.0.) THEN
         FN = 0.
         NE = 0
         RETURN

      END IF

      XTEMP = ABS(X)
      IF (X.LT.0. .AND. TTYPE.EQ.'O') THEN
         TTYPE = 'U'

      ELSE IF (X.LT.0. .AND. TTYPE.EQ.'U') THEN
         TTYPE = 'O'
      END IF

      E1 = ALOG10(XTEMP)
      JE = INT(E1)
      IF (JE.LE.0. .AND. XTEMP.LT.1.) THEN
         JE = JE - 1
      END IF

      IF (JE.LT.0) THEN
         E2 = 1./10.** (-JE)

      ELSE IF (JE.EQ.0) THEN
         E2 = 1.

      ELSE IF (JE.GT.0) THEN
         E2 = 10.**JE
      END IF

      F1 = XTEMP/E2
      IF (TTYPE.EQ.'O') THEN
         DO 2910 I = 1,11
            IF (F1/1.007.LE.FNICE(I)) THEN
               GO TO 2920

            END IF

 2910    CONTINUE
 2920    CONTINUE
      END IF

      IF (TTYPE.EQ.'U') THEN
         DO 2930 I = 11,1,-1
            IF (F1*1.007.GE.FNICE(I)) THEN
               GO TO 2940

            END IF

 2930    CONTINUE
 2940    CONTINUE
      END IF

      IF (I.GT.11) THEN
         I = 11
      END IF

      IF (I.LT.1) THEN
         I = 1
      END IF

      FN = FNICE(I)
      INTER = INTERA(I)
      NMIN = NMINA(I)
      NE = JE
      IF (X.LT.0.) THEN
         FN = -FN
      END IF

      IF (FN.EQ.0.) THEN
         NE = 0
      END IF

      RETURN

      END
      SUBROUTINE PLTNXY(X,Y,NUM,XLAB,XUNIT,YLAB,YUNIT)
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
      LOGICAL CPUIFC
      CHARACTER*(*) XLAB,XUNIT,YLAB,YUNIT
      DIMENSION X(1),Y(1)
      REAL INTERX,INTERY

 2950 CONTINUE
      XLENT = GRAPHP(3)
      YLENT = GRAPHP(4)
      CALL VECRGS(IABS(NUM),X,XMAX,XMIN)
      CALL VECRGS(IABS(NUM),Y,YMAX,YMIN)
      IF (GRAPHP(22).EQ.1.) THEN
         CALL PLTINO(XMIN,XMAX,FNLOWX,FNUPPX,INTERX,IEXPX,NMINX)
         XSTART = FNLOWX
         XEND = FNUPPX
         TNEXPX = 10.**IEXPX
         CALL PLTINO(YMIN,YMAX,FNLOWY,FNUPPY,INTERY,IEXPY,NMINY)
         YSTART = FNLOWY
         YEND = FNUPPY
         TNEXPY = 10.**IEXPY
         GRAPHP(24) = FNLOWX*TNEXPX
         GRAPHP(25) = FNUPPX*TNEXPX
         GRAPHP(26) = (FNUPPX-FNLOWX)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = FNLOWY*TNEXPY
         GRAPHP(29) = FNUPPY*TNEXPY
         GRAPHP(30) = (FNUPPY-FNLOWY)/INTERY
         GRAPHP(31) = NMINY
         GRAPHP(78) = FNLOWX*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = FNUPPX*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = FNLOWY*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = FNUPPY*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.2.) THEN
         DELX = XMAX - XMIN
         DELY = YMAX - YMIN
         CALL PLTINO(XMIN,XMAX,FNLOWX,FNUPPX,INTERX,IEXPX,NMINX)
         CALL PLTINO(YMIN,YMAX,FNLOWY,FNUPPY,INTERY,IEXPY,NMINY)
         IF (DELX.GT.DELY) THEN
            IF (IEXPX.NE.IEXPY) THEN
               CALL PLTNCF(YMIN,'u',FNLOWY,IEXPX)
               IEXPY = IEXPX
               INTERY = INTERX
               IF (INTERY.EQ.2. .AND. AMOD(ABS(FNLOWY),2.).EQ.1.) THEN
                  FNLOWY = FNLOWY - 1.
               END IF

               NMINY = NMINX
            END IF

            FNUPPY = FNLOWY + FNUPPX - FNLOWX
            IF (FNUPPY*10.**IEXPY.LT.YMAX) THEN
               FNUPPY = FNUPPY + 1.
               FNUPPX = FNUPPX + 1.
            END IF

         END IF

         IF (DELX.LT.DELY) THEN
            IF (IEXPX.NE.IEXPY) THEN
               CALL PLTNCF(XMIN,'u',FNLOWX,IEXPY)
               IEXPX = IEXPY
               INTERX = INTERY
               IF (INTERX.EQ.2. .AND. AMOD(ABS(FNLOWX),2.).EQ.1.) THEN
                  FNLOWX = FNLOWX - 1.
               END IF

               NMINX = NMINY
            END IF

            FNUPPX = FNLOWX + FNUPPY - FNLOWY
            IF (FNUPPX*10.**IEXPX.LT.XMAX) THEN
               FNUPPX = FNUPPX + 1.
               FNUPPY = FNUPPY + 1.
            END IF

         END IF

         IF (GRAPHP(3).NE.GRAPHP(4)) THEN
            XLENT = AMIN1(GRAPHP(3),GRAPHP(4))
            YLENT = XLENT
         END IF

         TNEXPX = 10.**IEXPX
         TNEXPY = 10.**IEXPY
         XSTART = FNLOWX
         XEND = FNUPPX
         YSTART = FNLOWY
         YEND = FNUPPY
         GRAPHP(24) = FNLOWX*TNEXPX
         GRAPHP(25) = FNUPPX*TNEXPX
         GRAPHP(26) = (FNUPPX-FNLOWX)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = FNLOWY*TNEXPY
         GRAPHP(29) = FNUPPY*TNEXPY
         GRAPHP(30) = (FNUPPY-FNLOWY)/INTERY
         GRAPHP(31) = NMINY
         GRAPHP(78) = FNLOWX*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = FNUPPX*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = FNLOWY*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = FNUPPY*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.3.) THEN
         TINT = (GRAPHP(25)-GRAPHP(24))/GRAPHP(26)
         IEXPX = NINT(ALOG10(ABS(TINT)))
         TNEXPX = 10.**IEXPX
         FNLOWX = GRAPHP(24)/TNEXPX
         FNUPPX = GRAPHP(25)/TNEXPX
         INTERX = (FNUPPX-FNLOWX)/NINT(GRAPHP(26))
         NMINX = INT(GRAPHP(27))
         XSTART = FNLOWX
         XEND = FNUPPX
         TINT = (GRAPHP(29)-GRAPHP(28))/GRAPHP(30)
         IEXPY = NINT(ALOG10(ABS(TINT)))
         TNEXPY = 10.**IEXPY
         FNLOWY = GRAPHP(28)/TNEXPY
         FNUPPY = GRAPHP(29)/TNEXPY
         INTERY = (FNUPPY-FNLOWY)/NINT(GRAPHP(30))
         NMINY = INT(GRAPHP(31))
         YSTART = FNLOWY
         YEND = FNUPPY
         GRAPHP(78) = XSTART*TNEXPX
         GRAPHP(79) = GRAPHP(78)
         GRAPHP(80) = XEND*TNEXPX
         GRAPHP(81) = INTERX*TNEXPX
         GRAPHP(82) = NMINX
         GRAPHP(83) = YSTART*TNEXPY
         GRAPHP(84) = GRAPHP(83)
         GRAPHP(85) = YEND*TNEXPY
         GRAPHP(86) = INTERY*TNEXPY
         GRAPHP(87) = NMINY

      ELSE IF (GRAPHP(22).EQ.4.) THEN
         IEXPX = NINT(ALOG10(ABS(GRAPHP(81))))
         TNEXPX = 10.**IEXPX
         XSTART = GRAPHP(78)/TNEXPX
         XEND = GRAPHP(80)/TNEXPX
         FNLOWX = GRAPHP(79)/TNEXPX
         INTERX = GRAPHP(81)/TNEXPX
         NMINX = INT(GRAPHP(82))
         IEXPY = NINT(ALOG10(ABS(GRAPHP(86))))
         TNEXPY = 10.**IEXPY
         YSTART = GRAPHP(83)/TNEXPY
         YEND = GRAPHP(85)/TNEXPY
         FNLOWY = GRAPHP(84)/TNEXPY
         INTERY = GRAPHP(86)/TNEXPY
         NMINY = INT(GRAPHP(87))
         GRAPHP(24) = XSTART*TNEXPX
         GRAPHP(25) = XEND*TNEXPX
         GRAPHP(26) = (XSTART-XEND)/INTERX
         GRAPHP(27) = NMINX
         GRAPHP(28) = YSTART*TNEXPY
         GRAPHP(29) = YEND*TNEXPY
         GRAPHP(30) = (YSTART-YEND)/INTERY
         GRAPHP(31) = NMINY
      END IF

      IF (GRAPHP(91).NE.-999999.) THEN
         FAC = 10.** (IEXPX-GRAPHP(91))
         IEXPX = GRAPHP(91)
         TNEXPX = 10.**IEXPX
         XSTART = XSTART*FAC
         XEND = XEND*FAC
         FNLOWX = FNLOWX*FAC
         INTERX = INTERX*FAC
      END IF

      IF (GRAPHP(90).NE.-999999.) THEN
         FAC = 10.** (IEXPY-GRAPHP(90))
         IEXPY = GRAPHP(90)
         TNEXPY = 10.**IEXPY
         YSTART = YSTART*FAC
         YEND = YEND*FAC
         FNLOWY = FNLOWY*FAC
         INTERY = INTERY*FAC
      END IF

      IF (GRAPHP(40).EQ.1.) THEN
         YSTART = YSTART*TNEXPY
         YEND = YEND*TNEXPY
         FNLOWY = FNLOWY*TNEXPY
         FNUPPY = FNUPPY*TNEXPY
         INTERY = INTERY*TNEXPY
         IEXPY = 0
         TNEXPY = 1.
      END IF

      IF (GRAPHP(40).EQ.2.) THEN
         XSTART = XSTART*TNEXPX
         XEND = XEND*TNEXPX
         FNLOWX = FNLOWX*TNEXPX
         FNUPPX = FNUPPX*TNEXPX
         INTERX = INTERX*TNEXPX
         IEXPX = 0
         TNEXPX = 1.
      END IF

      IF (GRAPHP(40).EQ.4.) THEN
         XSTART = XSTART*TNEXPX
         XEND = XEND*TNEXPX
         FNLOWX = FNLOWX*TNEXPX
         FNUPPX = FNUPPX*TNEXPX
         INTERX = INTERX*TNEXPX
         IEXPX = 0
         TNEXPX = 1.
         YSTART = YSTART*TNEXPY
         YEND = YEND*TNEXPY
         FNLOWY = FNLOWY*TNEXPY
         FNUPPY = FNUPPY*TNEXPY
         INTERY = INTERY*TNEXPY
         IEXPY = 0
         TNEXPY = 1.
      END IF

      CALL PLTGM2(XSTART*TNEXPX,XEND*TNEXPX,YSTART*TNEXPY,
     *            YEND*TNEXPY,GRAPHP(1),GRAPHP(1)+XLENT,GRAPHP(2),
     *            GRAPHP(2)+YLENT,GRAPHP(7))
      CALL PLTUWN(GRAPHP(7))
      CALL PLTAXS(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'x',XSTART,XEND,
     *            FNLOWX,INT(GRAPHP(41)),INTERX,NMINX,XLAB,XUNIT,IEXPX)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2970

      END IF

      CALL PLTAXS(GRAPHP(1),GRAPHP(2),XLENT,YLENT,'y',YSTART,YEND,
     *            FNLOWY,INT(GRAPHP(42)),INTERY,NMINY,YLAB,YUNIT,IEXPY)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2970

      END IF

      CALL PLTCUR(X,Y,NUM)
      IF (CPUIFC(.FALSE.)) THEN
         GO TO 2970

      END IF

 2960 IF (.NOT. (.TRUE.)) GO TO 2950
 2970 CONTINUE
      RETURN

      END
      SUBROUTINE PLTUWN(GMAP)
      DIMENSION GMAP(*)

      GMAP(7) = GMAP(7) - .0002
      GMAP(8) = GMAP(8) - .0002
      GMAP(9) = GMAP(9) + .0006
      GMAP(10) = GMAP(10) - .0002
      GMAP(11) = GMAP(11) + .0006
      GMAP(12) = GMAP(12) + .0006
      GMAP(13) = GMAP(13) - .0002
      GMAP(14) = GMAP(14) + .0006
      RETURN

      END
      SUBROUTINE PLTAV2(UMAP,N,X1,Y1,X2,Y2,TH,XL)
      REAL UMAP(*)
      INTEGER N
      REAL X1(*),Y1(*)
      REAL X2(*),Y2(*)
      REAL TH
      REAL XL
      REAL PX(32),PY(32),QX(32),QY(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      MASK = -1
      CALL PLTMV2(UMAP,N,MASK,X1,Y1,X2,Y2,PX,PY,QX,QY)
      DO 2000 I = 1,N
         IF (IAND(MASK,IZBIT(I)).NE.0) THEN
            CALL PLTARR(PX(I),PY(I),QX(I),QY(I),TH,XL)
         END IF

 2000 CONTINUE
      RETURN

      END
      SUBROUTINE PLTAV3(UMAP,N,X1,Y1,Z1,X2,Y2,Z2,TH,XL)
      REAL UMAP(*)
      INTEGER N
      REAL X1(*),Y1(*),Z1(*)
      REAL X2(*),Y2(*),Z2(*)
      REAL TH
      REAL XL
      REAL PX(32),PY(32),QX(32),QY(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      MASK = -1
      CALL PLTMV3(UMAP,N,MASK,X1,Y1,Z1,X2,Y2,Z2,PX,PY,QX,QY)
      DO 2020 I = 1,N
         IF (IAND(MASK,IZBIT(I)).NE.0) THEN
            CALL PLTARR(PX(I),PY(I),QX(I),QY(I),TH,XL)
         END IF

 2020 CONTINUE
      RETURN

      END
      SUBROUTINE PLTGM2(XL,XU,YL,YU,PXL,PXU,PYL,PYU,UMAP)
      REAL UMAP(*)

      DU = XU - XL
      DP = PXU - PXL
      IF (DU.EQ.0.) THEN
         DU = 1.
      END IF

      UMAP(1) = DP/DU
      UMAP(1+1) = 0.
      UMAP(1+2) = 0.
      DV = YU - YL
      DQ = PYU - PYL
      IF (DV.EQ.0.) THEN
         DV = 1.
      END IF

      UMAP(1+3) = DQ/DV
      UMAP(5) = PXL - XL*UMAP(1)
      UMAP(5+1) = PYL - YL*UMAP(1+3)
      UMAP(7) = PXL
      UMAP(7+1) = PYL
      UMAP(9) = PXU
      UMAP(9+1) = PYL
      UMAP(11) = PXU
      UMAP(11+1) = PYU
      UMAP(13) = PXL
      UMAP(13+1) = PYU
      RETURN

      END
      SUBROUTINE PLTGM3(PX,PY,PZ,S,UMAP)
      REAL UMAP(*)
      COMMON /CENBOD/XC,YC,ZC

      UMAP(17) = 1.
      UMAP(18) = PX
      UMAP(18+1) = PY
      UMAP(18+2) = PZ + S*2.5
      UMAP(21) = 1.
      UMAP(21+1) = 0.
      UMAP(21+2) = 0.
      UMAP(24) = 0.
      UMAP(24+1) = 1.
      UMAP(24+2) = 0.
      UMAP(27) = 0.
      UMAP(27+1) = 0.
      UMAP(27+2) = -1.
      UMAP(15) = S*.01
      UMAP(16) = S*100.
      UMAP(30) = 1.
      XC = PX
      YC = PY
      ZC = PZ
      RETURN

      END
      SUBROUTINE PLTMP2(UMP,N,MASK,PX,PY,QX,QY)
      DIMENSION UMP(*),MASK(*),PX(*),PY(*),QX(*),QY(*)

      AXX = UMP(1)
      AYY = UMP(4)
      AXY = UMP(3)
      AYX = UMP(2)
      BX = UMP(5)
      BY = UMP(6)
      DO 2040 I = 1,N
         QX(I) = AXX*PX(I) + AXY*PY(I) + BX
         QY(I) = AYX*PX(I) + AYY*PY(I) + BY
 2040 CONTINUE
      CALL PLTCP2(N,MASK,QX,QY,UMP(7),UMP(9))
      CALL PLTCP2(N,MASK,QX,QY,UMP(9),UMP(11))
      CALL PLTCP2(N,MASK,QX,QY,UMP(11),UMP(13))
      CALL PLTCP2(N,MASK,QX,QY,UMP(13),UMP(7))
      RETURN

      END
      SUBROUTINE PLTMP3(UMAP,N,MASK,PX,PY,PZ,QX,QY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),PZ(*),QX(*),QY(*)
      DIMENSION Q1(3),V1(3),Q2(3),V2(3),TPX(32),TPY(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      DO 2060 L = 1,3
         V1(L) = UMAP(18+L-1) + UMAP(15)*UMAP(27+L-1)
         Q1(L) = UMAP(L+27-1)
         V2(L) = UMAP(18+L-1) + UMAP(16)*UMAP(27+L-1)
         Q2(L) = -UMAP(L+27-1)
 2060 CONTINUE
      J = 0
      KM = 0
 2080 IF (.NOT. (J.LT.N)) GO TO 2090
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V1,Q1)
      CALL PLTCP3(JN,MASK(KM),PX(J1),PY(J1),PZ(J1),V2,Q2)
      IF (UMAP(17).EQ.1.) THEN
         M = MASK(KM)
         DO 2100 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,M).NE.0) THEN
               PMS = (PX(K+J1-1)-UMAP(18))*UMAP(27) +
     *               (PY(K+J1-1)-UMAP(19))*UMAP(28) +
     *               (PZ(K+J1-1)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TPX(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(21)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(22)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(23))
               TPY(K) = R* ((PX(K+J1-1)-UMAP(18))*UMAP(24)+
     *                  (PY(K+J1-1)-UMAP(19))*UMAP(25)+
     *                  (PZ(K+J1-1)-UMAP(20))*UMAP(26))
            END IF

 2100    CONTINUE
      END IF

      IF (UMAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMP2(UMAP(1),JN,MASK(KM),TPX,TPY,QX(J1),QY(J1))
      GO TO 2080

 2090 CONTINUE
      RETURN

      END
      SUBROUTINE PLTMV2(UMAP,N,MASK,PX,PY,QX,QY,PPX,PPY,QQX,QQY)
      DIMENSION UMAP(*),MASK(*),PX(*),PY(*),QX(*),QY(*),PPX(*),PPY(*),
     *          QQX(*),QQY(*)

      AXX = UMAP(1)
      AYY = UMAP(4)
      AXY = UMAP(3)
      AYX = UMAP(2)
      BX = UMAP(5)
      BY = UMAP(6)
      DO 2120 I = 1,N
         PPX(I) = AXX*PX(I) + AXY*PY(I) + BX
         QQX(I) = AXX*QX(I) + AXY*QY(I) + BX
         PPY(I) = AYX*PX(I) + AYY*PY(I) + BY
         QQY(I) = AYX*QX(I) + AYY*QY(I) + BY
 2120 CONTINUE
      J = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J + 1
      KM = 1 + J/32
      J = J + JN
      CALL PLTVWV(UMAP(7),UMAP(11),JN,MASK(KM),PPX(J1),PPY(J1),QQX(J1),
     *            QQY(J1))
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
      SUBROUTINE PLTMV3(UMAP,N,MASK,UX,UY,UZ,VX,VY,VZ,PX,PY,QX,QY)
      DIMENSION UMAP(*),MASK(*),UX(*),UY(*),UZ(*),VX(*),VY(*),VZ(*),
     *          PX(*),PY(*),QX(*),QY(*)
      DIMENSION TUX(32),TUY(32),TUZ(32),TVX(32),TVY(32),TVZ(32),
     *          TTUX(32),TTUY(32),TTUZ(32),TTVX(32),TTVY(32),TTVZ(32),
     *          V1(3),Q1(3),V2(3),Q2(3)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      DO 2160 L = 1,3
         V1(L) = UMAP(18+L-1) + UMAP(15)*UMAP(27+L-1)
         Q1(L) = UMAP(27+L-1)
         V2(L) = UMAP(18+L-1) + UMAP(16)*UMAP(27+L-1)
         Q2(L) = -UMAP(27+L-1)
 2160 CONTINUE
      J = 0
      KM = 0
 2180 IF (.NOT. (J.LT.N)) GO TO 2190
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      CALL PLTCV3(JN,MASK(KM),UX(J1),UY(J1),UZ(J1),VX(J1),VY(J1),VZ(J1),
     *            TUX,TUY,TUZ,TVX,TVY,TVZ,V1,Q1)
      CALL PLTCV3(JN,MASK(KM),TUX,TUY,TUZ,TVX,TVY,TVZ,TTUX,TTUY,TTUZ,
     *            TTVX,TTVY,TTVZ,V2,Q2)
      IF (UMAP(17).EQ.1.) THEN
         DO 2200 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK(KM)).NE.0) THEN
               PMS = (TTUX(K)-UMAP(18))*UMAP(27) +
     *               (TTUY(K)-UMAP(19))*UMAP(28) +
     *               (TTUZ(K)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TUX(K) = R* ((TTUX(K)-UMAP(18))*UMAP(21)+
     *                  (TTUY(K)-UMAP(19))*UMAP(22)+
     *                  (TTUZ(K)-UMAP(20))*UMAP(23))
               TUY(K) = R* ((TTUX(K)-UMAP(18))*UMAP(24)+
     *                  (TTUY(K)-UMAP(19))*UMAP(25)+
     *                  (TTUZ(K)-UMAP(20))*UMAP(26))
               PMS = (TTVX(K)-UMAP(18))*UMAP(27) +
     *               (TTVY(K)-UMAP(19))*UMAP(28) +
     *               (TTVZ(K)-UMAP(20))*UMAP(29)
               R = UMAP(30)/PMS
               TVX(K) = R* ((TTVX(K)-UMAP(18))*UMAP(21)+
     *                  (TTVY(K)-UMAP(19))*UMAP(22)+
     *                  (TTVZ(K)-UMAP(20))*UMAP(23))
               TVY(K) = R* ((TTVX(K)-UMAP(18))*UMAP(24)+
     *                  (TTVY(K)-UMAP(19))*UMAP(25)+
     *                  (TTVZ(K)-UMAP(20))*UMAP(26))
            END IF

 2200    CONTINUE

      ELSE IF (UMAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMV2(UMAP,JN,MASK(KM),TUX,TUY,TVX,TVY,PX(J1),PY(J1),QX(J1),
     *            QY(J1))
      GO TO 2180

 2190 CONTINUE
      RETURN

      END
      SUBROUTINE PLTMG2(MAP,N,XV,YV,NO,XVO,YVO)
      REAL MAP(*)
      INTEGER N
      REAL XV(*),YV(*)
      INTEGER NO
      REAL XVO(*),YVO(*)
      REAL XWORK(50),YWORK(50)
      INTEGER NWORK

      NOSAVE = NO
      AXX = MAP(1)
      AYY = MAP(4)
      AXY = MAP(3)
      AYX = MAP(2)
      BX = MAP(5)
      BY = MAP(6)
      DO 2220 I = 1,N
         XVO(I) = AXX*XV(I) + AXY*YV(I) + BX
         YVO(I) = AYX*XV(I) + AYY*YV(I) + BY
 2220 CONTINUE
      NWORK = 50
      CALL PLTCG2(N,XVO,YVO,NWORK,XWORK,YWORK,MAP(7),MAP(9))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(9),MAP(11))
      NWORK = 50
      CALL PLTCG2(NO,XVO,YVO,NWORK,XWORK,YWORK,MAP(11),MAP(13))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(13),MAP(7))
      RETURN

      END
      SUBROUTINE PLTMZM(FACT,UMAP)
      REAL FACT
      REAL UMAP(*)

      UMAP(30) = FACT
      RETURN

      END
      SUBROUTINE PLTMMV(DX,DY,DZ,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL DX,DY,DZ
      REAL UMAP(*)

      UMAP(18) = UMAP(18) + DX
      UMAP(18+1) = UMAP(18+1) + DY
      UMAP(18+2) = UMAP(18+2) + DZ
      XC = XC + DX
      YC = YC + DY
      ZC = ZC + DZ
      RETURN

      END
      SUBROUTINE PLTMAA(ALT,AZI,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL ALT
      REAL AZI
      REAL UMAP(*)
      REAL R

      R = (UMAP(18)-XC)**2 + (UMAP(18+1)-YC)**2 + (UMAP(18+2)-ZC)**2
      R = SQRT(R)
      UMAP(21) = 1.
      UMAP(21+1) = 0.
      UMAP(21+2) = 0.
      UMAP(24) = 0.
      UMAP(24+1) = 1.
      UMAP(24+2) = 0.
      UMAP(27) = 0.
      UMAP(27+1) = 0.
      UMAP(27+2) = -1.
      CALL PLTRTX(- (90-ALT)*3.1415927/180.,UMAP)
      CALL PLTRTZ(AZI*3.1415927/180.,UMAP)
      UMAP(18) = XC - R*UMAP(27)
      UMAP(18+1) = YC - R*UMAP(27+1)
      UMAP(18+2) = ZC - R*UMAP(27+2)
      RETURN

      END
      SUBROUTINE PLTRTZ(VAL,UMAP)
      REAL UMAP(*)
      REAL VAL
      REAL A(9),B(9)

      S = SIN(VAL)
      C = COS(VAL)
      DO 2240 I = 1,9
         A(I) = UMAP(21-1+I)
         B(I) = 0.
 2240 CONTINUE
      B(1) = C
      B(2) = S
      B(4) = -S
      B(5) = C
      B(9) = 1.
      CALL PLTROT(B,A,UMAP(21))
      RETURN

      END
      SUBROUTINE PLTRTY(VAL,UMAP)
      REAL UMAP(*)
      REAL VAL
      REAL A(9),B(9)

      S = SIN(VAL)
      C = COS(VAL)
      DO 2260 I = 1,9
         A(I) = UMAP(21-1+I)
         B(I) = 0.
 2260 CONTINUE
      B(1) = C
      B(3) = S
      B(7) = -S
      B(9) = C
      B(5) = 1.
      CALL PLTROT(B,A,UMAP(21))
      RETURN

      END
      SUBROUTINE PLTRTX(VAL,UMAP)
      REAL UMAP(*)
      REAL VAL
      REAL A(9),B(9)

      S = SIN(VAL)
      C = COS(VAL)
      DO 2280 I = 1,9
         A(I) = UMAP(21-1+I)
         B(I) = 0.
 2280 CONTINUE
      B(1) = 1.
      B(5) = C
      B(9) = C
      B(6) = S
      B(8) = -S
      CALL PLTROT(B,A,UMAP(21))
      RETURN

      END
      SUBROUTINE PLTROT(R,A,B)
      REAL R(3,3),A(3,3),B(3,3)

      DO 2300 K = 1,3
         DO 2320 I = 1,3
            S = 0.
            DO 2340 J = 1,3
               S = S + R(J,I)*A(J,K)
 2340       CONTINUE
            B(I,K) = S
 2320    CONTINUE
 2300 CONTINUE
      RETURN

      END
      SUBROUTINE PLTMOR(X,Y,Z,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL X,Y,Z
      REAL UMAP(*)

      DX = X - XC
      DY = Y - YC
      DZ = Z - ZC
      CALL PLTMMV(DX,DY,DZ,UMAP)
      RETURN

      END
      SUBROUTINE PLTMMO(FACT,UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      R = (UMAP(18)-XC)**2 + (UMAP(18+1)-YC)**2 + (UMAP(18+2)-ZC)**2
      R = SQRT(R)*FACT
      UMAP(18) = XC - R*UMAP(27)
      UMAP(18+1) = YC - R*UMAP(27+1)
      UMAP(18+2) = ZC - R*UMAP(27+2)
      RETURN

      END
      SUBROUTINE PLTDP2(MAP,N,PX,PY)
      REAL MAP(*),PX(*),PY(*)
      REAL XWORK(32),YWORK(32)

      J = 0
 2360 IF (.NOT. (J.LT.N)) GO TO 2370
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      I = 1
 2380 IF (.NOT. (I.LE.JN)) GO TO 2400
      XWORK(I) = MAP(1)*PX(J1+I) + MAP(3)*PY(J1+I) + MAP(5)
      YWORK(I) = MAP(2)*PX(J1+I) + MAP(4)*PY(J1+I) + MAP(6)
 2390 I = I + 1
      GO TO 2380

 2400 CONTINUE
      MASK = -1
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(7),MAP(9))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(9),MAP(11))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(11),MAP(13))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(13),MAP(7))
      CALL PLTPTM(JN,MASK,XWORK,YWORK)
      GO TO 2360

 2370 CONTINUE
      RETURN

      END
      SUBROUTINE PLTDP3(MAP,N,PX,PY,PZ)
      REAL MAP(*),PX(*),PY(*),PZ(*)
      DIMENSION Q1(3),V1(3),Q2(3),V2(3),TPX(32),TPY(32),QX(32),QY(32)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      DO 2410 L = 1,3
         V1(L) = MAP(18+L-1) + MAP(15)*MAP(27+L-1)
         Q1(L) = MAP(L+27-1)
         V2(L) = MAP(18+L-1) + MAP(16)*MAP(27+L-1)
         Q2(L) = -MAP(L+27-1)
 2410 CONTINUE
      J = 0
      KM = 0
 2430 IF (.NOT. (J.LT.N)) GO TO 2440
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      MASK = -1
      CALL PLTCP3(JN,MASK,PX(J1),PY(J1),PZ(J1),V1,Q1)
      CALL PLTCP3(JN,MASK,PX(J1),PY(J1),PZ(J1),V2,Q2)
      IF (MAP(17).EQ.1.) THEN
         DO 2450 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK).NE.0) THEN
               PMS = (PX(K+J1-1)-MAP(18))*MAP(27) +
     *               (PY(K+J1-1)-MAP(19))*MAP(28) +
     *               (PZ(K+J1-1)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TPX(K) = R* ((PX(K+J1-1)-MAP(18))*MAP(21)+
     *                  (PY(K+J1-1)-MAP(19))*MAP(22)+
     *                  (PZ(K+J1-1)-MAP(20))*MAP(23))
               TPY(K) = R* ((PX(K+J1-1)-MAP(18))*MAP(24)+
     *                  (PY(K+J1-1)-MAP(19))*MAP(25)+
     *                  (PZ(K+J1-1)-MAP(20))*MAP(26))
            END IF

 2450    CONTINUE
      END IF

      IF (MAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMP2(MAP,JN,MASK,TPX,TPY,QX,QY)
      CALL PLTPTM(JN,MASK,QX,QY)
      GO TO 2430

 2440 CONTINUE
      RETURN

      END
      SUBROUTINE PLTDV2(MAP,N,PX,PY,QX,QY)
      REAL MAP(*),PX(*),PY(*),QX(*),QY(*)
      DIMENSION PPX(32),PPY(32),QQX(32),QQY(32)

      J = 0
 2470 IF (.NOT. (J.LT.N)) GO TO 2480
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      DO 2490 I = 1,JN
         PPX(I) = MAP(1)*PX(I+J1) + MAP(3)*PY(I+J1) + MAP(5)
         QQX(I) = MAP(1)*QX(I+J1) + MAP(3)*QY(I+J1) + MAP(5)
         PPY(I) = MAP(2)*PX(I+J1) + MAP(4)*PY(I+J1) + MAP(6)
         QQY(I) = MAP(2)*QX(I+J1) + MAP(4)*QY(I+J1) + MAP(6)
 2490 CONTINUE
      MASK = -1
      CALL PLTVWV(MAP(7),MAP(11),JN,MASK,PPX,PPY,QQX,QQY)
      CALL PLTVCM(JN,MASK,PPX,PPY,QQX,QQY)
      GO TO 2470

 2480 CONTINUE
      RETURN

      END
      SUBROUTINE PLTDV3(MAP,N,UX,UY,UZ,VX,VY,VZ)
      REAL MAP(*),UX(*),UY(*),UZ(*),VX(*),VY(*),VZ(*)
      DIMENSION TUX(32),TUY(32),TUZ(32),TVX(32),TVY(32),TVZ(32),
     *          TTUX(32),TTUY(32),TTUZ(32),TTVX(32),TTVY(32),TTVZ(32),
     *          V1(3),Q1(3),V2(3),Q2(3)
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      DO 2510 L = 1,3
         V1(L) = MAP(18+L-1) + MAP(15)*MAP(27+L-1)
         Q1(L) = MAP(27+L-1)
         V2(L) = MAP(18+L-1) + MAP(16)*MAP(27+L-1)
         Q2(L) = -MAP(27+L-1)
 2510 CONTINUE
      J = 0
 2530 IF (.NOT. (J.LT.N)) GO TO 2540
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      MASK = -1
      CALL PLTCV3(JN,MASK,UX(J1),UY(J1),UZ(J1),VX(J1),VY(J1),VZ(J1),TUX,
     *            TUY,TUZ,TVX,TVY,TVZ,V1,Q1)
      CALL PLTCV3(JN,MASK,TUX,TUY,TUZ,TVX,TVY,TVZ,TTUX,TTUY,TTUZ,TTVX,
     *            TTVY,TTVZ,V2,Q2)
      IF (MAP(17).EQ.1.) THEN
         DO 2550 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK).NE.0) THEN
               PMS = (TTUX(K)-MAP(18))*MAP(27) +
     *               (TTUY(K)-MAP(19))*MAP(28) +
     *               (TTUZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TUX(K) = R* ((TTUX(K)-MAP(18))*MAP(21)+
     *                  (TTUY(K)-MAP(19))*MAP(22)+
     *                  (TTUZ(K)-MAP(20))*MAP(23))
               TUY(K) = R* ((TTUX(K)-MAP(18))*MAP(24)+
     *                  (TTUY(K)-MAP(19))*MAP(25)+
     *                  (TTUZ(K)-MAP(20))*MAP(26))
               PMS = (TTVX(K)-MAP(18))*MAP(27) +
     *               (TTVY(K)-MAP(19))*MAP(28) +
     *               (TTVZ(K)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TVX(K) = R* ((TTVX(K)-MAP(18))*MAP(21)+
     *                  (TTVY(K)-MAP(19))*MAP(22)+
     *                  (TTVZ(K)-MAP(20))*MAP(23))
               TVY(K) = R* ((TTVX(K)-MAP(18))*MAP(24)+
     *                  (TTVY(K)-MAP(19))*MAP(25)+
     *                  (TTVZ(K)-MAP(20))*MAP(26))
            END IF

 2550    CONTINUE

      ELSE IF (MAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMV2(MAP,JN,MASK,TUX,TUY,TVX,TVY,TTUX,TTUY,TTVX,TTVY)
      CALL PLTVCM(JN,MASK,TTUX,TTUY,TTVX,TTVY)
      GO TO 2530

 2540 CONTINUE
      RETURN

      END
      SUBROUTINE PLTDG2(MAP,N,XV,YV)
      REAL MAP(*)
      INTEGER N
      REAL XV(*),YV(*)
      REAL XWORK(50),YWORK(50),XWORK1(50),YWORK1(50)
      INTEGER NWORK

      AXX = MAP(1)
      AYY = MAP(4)
      AXY = MAP(3)
      AYX = MAP(2)
      BX = MAP(5)
      BY = MAP(6)
      DO 2570 I = 1,N
         XWORK1(I) = AXX*XV(I) + AXY*YV(I) + BX
         YWORK1(I) = AYX*XV(I) + AYY*YV(I) + BY
 2570 CONTINUE
      NWORK = 50
      CALL PLTCG2(N,XWORK1,YWORK1,NWORK,XWORK,YWORK,MAP(7),MAP(9))
      NO = 50
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XWORK1,YWORK1,MAP(9),MAP(11))
      NWORK = 50
      CALL PLTCG2(NO,XWORK1,YWORK1,NWORK,XWORK,YWORK,MAP(11),MAP(13))
      NO = 50
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XWORK1,YWORK1,MAP(13),MAP(7))
      CALL PLTPLY(NO,XWORK1,YWORK1)
      RETURN

      END
      SUBROUTINE PLTMIX(UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      UMAP(21) = -UMAP(21)
      UMAP(24) = -UMAP(24)
      UMAP(27) = -UMAP(27)
      XC = -XC
      UMAP(18) = -UMAP(18)
      RETURN

      END
      SUBROUTINE PLTMIY(UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      UMAP(21+1) = -UMAP(21+1)
      UMAP(24+1) = -UMAP(24+1)
      UMAP(27+1) = -UMAP(27+1)
      YC = -YC
      UMAP(18+1) = -UMAP(18+1)
      RETURN

      END
      SUBROUTINE PLTMIZ(UMAP)
      COMMON /CENBOD/XC,YC,ZC
      REAL UMAP(*)

      UMAP(21+2) = -UMAP(21+2)
      UMAP(24+2) = -UMAP(24+2)
      UMAP(27+2) = -UMAP(27+2)
      ZC = -ZC
      UMAP(18+2) = -UMAP(18+2)
      RETURN

      END
      SUBROUTINE PLTIX2(UMAP)
      REAL UMAP(*)

      UMAP(1+3) = -UMAP(1+3)
      RETURN

      END
      SUBROUTINE PLTIY2(UMAP)
      REAL UMAP(*)

      UMAP(1) = -UMAP(1)
      RETURN

      END
      SUBROUTINE PLTSVD
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      IF (IPOPD.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSVD','Too many calls to PLTSVD.',3)
         RETURN

      END IF

      IPOPD = IPOPD + 1
      DO 2000 I = 1,5
         TDEVP(I,IPOPD) = DEVP(I)
 2000 CONTINUE
      RETURN

      END
      SUBROUTINE PLTSVG
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      IF (IPOPG.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSVG','Too many calls to PLTSVG.',3)
         RETURN

      END IF

      IPOPG = IPOPG + 1
      DO 2020 I = 1,100
         TGRAPH(I,IPOPG) = GRAPHP(I)
 2020 CONTINUE
      RETURN

      END
      SUBROUTINE PLTSVM
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      IF (IPOPM.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSVM','Too many calls to PLTSVM.',3)
         RETURN

      END IF

      IPOPM = IPOPM + 1
      DO 2040 I = 1,11
         TMAPP(I,IPOPM) = MAPP(I)
 2040 CONTINUE
      RETURN

      END
      SUBROUTINE PLTSVT
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      IF (IPOPT.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSVT','Too many calls to PLTSVT.',3)
         RETURN

      END IF

      IPOPT = IPOPT + 1
      DO 2060 I = 1,40
         TTEXTP(I,IPOPT) = TEXTP(I)
 2060 CONTINUE
      RETURN

      END
      SUBROUTINE PLTSVV
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      IF (IPOPV.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSVV','Too many calls to PLTSVV.',3)
         RETURN

      END IF

      IPOPV = IPOPV + 1
      DO 2080 I = 1,5
         TVECTP(I,IPOPV) = VECTP(I)
 2080 CONTINUE
      RETURN

      END
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

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
      SUBROUTINE PLTREG
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      DO 2120 I = 1,6
         GRAPHP(I) = TGRAPH(I,IPOPG)
 2120 CONTINUE
      DO 2140 I = 21,100
         GRAPHP(I) = TGRAPH(I,IPOPG)
 2140 CONTINUE
      IF (IPOPG.NE.1) THEN
         IPOPG = IPOPG - 1
      END IF

      RETURN

      END
      SUBROUTINE PLTREM
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      DO 2160 I = 1,11
         MAPP(I) = TMAPP(I,IPOPM)
 2160 CONTINUE
      IF (IPOPM.NE.1) THEN
         IPOPM = IPOPM - 1
      END IF

      RETURN

      END
      SUBROUTINE PLTRET
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

      DO 2180 I = 1,40
         TEXTP(I) = TTEXTP(I,IPOPT)
 2180 CONTINUE
      IF (TEXTP(35).NE.TTEXTP(35,IPOPT)) THEN
         CALL PLTSTT(1,TEXTP(35))
      END IF

      IF (IPOPT.NE.1) THEN
         IPOPT = IPOPT - 1
      END IF

      RETURN

      END
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
     *       TGRAPH(100,10),TMAPP(5,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM

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
      LOGICAL FUNCTION PLTSTC(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR

      PLTSTC = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSC

      ELSE IF (INDX.EQ.1) THEN
         BTEMP = COLP(1)
         COLP(1) = BUFF(1)
         IF (BUFF(1).EQ.-1.) THEN
            CALL PLTSPC(0.,.706,0.,.706,.1875,0.,0.,1.)
            CALL PLTSPC(.1875,0.,0.,1.,.3708,0.,1.,0.)
            CALL PLTSPC(.3708,0.,1.,0.,.6208,1.,1.,0.)
            CALL PLTSPC(.6208,1.,1.,0.,.8292,1.,.659,0.)
            CALL PLTSPC(.8292,1.,.659,0.,1.,1.,0.,0.)

         ELSE IF (BUFF(1).EQ.0.) THEN
            CALL PLTSPC(0.,1.,0.,0.,.1708,1.,.659,0.)
            CALL PLTSPC(.1708,1.,.659,0.,.3792,1.,1.,0.)
            CALL PLTSPC(.3792,1.,1.,0.,.6292,0.,1.,0.)
            CALL PLTSPC(.6292,0.,1.,0.,.8125,0.,0.,1.)
            CALL PLTSPC(.8125,0.,0.,1.,1.,.706,0.,.706)

         ELSE IF (BUFF(1).EQ.1.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,2),PALETT(2,2),
     *                  PALETT(3,2))

         ELSE IF (BUFF(1).EQ.2.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,3),PALETT(2,3),
     *                  PALETT(3,3))

         ELSE IF (BUFF(1).EQ.3.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,4),PALETT(2,4),
     *                  PALETT(3,4))

         ELSE IF (BUFF(1).EQ.4.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,5),PALETT(2,5),
     *                  PALETT(3,5))

         ELSE IF (BUFF(1).EQ.5.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,6),PALETT(2,6),
     *                  PALETT(3,6))

         ELSE IF (BUFF(1).EQ.6.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,7),PALETT(2,7),
     *                  PALETT(3,7))

         ELSE IF (BUFF(1).EQ.7.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,8),PALETT(2,8),
     *                  PALETT(3,8))

         ELSE IF (BUFF(1).EQ.8.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,9),PALETT(2,9),
     *                  PALETT(3,9))

         ELSE IF (BUFF(1).EQ.9.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,10),PALETT(2,10),
     *                  PALETT(3,10))

         ELSE IF (BUFF(1).EQ.10.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,11),PALETT(2,11),
     *                  PALETT(3,11))

         ELSE IF (BUFF(1).EQ.11.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,12),PALETT(2,12),
     *                  PALETT(3,12))

         ELSE IF (BUFF(1).EQ.12.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,13),PALETT(2,13),
     *                  PALETT(3,13))

         ELSE IF (BUFF(1).EQ.13.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,14),PALETT(2,14),
     *                  PALETT(3,14))

         ELSE IF (BUFF(1).EQ.14.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,15),PALETT(2,15),
     *                  PALETT(3,15))

         ELSE IF (BUFF(1).EQ.15.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,16),PALETT(2,16),
     *                  PALETT(3,16))

         ELSE
            COLP(1) = BTEMP
            CALL CHRRVC(BUFF(1),IERROR,L)
            CALL PLTFLU
            CALL SIORPT('PLTSTC','Illegal buffer '//IERROR(1:L)//
     *                  ' passed to PLTSTC.',2)
            PLTSTC = .FALSE.
         END IF

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTC','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTC = .FALSE.
      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTC(INDX,BUFF)
      CHARACTER*16 IERROR
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
      INTEGER INDX
      REAL BUFF(*)

      PLTGTC = .TRUE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = COLP(1)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTC','Illegal index '//IERROR(1:L)//'.',2)
         PLTGTC = .FALSE.
         RETURN

      END IF

      RETURN

      END
      SUBROUTINE PLTRSC
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

      COLP(2) = MIN(16.,DEVCAP(4))
      COLP(3) = DEVCAP(4) - COLP(2)
      PALETT(1,1) = 0.
      PALETT(2,1) = 0.
      PALETT(3,1) = 0.
      CALL VDSTCO(1,0,PALETT(1,1),0)
      PALETT(1,2) = 1.
      PALETT(2,2) = 0.
      PALETT(2,3) = 0.
      CALL VDSTCO(1,1,PALETT(1,2),0)
      PALETT(1,3) = 0.
      PALETT(2,3) = 1.
      PALETT(3,3) = 0.
      CALL VDSTCO(1,2,PALETT(1,3),0)
      PALETT(1,4) = 1.
      PALETT(2,4) = 1.
      PALETT(3,4) = 0.
      CALL VDSTCO(1,3,PALETT(1,4),0)
      PALETT(1,5) = 0.
      PALETT(2,5) = 0.
      PALETT(3,5) = 1.
      CALL VDSTCO(1,4,PALETT(1,5),0)
      PALETT(1,6) = 1.
      PALETT(2,6) = 0.
      PALETT(3,6) = 1.
      CALL VDSTCO(1,5,PALETT(1,6),0)
      PALETT(1,7) = 0.
      PALETT(2,7) = 1.
      PALETT(3,7) = 1.
      CALL VDSTCO(1,6,PALETT(1,7),0)
      PALETT(1,8) = 1.
      PALETT(2,8) = 1.
      PALETT(3,8) = 1.
      CALL VDSTCO(1,7,PALETT(1,8),0)
      PALETT(1,9) = .4
      PALETT(2,9) = .4
      PALETT(3,9) = .4
      CALL VDSTCO(1,8,PALETT(1,9),0)
      PALETT(1,10) = .7
      PALETT(2,10) = .7
      PALETT(3,10) = .7
      CALL VDSTCO(1,9,PALETT(1,10),0)
      PALETT(1,11) = .225
      PALETT(2,11) = .225
      PALETT(3,11) = .225
      CALL VDSTCO(1,10,PALETT(1,11),0)
      PALETT(1,12) = 1.
      PALETT(2,12) = .35
      PALETT(3,12) = .45
      CALL VDSTCO(1,11,PALETT(1,12),0)
      PALETT(1,13) = .5
      PALETT(2,13) = 1.
      PALETT(3,13) = .8
      CALL VDSTCO(1,12,PALETT(1,13),0)
      PALETT(1,14) = .4
      PALETT(2,14) = .7
      PALETT(3,14) = 1.
      CALL VDSTCO(1,13,PALETT(1,14),0)
      PALETT(1,15) = .706
      PALETT(2,15) = 0.
      PALETT(3,15) = .706
      CALL VDSTCO(1,14,PALETT(1,15),0)
      PALETT(1,16) = 1.
      PALETT(2,16) = .659
      PALETT(3,16) = 0.
      CALL VDSTCO(1,15,PALETT(1,16),0)
      CALL PLTSTC(1,0.)
      RETURN

      END
      LOGICAL FUNCTION PLTSTD(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR
      DATA ZZZFC/-1./,ZZZBC/-1./,ZZZIN/-1./

      PLTSTD = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSD

      ELSE IF (INDX.EQ.1) THEN
         IF (ZZZFC.EQ.BUFF(1)) THEN
            RETURN

         END IF

         ZZZFC = BUFF(1)
         J = NINT(BUFF(1))
         CALL VDSTFC(J)
         DEVP(1) = BUFF(1)

      ELSE IF (INDX.EQ.2) THEN
         IF (ZZZBC.EQ.BUFF(1)) THEN
            RETURN

         END IF

         ZZZBC = BUFF(1)
         CALL VDSTBC(IFIX(BUFF(1)))
         DEVP(2) = BUFF(1)

      ELSE IF (INDX.EQ.3) THEN
         IF (ZZZIN.EQ.BUFF(1)) THEN
            RETURN

         END IF

         ZZZIN = BUFF(1)
         CALL VDSTIN(BUFF(1)/100.)
         DEVP(3) = BUFF(1)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTD','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTD = .FALSE.
         RETURN

      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTD(INDX,BUFF)
      CHARACTER*16 IERROR
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
      INTEGER INDX
      REAL BUFF(*)

      PLTGTD = .TRUE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = DEVP(1)

      ELSE IF (INDX.EQ.2) THEN
         BUFF(1) = DEVP(2)

      ELSE IF (INDX.EQ.3) THEN
         BUFF(1) = DEVP(3)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTD','Illegal index '//IERROR(1:L)//'.',2)
         PLTGTD = .FALSE.
         RETURN

      END IF

      RETURN

      END
      SUBROUTINE PLTRSD
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

      CALL VDSTFC(IFIX(DEFOUT(1)))
      CALL VDSTBC(IFIX(DEFOUT(2)))
      CALL VDSTIN(DEFOUT(3))
      DEVP(1) = DEFOUT(1)
      DEVP(2) = DEFOUT(2)
      DEVP(3) = DEFOUT(3)
      RETURN

      END
      LOGICAL FUNCTION PLTSTG(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR

      PLTSTG = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSG

      ELSE IF (INDX.EQ.1) THEN
         GRAPHP(1) = BUFF(1)

      ELSE IF (INDX.EQ.2) THEN
         GRAPHP(2) = BUFF(1)

      ELSE IF (INDX.EQ.3) THEN
         GRAPHP(3) = BUFF(1)

      ELSE IF (INDX.EQ.4) THEN
         GRAPHP(4) = BUFF(1)

      ELSE IF (INDX.EQ.5) THEN
         GRAPHP(5) = BUFF(1)

      ELSE IF (INDX.EQ.6) THEN
         GRAPHP(38) = BUFF(1)

      ELSE IF (INDX.EQ.7) THEN
         IF (BUFF(1).EQ.0.) THEN
            GRAPHP(6) = 0.

         ELSE
            GRAPHP(6) = 1.
            GRAPHP(47) = BUFF(1) + 4.
         END IF

      ELSE IF (INDX.EQ.8) THEN
         IF (BUFF(1).LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTSTG',
     *                  'Symbol increment must be greater than zero.',2)
            PLTSTG = .FALSE.
            RETURN

         END IF

         GRAPHP(23) = BUFF(1)

      ELSE IF (INDX.EQ.9) THEN
         GRAPHP(21) = BUFF(1)

      ELSE IF (INDX.EQ.10) THEN
         GRAPHP(37) = BUFF(1)

      ELSE IF (INDX.EQ.11) THEN
         IF (BUFF(1).EQ.1.) THEN
            GRAPHP(22) = BUFF(1)

         ELSE IF (BUFF(1).EQ.2.) THEN
            GRAPHP(22) = BUFF(1)

         ELSE IF (BUFF(1).EQ.3.) THEN
            IF (BUFF(2).EQ.BUFF(3)) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG','XMIN cannot be equal to XMAX.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(6).EQ.BUFF(7)) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG','YMIN cannot be equal to YMAX.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(4).EQ.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                  'Number of major x intervals cannot equal zero.'
     *                     ,2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(8).EQ.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                  'Number of major y intervals cannot equal zero.'
     *                     ,2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            GRAPHP(22) = BUFF(1)
            DO 2220 I = 0,7
               GRAPHP(I+24) = BUFF(I+2)
 2220       CONTINUE

         ELSE IF (BUFF(1).EQ.4.) THEN
            IF (BUFF(4).EQ.BUFF(3)) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                     'XMAX cannot equal first nice X number.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(5).EQ.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG','X interval cannot equal zero.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(9).EQ.BUFF(8)) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                     'YMAX cannot equal first nice Y number.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(10).EQ.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG','Y interval cannot equal zero.',2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF ((BUFF(3).LT.BUFF(2).AND.BUFF(4).GT.BUFF(2)) .OR.
     *          (BUFF(3).GT.BUFF(2).AND.BUFF(4).LT.BUFF(3))) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                  'Invalid specification of XMIN, XNICE and XMAX.'
     *                     ,2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF ((BUFF(8).LT.BUFF(7).AND.BUFF(9).GT.BUFF(7)) .OR.
     *          (BUFF(8).GT.BUFF(7).AND.BUFF(9).LT.BUFF(8))) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                  'Invalid specification of YMIN, YNICE and YMAX.'
     *                     ,2)
               PLTSTG = .FALSE.
               RETURN

            END IF

            IF (BUFF(4).LT.BUFF(2) .AND. BUFF(5).GT.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                     'Setting X interval negative as XMAX < XMIN.'
     *                     ,2)
               BUFF(5) = -BUFF(5)
            END IF

            IF (BUFF(9).LT.BUFF(7) .AND. BUFF(10).GT.0.) THEN
               CALL PLTFLU
               CALL SIORPT('PLTSTG',
     *                     'Setting Y interval negative as YMAX < YMIN.'
     *                     ,2)
               BUFF(10) = -BUFF(10)
            END IF

            GRAPHP(22) = BUFF(1)
            DO 2240 I = 0,9
               GRAPHP(I+78) = BUFF(I+2)
 2240       CONTINUE

         ELSE
            CALL CHRRVC(BUFF(1),IERROR,L)
            CALL PLTFLU
            CALL SIORPT('PLTSTG','Illegal buffer value '//IERROR(1:L)//
     *                  ' in graph scaling type.',2)
            PLTSTG = .FALSE.
            RETURN

         END IF

      ELSE IF (INDX.EQ.12) THEN
         GRAPHP(32) = BUFF(1)

      ELSE IF (INDX.EQ.13) THEN
         GRAPHP(91) = BUFF(1)

      ELSE IF (INDX.EQ.14) THEN
         GRAPHP(90) = BUFF(1)

      ELSE IF (INDX.EQ.15) THEN
         GRAPHP(35) = BUFF(1)

      ELSE IF (INDX.EQ.16) THEN
         GRAPHP(36) = BUFF(1)

      ELSE IF (INDX.EQ.17) THEN
         GRAPHP(39) = BUFF(1)

      ELSE IF (INDX.EQ.18) THEN
         GRAPHP(40) = BUFF(1)

      ELSE IF (INDX.EQ.19) THEN
         GRAPHP(41) = BUFF(1)

      ELSE IF (INDX.EQ.20) THEN
         GRAPHP(42) = BUFF(1)

      ELSE IF (INDX.EQ.21) THEN
         GRAPHP(92) = BUFF(1)

      ELSE IF (INDX.EQ.22) THEN
         GRAPHP(44) = BUFF(1)

      ELSE IF (INDX.EQ.23) THEN
         GRAPHP(45) = BUFF(1)

      ELSE IF (INDX.EQ.47) THEN
         GRAPHP(88) = BUFF(1)

      ELSE IF (INDX.EQ.48) THEN
         GRAPHP(89) = BUFF(1)

      ELSE IF (INDX.EQ.24) THEN
         GRAPHP(46) = BUFF(1)

      ELSE IF (INDX.EQ.25) THEN
         GRAPHP(48) = BUFF(1)

      ELSE IF (INDX.EQ.26) THEN
         GRAPHP(49) = BUFF(1)

      ELSE IF (INDX.EQ.27) THEN
         DO 2260 I = 0,13
            GRAPHP(7+I) = BUFF(I+1)
 2260    CONTINUE

      ELSE IF (INDX.EQ.28) THEN
         GRAPHP(62) = BUFF(1)

      ELSE IF (INDX.EQ.29) THEN
         GRAPHP(63) = BUFF(1)

      ELSE IF (INDX.EQ.30) THEN
         GRAPHP(64) = BUFF(1)

      ELSE IF (INDX.EQ.31) THEN
         GRAPHP(65) = BUFF(1)

      ELSE IF (INDX.EQ.32) THEN
         GRAPHP(66) = BUFF(1)

      ELSE IF (INDX.EQ.33) THEN
         GRAPHP(67) = BUFF(1)

      ELSE IF (INDX.EQ.34) THEN
         GRAPHP(68) = BUFF(1)

      ELSE IF (INDX.EQ.35) THEN
         GRAPHP(69) = BUFF(1)

      ELSE IF (INDX.EQ.36) THEN
         GRAPHP(70) = BUFF(1)

      ELSE IF (INDX.EQ.37) THEN
         GRAPHP(71) = BUFF(1)

      ELSE IF (INDX.EQ.38) THEN
         GRAPHP(72) = BUFF(1)

      ELSE IF (INDX.EQ.39) THEN
         GRAPHP(73) = BUFF(1)

      ELSE IF (INDX.EQ.43) THEN
         GRAPHP(74) = BUFF(1)

      ELSE IF (INDX.EQ.44) THEN
         GRAPHP(75) = BUFF(1)

      ELSE IF (INDX.EQ.45) THEN
         GRAPHP(76) = BUFF(1)

      ELSE IF (INDX.EQ.46) THEN
         GRAPHP(77) = BUFF(1)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTG','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTG = .FALSE.
         RETURN

      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTG(INDX,BUFF)
      CHARACTER*16 IERROR
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
      REAL BUFF(*)

      PLTGTG = .TRUE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = GRAPHP(1)

      ELSE IF (INDX.EQ.2) THEN
         BUFF(1) = GRAPHP(2)

      ELSE IF (INDX.EQ.3) THEN
         BUFF(1) = GRAPHP(3)

      ELSE IF (INDX.EQ.4) THEN
         BUFF(1) = GRAPHP(4)

      ELSE IF (INDX.EQ.5) THEN
         BUFF(1) = GRAPHP(5)

      ELSE IF (INDX.EQ.6) THEN
         BUFF(1) = GRAPHP(38)

      ELSE IF (INDX.EQ.7) THEN
         IF (GRAPHP(6).EQ.1.) THEN
            BUFF(1) = GRAPHP(47) - 4.

         ELSE
            BUFF(1) = 0.
         END IF

      ELSE IF (INDX.EQ.8) THEN
         BUFF(1) = GRAPHP(23)

      ELSE IF (INDX.EQ.9) THEN
         BUFF(1) = GRAPHP(21)

      ELSE IF (INDX.EQ.10) THEN
         BUFF(1) = GRAPHP(37)

      ELSE IF (INDX.EQ.11) THEN
         IF (BUFF(1).EQ.3.) THEN
            BUFF(2) = GRAPHP(24)
            BUFF(3) = GRAPHP(25)
            BUFF(4) = GRAPHP(26)
            BUFF(5) = GRAPHP(27)
            BUFF(6) = GRAPHP(28)
            BUFF(7) = GRAPHP(29)
            BUFF(8) = GRAPHP(30)
            BUFF(9) = GRAPHP(31)
         END IF

         IF (BUFF(1).EQ.4.) THEN
            BUFF(2) = GRAPHP(78)
            BUFF(3) = GRAPHP(79)
            BUFF(4) = GRAPHP(80)
            BUFF(5) = GRAPHP(81)
            BUFF(6) = GRAPHP(82)
            BUFF(7) = GRAPHP(83)
            BUFF(8) = GRAPHP(84)
            BUFF(9) = GRAPHP(85)
            BUFF(10) = GRAPHP(86)
            BUFF(11) = GRAPHP(87)
         END IF

         BUFF(1) = GRAPHP(22)

      ELSE IF (INDX.EQ.12) THEN
         BUFF(1) = GRAPHP(32)

      ELSE IF (INDX.EQ.13) THEN
         BUFF(1) = GRAPHP(91)

      ELSE IF (INDX.EQ.14) THEN
         BUFF(1) = GRAPHP(90)

      ELSE IF (INDX.EQ.15) THEN
         BUFF(1) = GRAPHP(35)

      ELSE IF (INDX.EQ.16) THEN
         BUFF(1) = GRAPHP(36)

      ELSE IF (INDX.EQ.17) THEN
         BUFF(1) = GRAPHP(39)

      ELSE IF (INDX.EQ.18) THEN
         BUFF(1) = GRAPHP(40)

      ELSE IF (INDX.EQ.19) THEN
         BUFF(1) = GRAPHP(41)

      ELSE IF (INDX.EQ.20) THEN
         BUFF(1) = GRAPHP(42)

      ELSE IF (INDX.EQ.21) THEN
         BUFF(1) = GRAPHP(92)

      ELSE IF (INDX.EQ.22) THEN
         BUFF(1) = GRAPHP(44)

      ELSE IF (INDX.EQ.23) THEN
         BUFF(1) = GRAPHP(45)

      ELSE IF (INDX.EQ.47) THEN
         BUFF(1) = GRAPHP(88)

      ELSE IF (INDX.EQ.48) THEN
         BUFF(1) = GRAPHP(89)

      ELSE IF (INDX.EQ.24) THEN
         BUFF(1) = GRAPHP(46)

      ELSE IF (INDX.EQ.27) THEN
         DO 2280 I = 0,13
            BUFF(I+1) = GRAPHP(7+I)
 2280    CONTINUE

      ELSE IF (INDX.EQ.28) THEN
         BUFF(1) = GRAPHP(62)

      ELSE IF (INDX.EQ.29) THEN
         BUFF(1) = GRAPHP(63)

      ELSE IF (INDX.EQ.30) THEN
         BUFF(1) = GRAPHP(64)

      ELSE IF (INDX.EQ.31) THEN
         BUFF(1) = GRAPHP(65)

      ELSE IF (INDX.EQ.32) THEN
         BUFF(1) = GRAPHP(66)

      ELSE IF (INDX.EQ.33) THEN
         BUFF(1) = GRAPHP(67)

      ELSE IF (INDX.EQ.34) THEN
         BUFF(1) = GRAPHP(68)

      ELSE IF (INDX.EQ.35) THEN
         BUFF(1) = GRAPHP(69)

      ELSE IF (INDX.EQ.36) THEN
         BUFF(1) = GRAPHP(70)

      ELSE IF (INDX.EQ.37) THEN
         BUFF(1) = GRAPHP(71)

      ELSE IF (INDX.EQ.38) THEN
         BUFF(1) = GRAPHP(72)

      ELSE IF (INDX.EQ.39) THEN
         BUFF(1) = GRAPHP(73)

      ELSE IF (INDX.EQ.43) THEN
         BUFF(1) = GRAPHP(74)

      ELSE IF (INDX.EQ.44) THEN
         BUFF(1) = GRAPHP(75)

      ELSE IF (INDX.EQ.45) THEN
         BUFF(1) = GRAPHP(76)

      ELSE IF (INDX.EQ.46) THEN
         BUFF(1) = GRAPHP(77)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTG','Illegal index '//IERROR(1:L)//'.',2)
         PLTGTG = .FALSE.
         RETURN

      END IF

      RETURN

      END
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
      LOGICAL FUNCTION PLTSTM(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR
      REAL LEFT
      CHARACTER*6 SUBNAM
      DATA SUBNAM/'PLTSTM'/

      PLTSTM = .FALSE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSM

      ELSE IF (INDX.EQ.1) THEN
         IF (BUFF(1).EQ.0.) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,'You cannot set the X scale factor to 0.'
     *                  ,2)
            RETURN

         END IF

         MAPP(1) = BUFF(1)

      ELSE IF (INDX.EQ.2) THEN
         IF (BUFF(1).EQ.0.) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,'You cannot set the Y scale factor to 0.'
     *                  ,2)
            RETURN

         END IF

         MAPP(2) = BUFF(1)

      ELSE IF (INDX.EQ.3) THEN
         MAPP(3) = BUFF(1)

      ELSE IF (INDX.EQ.4) THEN
         MAPP(4) = BUFF(1)

      ELSE IF (INDX.EQ.5) THEN
         MAPP(5) = BUFF(1)

      ELSE IF (INDX.EQ.6) THEN
         LEFT = BUFF(1)
         RIGHT = BUFF(2)
         BOTTOM = BUFF(3)
         TOP = BUFF(4)
         IF (LEFT.GE.RIGHT) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,
     *  'You cannot specify the left viewport edge >= to the right edge'
     *                  ,2)
            RETURN

         END IF

         IF (TOP.LE.BOTTOM) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,
     *  'You cannot specify the top viewport edge <= to the bottom edge'
     *                  ,2)
            RETURN

         END IF

         IF (TOP.LT.0 .OR. TOP.GT.DEVP(5)) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,'Top viewport specification out of range'
     *                  ,2)
            RETURN

         END IF

         IF (BOTTOM.LT.0 .OR. BOTTOM.GT.DEVP(5)) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,
     *                  'Bottom viewport specification out of range',2)
            RETURN

         END IF

         IF (LEFT.LT.0 .OR. LEFT.GT.DEVP(4)) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,
     *                  'Left viewport specification out of range',2)
            RETURN

         END IF

         IF (RIGHT.LT.0 .OR. RIGHT.GT.DEVP(4)) THEN
            CALL PLTFLU
            CALL SIORPT(SUBNAM,
     *                  'Right viewport specification out of range',2)
            RETURN

         END IF

         MAPP(6) = LEFT
         MAPP(8) = RIGHT
         MAPP(7) = BOTTOM
         MAPP(9) = TOP

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTM','Illegal index '//IERROR(1:L)//'.',2)
         RETURN

      END IF

      PLTSTM = .TRUE.
      MAPP(10) = 0.
      IF (MAPP(1).NE.1. .OR. MAPP(2).NE.1. .OR. MAPP(3).NE.0. .OR.
     *    MAPP(4).NE.0. .OR. MAPP(5).NE.0.) THEN
         MAPP(10) = 1.
      END IF

      IF (MAPP(6).NE.0. .OR. MAPP(8).NE.DEVP(4) .OR. MAPP(7).NE.0. .OR.
     *    MAPP(9).NE.DEVP(5)) THEN
         MAPP(10) = 1.
      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTM(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR

      PLTGTM = .FALSE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = MAPP(1)

      ELSE IF (INDX.EQ.2) THEN
         BUFF(1) = MAPP(2)

      ELSE IF (INDX.EQ.3) THEN
         BUFF(1) = MAPP(3)

      ELSE IF (INDX.EQ.4) THEN
         BUFF(1) = MAPP(4)

      ELSE IF (INDX.EQ.5) THEN
         BUFF(1) = MAPP(5)

      ELSE IF (INDX.EQ.6) THEN
         BUFF(1) = MAPP(6)
         BUFF(2) = MAPP(8)
         BUFF(3) = MAPP(7)
         BUFF(4) = MAPP(9)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTM','Illegal index '//IERROR(1:L)//'.',2)
         RETURN

      END IF

      PLTGTM = .TRUE.
      RETURN

      END
      SUBROUTINE PLTRSM
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

      MAPP(1) = 1.
      MAPP(2) = 1.
      MAPP(3) = 0.
      MAPP(4) = 0.
      MAPP(5) = 0.
      MAPP(6) = 0.
      MAPP(8) = DEVP(4)
      MAPP(7) = 0.
      MAPP(9) = DEVP(5)
      MAPP(10) = 0.
      RETURN

      END
      LOGICAL FUNCTION PLTSTT(INDX,BUFF)
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
      DIMENSION BUFF(*),DUMMY(7)
      CHARACTER*16 IERROR

      PLTSTT = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRST

      ELSE IF (INDX.EQ.1) THEN
         CALL VDSTCS(BUFF(1))
         CALL VDIQOS(DUMMY)
         TEXTP(35) = DUMMY(6)
         TEXTP(36) = DUMMY(7)

      ELSE IF (INDX.EQ.2) THEN
         TEXTP(1) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.3) THEN
         TEXTP(2) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.4) THEN
         TEXTP(3) = BUFF(1)
         CALL PLTITM

      ELSE IF (INDX.EQ.5) THEN
         DO 2300 I = 20,27
            TEXTP(I) = BUFF(I-19)
 2300    CONTINUE

      ELSE IF (INDX.EQ.6) THEN
         TEXTP(30) = BUFF(1)

      ELSE IF (INDX.EQ.7) THEN
         TEXTP(31) = BUFF(1)

      ELSE IF (INDX.EQ.8) THEN
         TEXTP(32) = BUFF(1)

      ELSE IF (INDX.EQ.9) THEN
         TEXTP(33) = BUFF(1)

      ELSE IF (INDX.EQ.10) THEN
         TEXTP(34) = BUFF(1)

      ELSE IF (INDX.EQ.11) THEN
         TEXTP(37) = BUFF(1)

      ELSE IF (INDX.EQ.12) THEN
         IF (BUFF(1).EQ.1.) THEN
            CALL PLTFNT('ROMFNT')
            TEXTP(38) = 1.

         ELSE IF (BUFF(1).EQ.2.) THEN
            CALL PLTFNT('STKFNT')
            TEXTP(38) = 2.

         ELSE IF (BUFF(1).EQ.3.) THEN
            CALL PLTFNT('SSRFNT')
            TEXTP(38) = 3.

         ELSE
            CALL CHRRVC(BUFF(1),IERROR,L)
            CALL PLTFLU
            CALL SIORPT('PLTSTT','Illegal buffer '//IERROR(1:L)//
     *                  ' passed to PLTSTT.',2)
            PLTSTT = .FALSE.
         END IF

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTT','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTT = .FALSE.
         RETURN

      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTT(INDX,BUFF)
      CHARACTER*16 IERROR
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
      REAL BUFF(*)

      PLTGTT = .TRUE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = TEXTP(35)

      ELSE IF (INDX.EQ.2) THEN
         BUFF(1) = TEXTP(1)

      ELSE IF (INDX.EQ.3) THEN
         BUFF(1) = TEXTP(2)

      ELSE IF (INDX.EQ.4) THEN
         BUFF(1) = TEXTP(3)

      ELSE IF (INDX.EQ.5) THEN
         DO 2320 I = 20,27
            BUFF(I-19) = TEXTP(I)
 2320    CONTINUE

      ELSE IF (INDX.EQ.6) THEN
         BUFF(1) = TEXTP(30)

      ELSE IF (INDX.EQ.7) THEN
         BUFF(1) = TEXTP(31)

      ELSE IF (INDX.EQ.8) THEN
         BUFF(1) = TEXTP(32)

      ELSE IF (INDX.EQ.9) THEN
         BUFF(1) = TEXTP(33)

      ELSE IF (INDX.EQ.10) THEN
         BUFF(1) = TEXTP(34)

      ELSE IF (INDX.EQ.11) THEN
         BUFF(1) = TEXTP(37)

      ELSE IF (INDX.EQ.12) THEN
         BUFF(1) = TEXTP(38)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTT','Illegal index '//IERROR(1:L)//'.',2)
         PLTGTT = .FALSE.
         RETURN

      END IF

      RETURN

      END
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
      SUBROUTINE PLTITM
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
      DATA DPR/57.29577951/

      COSD(X) = COS(X/DPR)
      SIND(X) = SIN(X/DPR)
      TEXTP(12) = COSD(TEXTP(2)+TEXTP(3))
      TEXTP(13) = SIND(TEXTP(2)+TEXTP(3))
      TEXTP(14) = TEXTP(12)*TEXTP(1)
      TEXTP(15) = TEXTP(13)*TEXTP(1)
      TEXTP(16) = -TEXTP(13)*TEXTP(1)
      TEXTP(17) = TEXTP(12)*TEXTP(1)
      TEXTP(28) = COSD(TEXTP(2))
      TEXTP(29) = SIND(TEXTP(2))
      RETURN

      END
      LOGICAL FUNCTION PLTSTV(INDX,BUFF)
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
      DIMENSION BUFF(*)
      CHARACTER*16 IERROR
      DATA ZZZLS/-1./

      PLTSTV = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSV

      ELSE IF (INDX.EQ.1) THEN
         IF (ZZZLS.EQ.BUFF(1)) THEN
            RETURN

         END IF

         ZZZLS = BUFF(1)
         VECTP(1) = BUFF(1)
         IF (BUFF(1).NE.0.) THEN
            IBUFFT = IFIX(BUFF(1)) - 1
            CALL VDSTLS(IBUFFT)
         END IF

      ELSE IF (INDX.EQ.2) THEN
         CALL VDSTLW(BUFF(1)/1000.)
         VECTP(2) = BUFF(1)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTV','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTV = .FALSE.
         RETURN

      END IF

      RETURN

      END
      LOGICAL FUNCTION PLTGTV(INDX,BUFF)
      CHARACTER*16 IERROR
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
      DIMENSION BUFF(*)

      PLTGTV = .TRUE.
      IF (INDX.EQ.1) THEN
         BUFF(1) = VECTP(1)

      ELSE IF (INDX.EQ.2) THEN
         BUFF(1) = VECTP(2)

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTGTV','Illegal index '//IERROR(1:L)//'.',2)
         PLTGTV = .FALSE.
         RETURN

      END IF

      RETURN

      END
      SUBROUTINE PLTRSV
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

      IBUFF = DEFOUT(4)
      CALL VDSTLS(IBUFF)
      VECTP(1) = DEFOUT(4) + 1.
      CALL VDSTLW(DEFOUT(5))
      VECTP(2) = DEFOUT(5)
      RETURN

      END
      SUBROUTINE MP2PG(N,XV,YV,NO,XVO,YVO)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),XV(*),YV(*),XVO(*),YVO(*)

      IF (N.GT.32) THEN
         CALL PLTFLU
         CALL SIORPT('MP2PG','Too many vertices specified',2)
         RETURN

      END IF

      NOSAVE = NO
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      IF (NOCLIP) THEN
         CALL MPMUL2(N,XV,YV,MVP,TARR1,TARR2,TARR3,TARR4)
         NT = N

      ELSE
         CALL MPMUL2(N,XV,YV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2000 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            NT = 32
            CALL PLTCG2(N,TARR1,TARR2,NT,TARR5,TARR6,C1,C2)
 2000    CONTINUE
         IF (NT.GT.N) THEN
            DO 2020 J = N + 1,NT
               TARR4(J) = 1.
 2020       CONTINUE
         END IF

         CALL MPMUL4(NT,-1,TARR5,TARR6,TARR3,TARR4,VP,TARR1,TARR2,TARR3,
     *               TARR4)
      END IF

      DO 2040 I = 1,NT
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2040 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NO = NOSAVE
      CALL PLTVWG(C1,C2,NT,TARR1,TARR2,TARR3,NO,XVO,YVO,TARR3)
      RETURN

      END
      SUBROUTINE MP2PT(N,X0,Y0,PX,PY,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*),PX(*),PY(*),MASK(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      KM = 0
      J = 0
 2060 IF (.NOT. (J.LT.N)) GO TO 2070
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2080 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCP2(JN,MASK(KM),TARR1,TARR2,C1,C2)
 2080    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      DO 2100 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2100 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK(KM),TARR1,TARR2)
      DO 2120 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
 2120 CONTINUE
      GO TO 2060

 2070 CONTINUE
      RETURN

      END
      SUBROUTINE MP2VC(N,X0,Y0,X1,Y1,PX,PY,QX,QY,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*),X1(*),Y1(*),PX(*),PY(*),QX(*),
     *          QY(*),MASK(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      KM = 0
      J = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MVP,TARR5,TARR6,TARR7,TARR8)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MODEL,TARR5,TARR6,TARR7,TARR8)
         DO 2160 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCV2(JN,MASK(KM),TARR1,TARR2,TARR5,TARR6,TARR1,TARR2,
     *                  TARR5,TARR6,C1,C2)
 2160    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK(KM),TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      DO 2180 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
 2180 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK(KM),TARR1,TARR2,TARR5,TARR6)
      DO 2200 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
         QX(J1+I) = TARR5(I)
         QY(J1+I) = TARR6(I)
 2200 CONTINUE
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
      SUBROUTINE MP3PG(NV,XV,YV,ZV,NO,XVO,YVO,ZVO)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),XV(NV),YV(NV),ZV(NV),XVO(NO),YVO(NO),ZVO(NO)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      NOSAVE = NO
      NO = 0
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(NV,XV,YV,ZV,MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL3(NV,XV,YV,ZV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2220 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(NV,MASK,TARR1,TARR2,TARR3,C1,C2)
            IF (MASK.NE.-1) THEN
               GO TO 2220

            END IF

 2220    CONTINUE
         CALL MPMUL4(NV,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,NV,MASK,TARR4)
      IF (MASK.NE.-1) THEN
         RETURN

      END IF

      DO 2240 I = 1,NV
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR3(I) = TARR3(I)/TARR4(I)
 2240 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NO = NOSAVE
      CALL PLTVWG(C1,C2,NV,TARR1,TARR2,TARR3,NO,XVO,YVO,ZVO)
      RETURN

      END
      SUBROUTINE MP3PT(N,X0,Y0,Z0,PX,PY,PZ,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION PX(*),PY(*),PZ(*),MASK(*),C1(3),C2(3),X0(*),Y0(*),Z0(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      KM = 0
      J = 0
 2260 IF (.NOT. (J.LT.N)) GO TO 2270
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         DO 2280 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(JN,MASK(KM),TARR1,TARR2,TARR3,C1,C2)
 2280    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,JN,MASK(KM),TARR4)
      DO 2300 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR3(I) = TARR3(I)/TARR4(I)
 2300 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK(KM),TARR1,TARR2)
      DO 2320 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
         PZ(J1+I) = TARR3(I)
 2320 CONTINUE
      GO TO 2260

 2270 CONTINUE
      RETURN

      END
      SUBROUTINE MP3VC(N,X0,Y0,Z0,X1,Y1,Z1,PX,PY,PZ,QX,QY,QZ,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION PX(*),PY(*),PZ(*),QX(*),QY(*),QZ(*),MASK(*)
      DIMENSION C1(3),C2(3),X0(*),Y0(*),Z0(*),X1(*),Y1(*),Z1(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      J = 0
      KM = 0
 2340 IF (.NOT. (J.LT.N)) GO TO 2350
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MVP,TARR5,TARR6,
     *               TARR7,TARR8)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MODEL,TARR5,TARR6,
     *               TARR7,TARR8)
         MASK(KM) = -1
         DO 2360 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCV3(JN,MASK(KM),TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,
     *                  TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,C1,C2)
 2360    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK(KM),TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      CALL PLTZCV(CPNEAR,CPFAR,JN,MASK(KM),TARR1,TARR2,TARR4,TARR5,
     *            TARR6,TARR8)
      DO 2380 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR3(I) = TARR3(I)/TARR4(I)
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
         TARR7(I) = TARR7(I)/TARR8(I)
 2380 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK(KM),TARR1,TARR2,TARR5,TARR6)
      DO 2400 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
         PZ(J1+I) = TARR3(I)
         QX(J1+I) = TARR5(I)
         QY(J1+I) = TARR6(I)
         QZ(J1+I) = TARR7(I)
 2400 CONTINUE
      GO TO 2340

 2350 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION MPCLP2(N,C1X,C1Y,C2X,C2Y)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      DIMENSION C1X(*),C1Y(*),C2X(*),C2Y(*)
      PARAMETER (SUBNAM='MPCLP2')

      MPCLP2 = .FALSE.
      IF (N.GT.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'Too many clipping lines specified; max is 10',2)
         RETURN

      END IF

      IF (N.LT.0) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'You cannot specify less than zero clipping lines',
     *               2)
         RETURN

      END IF

      MPCLP2 = .TRUE.
      NCPLIN = N
      DO 2420 I = 1,N
         CPLINE(1,1,I) = C1X(I)
         CPLINE(1,2,I) = C1Y(I)
         CPLINE(2,1,I) = C2X(I)
         CPLINE(2,2,I) = C2Y(I)
 2420 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION MPCLP3(N,PX,PY,PZ,VX,VY,VZ)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      DIMENSION PX(*),PY(*),PZ(*),VX(*),VY(*),VZ(*)
      PARAMETER (SUBNAM='MPCLP3')

      MPCLP3 = .FALSE.
      IF (N.GT.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'Too many clipping planes specified; max is 10',2)
         RETURN

      END IF

      IF (N.LT.0) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'You cannot specify less than zero clipping planes'
     *               ,2)
         RETURN

      END IF

      MPCLP3 = .TRUE.
      NCPLAN = N
      DO 2440 I = 1,N
         CPPLAN(1,1,I) = PX(I)
         CPPLAN(1,2,I) = PY(I)
         CPPLAN(1,3,I) = PZ(I)
         CPPLAN(2,1,I) = VX(I)
         CPPLAN(2,2,I) = VY(I)
         CPPLAN(2,3,I) = VZ(I)
 2440 CONTINUE
      RETURN

      END
      SUBROUTINE MPD2PG(N,XV,YV,MODE)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),XV(*),YV(*)
      CHARACTER MODE* (*),TMODE*1

      IF (N.GT.32) THEN
         CALL PLTFLU
         CALL SIORPT('MPD2PG','Too many vertices specified',2)
         RETURN

      END IF

      TMODE = MODE
      CALL CHRDN(TMODE,TMODE)
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      IF (NOCLIP) THEN
         CALL MPMUL2(N,XV,YV,MVP,TARR1,TARR2,TARR3,TARR4)
         NO = N

      ELSE
         CALL MPMUL2(N,XV,YV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2460 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            NO = 32
            CALL PLTCG2(N,TARR1,TARR2,NO,TARR5,TARR6,C1,C2)
 2460    CONTINUE
         IF (NO.GT.N) THEN
            DO 2480 J = N + 1,NO
               TARR4(J) = 1.
 2480       CONTINUE
         END IF

         CALL MPMUL4(NO,-1,TARR5,TARR6,TARR3,TARR4,VP,TARR1,TARR2,TARR3,
     *               TARR4)
      END IF

      DO 2500 I = 1,NO
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2500 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NVO = 32
      CALL PLTVWG(C1,C2,NO,TARR1,TARR2,TARR3,NVO,TARR5,TARR6,TARR7)
      IF (TMODE.EQ.'s') THEN
         CALL PLTPLY(NVO,TARR5,TARR6)

      ELSE IF (TMODE.EQ.'o') THEN
         CALL PLTMOV(TARR5(1),TARR6(1))
         DO 2520 J = 2,NVO
            CALL PLTDRW(TARR5(J),TARR6(J))
 2520    CONTINUE
         CALL PLTDRW(TARR5(1),TARR6(1))

      ELSE
         CALL PLTFLU
         CALL SIORPT('MPD2PG','Unrecognized drawing mode: '//TMODE,2)
      END IF

      RETURN

      END
      SUBROUTINE MPD2PT(N,X0,Y0)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      J = 0
 2540 IF (.NOT. (J.LT.N)) GO TO 2550
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2560 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCP2(JN,MASK,TARR1,TARR2,C1,C2)
 2560    CONTINUE
         CALL MPMUL4(JN,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      DO 2580 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2580 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK,TARR1,TARR2)
      CALL PLTPTM(JN,MASK,TARR1,TARR2)
      GO TO 2540

 2550 CONTINUE
      RETURN

      END
      SUBROUTINE MPD2SY(N,X0,Y0,SYM)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*(*) SYM
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      J = 0
 2600 IF (.NOT. (J.LT.N)) GO TO 2610
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2620 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCP2(JN,MASK,TARR1,TARR2,C1,C2)
 2620    CONTINUE
         CALL MPMUL4(JN,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      DO 2640 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2640 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK,TARR1,TARR2)
      CALL PLTSBM(JN,MASK,TARR1,TARR2,SYM)
      GO TO 2600

 2610 CONTINUE
      RETURN

      END
      SUBROUTINE MPD2VC(N,X0,Y0,X1,Y1)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*),X1(*),Y1(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      J = 0
 2660 IF (.NOT. (J.LT.N)) GO TO 2670
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MVP,TARR5,TARR6,TARR7,TARR8)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MODEL,TARR5,TARR6,TARR7,TARR8)
         DO 2680 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCV2(JN,MASK,TARR1,TARR2,TARR5,TARR6,TARR1,TARR2,
     *                  TARR5,TARR6,C1,C2)
 2680    CONTINUE
         CALL MPMUL4(JN,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK,TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      DO 2700 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
 2700 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK,TARR1,TARR2,TARR5,TARR6)
      CALL PLTVCM(JN,MASK,TARR1,TARR2,TARR5,TARR6)
      GO TO 2660

 2670 CONTINUE
      RETURN

      END
      SUBROUTINE MPD3PG(NV,XV,YV,ZV,MODE)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),XV(NV),YV(NV),ZV(NV),XVO(50),YVO(50)
      CHARACTER MODE* (*),TMODE*1

      TMODE = MODE
      CALL CHRDN(TMODE,TMODE)
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(NV,XV,YV,ZV,MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL3(NV,XV,YV,ZV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2720 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(NV,MASK,TARR1,TARR2,TARR3,C1,C2)
            IF (MASK.NE.-1) THEN
               GO TO 2720

            END IF

 2720    CONTINUE
         CALL MPMUL4(NV,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,NV,MASK,TARR4)
      IF (MASK.NE.-1) THEN
         RETURN

      END IF

      DO 2740 I = 1,NV
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2740 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NO = 50
      CALL PLTVWG(C1,C2,NV,TARR1,TARR2,TARR3,NO,XVO,YVO,TARR4)
      IF (TMODE.EQ.'s') THEN
         CALL PLTPLY(NO,XVO,YVO)

      ELSE IF (TMODE.EQ.'o') THEN
         CALL PLTMOV(XVO(1),YVO(1))
         DO 2760 J = 2,NO
            CALL PLTDRW(XVO(J),YVO(J))
 2760    CONTINUE
         CALL PLTDRW(XVO(1),YVO(1))

      ELSE
         CALL PLTFLU
         CALL SIORPT('MPD3PG','Unrecognized drawing mode: '//TMODE,2)
         RETURN

      END IF

      RETURN

      END
      SUBROUTINE MPD3PT(N,X0,Y0,Z0)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),X0(*),Y0(*),Z0(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      J = 0
 2780 IF (.NOT. (J.LT.N)) GO TO 2790
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         DO 2800 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(JN,MASK,TARR1,TARR2,TARR3,C1,C2)
 2800    CONTINUE
         CALL MPMUL4(JN,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,JN,MASK,TARR4)
      DO 2820 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2820 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWP(C1,C2,JN,MASK,TARR1,TARR2)
      CALL PLTPTM(JN,MASK,TARR1,TARR2)
      GO TO 2780

 2790 CONTINUE
      RETURN

      END
      SUBROUTINE MPD3VC(N,X0,Y0,Z0,X1,Y1,Z1)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),X0(*),Y0(*),Z0(*),X1(*),Y1(*),Z1(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      J = 0
      KM = 0
 2840 IF (.NOT. (J.LT.N)) GO TO 2850
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MVP,TARR5,TARR6,
     *               TARR7,TARR8)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MODEL,TARR5,TARR6,
     *               TARR7,TARR8)
         MASK = -1
         DO 2860 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCV3(JN,MASK,TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,
     *                  TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,C1,C2)
 2860    CONTINUE
         CALL MPMUL4(JN,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK,TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      CALL PLTZCV(CPNEAR,CPFAR,JN,MASK,TARR1,TARR2,TARR4,TARR5,TARR6,
     *            TARR8)
      DO 2880 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
 2880 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK,TARR1,TARR2,TARR5,TARR6)
      CALL PLTVCM(JN,MASK,TARR1,TARR2,TARR5,TARR6)
      GO TO 2840

 2850 CONTINUE
      RETURN

      END
      SUBROUTINE MPGETM(TMODEL,TVIEW,TPROJ,TVWPT)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      DIMENSION TMODEL(4,4),TVIEW(4,4),TPROJ(4,4),TVWPT(4)

      DO 2900 I = 1,4
         DO 2920 J = 1,4
            TMODEL(J,I) = MODEL(J,I)
            TPROJ(J,I) = PROJ(J,I)
            TVIEW(J,I) = VIEW(J,I)
 2920    CONTINUE
         TVWPT(I) = VWPORT(I)
 2900 CONTINUE
      RETURN

      END
      SUBROUTINE MPINIT()
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
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

      CALL MXIDEN(4,MODEL)
      CALL MXIDEN(4,VIEW)
      CALL MXIDEN(4,PROJ)
      CALL MPVIEW(0.,DEVP(4),0.,DEVP(5))
      NCPLIN = 0
      NCPLAN = 0
      CALL MXIDEN(4,MVP)
      CALL MXIDEN(4,VP)
      RETURN

      END
      LOGICAL FUNCTION MPLOOK(VX,VY,VZ,PX,PY,PZ,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPLOOK')
      PARAMETER (DPR=57.2958)

      MPLOOK = .FALSE.
      DENTHE = SQRT((PX-VX)**2+ (PZ-VZ)**2)
      IF (DENTHE.EQ.0.) THEN
         SINTHE = 0.
         COSTHE = 1.

      ELSE
         SINTHE = (PX-VX)/DENTHE
         COSTHE = (VZ-PZ)/DENTHE
      END IF

      DENPHI = SQRT((PX-VX)**2+ (PY-VY)**2+ (PZ-VZ)**2)
      IF (DENPHI.EQ.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the same eye position as the reference positio
     *n',2)
         RETURN

      END IF

      MPLOOK = .TRUE.
      SINPHI = (VY-PY)/DENPHI
      COSPHI = DENTHE/DENPHI
      CALL LDTRAN(-VX,-VY,-VZ,TMAT1)
      CALL LDROTA('y',COSTHE,SINTHE,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,TMAT3)
      CALL LDROTA('x',COSPHI,SINPHI,TMAT1)
      CALL MXMULT(4,TMAT3,TMAT1,TMAT2)
      ANG = -TWIST/DPR
      CALL LDROTA('z',COS(ANG),SIN(ANG),TMAT1)
      CALL MXMULT(4,TMAT2,TMAT1,VIEW)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      PEYE(1) = VX
      PEYE(2) = VY
      PEYE(3) = VZ
      PLOOK(1) = PX
      PLOOK(2) = PY
      PLOOK(3) = PZ
      ETWIST = TWIST
      RETURN

      END
      LOGICAL FUNCTION MPORT2(LEFT,RIGHT,BOTTOM,TOP)
      REAL LEFT
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPORT2')
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      MPORT2 = .FALSE.
      IF (RIGHT.EQ.LEFT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the right and left edges of the clipping recta
     *ngle as equal',2)
         RETURN

      END IF

      IF (TOP.EQ.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the top and bottom edges of the clipping recta
     *ngle as equal',2)
         RETURN

      END IF

      MPORT2 = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(1,1) = 2./ (RIGHT-LEFT)
      PROJ(2,2) = 2./ (TOP-BOTTOM)
      PROJ(3,3) = -1.
      PROJ(4,4) = 1.
      PROJ(4,1) = - (RIGHT+LEFT)/ (RIGHT-LEFT)
      PROJ(4,2) = - (TOP+BOTTOM)/ (TOP-BOTTOM)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      LOGICAL FUNCTION MPORT3(LEFT,RIGHT,BOTTOM,TOP,NEAR,FAR)
      REAL LEFT,NEAR
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPORT3')
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      MPORT3 = .FALSE.
      IF (RIGHT.EQ.LEFT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the right and left edges of the clipping recta
     *ngle as equal',2)
         RETURN

      END IF

      IF (TOP.EQ.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the top and bottom edges of the clipping recta
     *ngle as equal',2)
         RETURN

      END IF

      IF (NEAR.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *         'You cannot specify the near clipping plane less than 0.'
     *               ,2)
         RETURN

      END IF

      IF (NEAR.GE.FAR) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the near clipping plane >= to the far clipping
     * plane',2)
         RETURN

      END IF

      MPORT3 = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(1,1) = 2./ (RIGHT-LEFT)
      PROJ(2,2) = 2./ (TOP-BOTTOM)
      PROJ(3,3) = -2./ (FAR-NEAR)
      PROJ(4,4) = 1.
      PROJ(4,1) = - (RIGHT+LEFT)/ (RIGHT-LEFT)
      PROJ(4,2) = - (TOP+BOTTOM)/ (TOP-BOTTOM)
      PROJ(4,3) = - (FAR+NEAR)/ (FAR-NEAR)
      CPNEAR = NEAR
      CPFAR = FAR
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      LOGICAL FUNCTION MPPERS(FOVY,ASPECT,NEAR,FAR)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPPERS')
      PARAMETER (DPR=57.2958)
      REAL NEAR
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      MPPERS = .FALSE.
      IF (NEAR.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *         'You cannot specify the near clipping plane less than 0.'
     *               ,2)
         RETURN

      END IF

      IF (NEAR.GE.FAR) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the near clipping plane >= to the far clipping
     * plane',2)
         RETURN

      END IF

      MPPERS = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(2,2) = 1./TAN((FOVY/DPR)/2.)
      PROJ(1,1) = PROJ(2,2)/ASPECT
      PROJ(3,3) = - (FAR+NEAR)/ (FAR-NEAR)
      PROJ(4,3) = - (2.*FAR*NEAR)/ (FAR-NEAR)
      PROJ(3,4) = -1.
      CPNEAR = NEAR
      CPFAR = FAR
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      SUBROUTINE MPPOLA(DIST,AZIM,INC,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      REAL INC
      PARAMETER (DPR=57.2958)

      ANG = -AZIM/DPR
      CALL LDROTA('y',COS(ANG),SIN(ANG),TMAT1)
      ANG = -INC/DPR
      CALL LDROTA('x',COS(ANG),SIN(ANG),TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,TMAT3)
      ANG = -TWIST/DPR
      CALL LDROTA('z',COS(ANG),SIN(ANG),TMAT1)
      CALL MXMULT(4,TMAT3,TMAT1,TMAT2)
      CALL LDTRAN(0.,0.,-DIST,TMAT1)
      CALL MXMULT(4,TMAT2,TMAT1,VIEW)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RINC = INC/DPR
      RAZIM = AZIM/DPR
      PEYE(1) = DIST*COS(RINC)*SIN(RAZIM)
      PEYE(2) = -DIST*SIN(RINC)
      PEYE(3) = DIST*COS(RINC)*COS(RAZIM)
      PLOOK(1) = 0.
      PLOOK(2) = 0.
      PLOOK(3) = 0.
      ETWIST = TWIST
      RETURN

      END
      SUBROUTINE MPPOPM()
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      DIMENSION EQMAP(193)
      EQUIVALENCE (EQMAP(1),MODEL(1,1))
      COMMON /MPSTCK/SVMAP(195,10),MAPDEP

      DO 2940 I = 1,193
         EQMAP(I) = SVMAP(I,MAPDEP)
 2940 CONTINUE
      NCPLIN = INT(SVMAP(194,MAPDEP))
      NCPLAN = INT(SVMAP(195,MAPDEP))
      IF (MAPDEP.NE.1) THEN
         MAPDEP = MAPDEP - 1
      END IF

      RETURN

      END
      LOGICAL FUNCTION MPPSHM()
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPPSHM')
      DIMENSION EQMAP(193)
      EQUIVALENCE (EQMAP(1),MODEL(1,1))
      COMMON /MPSTCK/SVMAP(195,10),MAPDEP

      MPPSHM = .FALSE.
      IF (MAPDEP.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Stack is full; depth is ten',2)
         RETURN

      END IF

      MPPSHM = .TRUE.
      MAPDEP = MAPDEP + 1
      DO 2960 I = 1,193
         SVMAP(I,MAPDEP) = EQMAP(I)
 2960 CONTINUE
      SVMAP(194,MAPDEP) = FLOAT(NCPLIN)
      SVMAP(195,MAPDEP) = FLOAT(NCPLAN)
      RETURN

      END
      SUBROUTINE MPPUTM(TMODEL,TVIEW,TPROJ,TVWPT)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      DIMENSION TMODEL(4,4),TVIEW(4,4),TPROJ(4,4),TVWPT(4)

      DO 2980 I = 1,4
         DO 3000 J = 1,4
            MODEL(J,I) = TMODEL(J,I)
            VIEW(J,I) = TVIEW(J,I)
            PROJ(J,I) = TPROJ(J,I)
 3000    CONTINUE
         VWPORT(I) = TVWPT(I)
 2980 CONTINUE
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      SUBROUTINE MPRLOC(EYE,LOOKAT,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      DIMENSION EYE(3),LOOKAT(3)
      REAL LOOKAT

      DO 3020 I = 1,3
         EYE(I) = PEYE(I)
         LOOKAT(I) = PLOOK(I)
 3020 CONTINUE
      TWIST = ETWIST
      RETURN

      END
      SUBROUTINE MPRESE()
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      CALL MXIDEN(4,MODEL)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      LOGICAL FUNCTION MPROTA(ANGLE,AXIS)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPROTA')
      CHARACTER*1 AXIS,TAXIS
      PARAMETER (DPR=57.2958)

      MPROTA = .FALSE.
      CALL CHRUP(AXIS,TAXIS)
      IF (TAXIS.NE.'X' .AND. TAXIS.NE.'Y' .AND. TAXIS.NE.'Z') THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Illegal axis '//AXIS//' specified',2)
         RETURN

      END IF

      MPROTA = .TRUE.
      CALL LDROTA(AXIS,COS(ANGLE/DPR),SIN(ANGLE/DPR),TMAT1)
      CALL MXCOPY(4,MODEL,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,MODEL)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      SUBROUTINE MPSCAL(X,Y,Z)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      CALL LDSCAL(X,Y,Z,TMAT1)
      CALL MXCOPY(4,MODEL,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,MODEL)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      SUBROUTINE MPTRAN(X,Y,Z)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      CALL LDTRAN(X,Y,Z,TMAT1)
      CALL MXCOPY(4,MODEL,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,MODEL)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      LOGICAL FUNCTION MPVIEW(LEFT,RIGHT,BOTTOM,TOP)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
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
      REAL LEFT
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPVIEW')

      MPVIEW = .FALSE.
      IF (LEFT.GE.RIGHT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *  'You cannot specify the left viewport edge >= to the right edge'
     *               ,2)
         RETURN

      END IF

      IF (TOP.LE.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *  'You cannot specify the top viewport edge <= to the bottom edge'
     *               ,2)
         RETURN

      END IF

      IF (TOP.LT.0 .OR. TOP.GT.DEVP(5)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Top viewport specification out of range',2)
         RETURN

      END IF

      IF (BOTTOM.LT.0 .OR. BOTTOM.GT.DEVP(5)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Bottom viewport specification out of range'
     *               ,2)
         RETURN

      END IF

      IF (LEFT.LT.0 .OR. LEFT.GT.DEVP(4)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Left viewport specification out of range',
     *               2)
         RETURN

      END IF

      IF (RIGHT.LT.0 .OR. RIGHT.GT.DEVP(4)) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Right viewport specification out of range',
     *               2)
         RETURN

      END IF

      MPVIEW = .TRUE.
      VWPORT(1) = LEFT
      VWPORT(2) = RIGHT
      VWPORT(3) = BOTTOM
      VWPORT(4) = TOP
      RETURN

      END
      LOGICAL FUNCTION MPWIND(LEFT,RIGHT,BOTTOM,TOP,NEAR,FAR)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      REAL LEFT,NEAR
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPWIND')

      MPWIND = .FALSE.
      IF (NEAR.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *         'You cannot specify the near clipping plane less than 0.'
     *               ,2)
         RETURN

      END IF

      IF (NEAR.GE.FAR) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the near clipping plane >= to the far clipping
     * plane',2)
         RETURN

      END IF

      IF (RIGHT.EQ.LEFT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the right and left edges of the viewing frustr
     *um as equal',2)
         RETURN

      END IF

      IF (TOP.EQ.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the top and bottom edges of the viewing frustr
     *um as equal',2)
         RETURN

      END IF

      MPWIND = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(1,1) = 2.*NEAR/ (RIGHT-LEFT)
      PROJ(2,2) = 2.*NEAR/ (TOP-BOTTOM)
      PROJ(3,1) = (RIGHT+LEFT)/ (RIGHT-LEFT)
      PROJ(3,2) = (TOP+BOTTOM)/ (TOP-BOTTOM)
      PROJ(3,3) = - (FAR+NEAR)/ (FAR-NEAR)
      PROJ(3,4) = -1.
      PROJ(4,3) = - (2.*FAR*NEAR)/ (FAR-NEAR)
      CPNEAR = NEAR
      CPFAR = FAR
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
      SUBROUTINE LDROTA(AXIS,COSANG,SINANG,MAT)
      REAL MAT(4,4)
      CHARACTER*(*) AXIS
      CHARACTER*1 TAXIS

      CALL MXIDEN(4,MAT)
      TAXIS = AXIS
      CALL CHRUP(TAXIS,TAXIS)
      IF (TAXIS.EQ.'X') THEN
         MAT(2,2) = COSANG
         MAT(2,3) = SINANG
         MAT(3,2) = -SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Y') THEN
         MAT(1,1) = COSANG
         MAT(1,3) = -SINANG
         MAT(3,1) = SINANG
         MAT(3,3) = COSANG

      ELSE IF (TAXIS.EQ.'Z') THEN
         MAT(1,1) = COSANG
         MAT(1,2) = SINANG
         MAT(2,1) = -SINANG
         MAT(2,2) = COSANG
      END IF

      RETURN

      END
      SUBROUTINE LDSCAL(X,Y,Z,MAT)
      REAL MAT(4,4)

      CALL MXIDEN(4,MAT)
      MAT(1,1) = X
      MAT(2,2) = Y
      MAT(3,3) = Z
      RETURN

      END
      SUBROUTINE LDTRAN(X,Y,Z,MAT)
      REAL MAT(4,4)

      CALL MXIDEN(4,MAT)
      MAT(4,1) = X
      MAT(4,2) = Y
      MAT(4,3) = Z
      RETURN

      END
      SUBROUTINE MPMUL2(N,X0,Y0,MAT,RES1,RES2,RES3,RES4)
      DIMENSION X0(*),Y0(*),MAT(4,4),RES1(*),RES2(*),RES3(*),RES4(*)
      DIMENSION VEC(4),VECR(4)
      REAL MAT
      DATA VEC(3)/0./,VEC(4)/1./

      DO 3040 I = 1,N
         VEC(1) = X0(I)
         VEC(2) = Y0(I)
         DO 3060 J = 1,4
            VECR(J) = 0.0
            DO 3080 K = 1,4
               VECR(J) = VECR(J) + MAT(K,J)*VEC(K)
 3080       CONTINUE
 3060    CONTINUE
         RES1(I) = VECR(1)
         RES2(I) = VECR(2)
         RES3(I) = VECR(3)
         RES4(I) = VECR(4)
 3040 CONTINUE
      RETURN

      END
      SUBROUTINE MPMUL3(N,X0,Y0,Z0,MAT,RES1,RES2,RES3,RES4)
      DIMENSION X0(*),Y0(*),Z0(*),MAT(4,4),RES1(*),RES2(*),RES3(*),
     *          RES4(*),VEC(4),VECR(4)
      REAL MAT
      DATA VEC(4)/1./

      DO 3100 I = 1,N
         VEC(1) = X0(I)
         VEC(2) = Y0(I)
         VEC(3) = Z0(I)
         DO 3120 J = 1,4
            VECR(J) = 0.0
            DO 3140 K = 1,4
               VECR(J) = VECR(J) + MAT(K,J)*VEC(K)
 3140       CONTINUE
 3120    CONTINUE
         RES1(I) = VECR(1)
         RES2(I) = VECR(2)
         RES3(I) = VECR(3)
         RES4(I) = VECR(4)
 3100 CONTINUE
      RETURN

      END
      SUBROUTINE MPMUL4(N,MASK,ARR1,ARR2,ARR3,ARR4,MAT,RES1,RES2,RES3,
     *                  RES4)
      DIMENSION VEC(4),VECR(4),ARR1(*),ARR2(*),ARR3(*),ARR4(*),RES1(*),
     *          RES2(*),RES3(*),RES4(*),MAT(4,4)
      REAL MAT
      INTEGER IZBIT(32)
      DATA IZBIT/1,          2,          4,          8,
     *          16,         32,         64,        128,
     *         256,        512,       1024,       2048,
     *        4096,       8192,      16384,      32768,
     *       65536,     131072,     262144,     524288,
     *     1048576,    2097152,    4194304,    8388608,
     *    16777216,   33554432,   67108864,  134217728,
     *   268435456,  536870912, 1073741824, 2147483648/

      IF (MASK.EQ.0) THEN
         RETURN

      END IF

      DO 3160 I = 1,N
         IF (IAND(MASK,IZBIT(I)).NE.0) THEN
            VEC(1) = ARR1(I)
            VEC(2) = ARR2(I)
            VEC(3) = ARR3(I)
            VEC(4) = ARR4(I)
            DO 3180 J = 1,4
               VECR(J) = 0.0
               DO 3200 K = 1,4
                  VECR(J) = VECR(J) + MAT(K,J)*VEC(K)
 3200          CONTINUE
 3180       CONTINUE
            RES1(I) = VECR(1)
            RES2(I) = VECR(2)
            RES3(I) = VECR(3)
            RES4(I) = VECR(4)
         END IF

 3160 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION CHRCI(LINE,INTE)
      CHARACTER LINE* (*),FORM*10,CL*2

      CALL CHRTRM(LINE,LL)
      CALL CHRIC(LL,CL,NL)
      FORM = '(i'//CL(1:NL)//')'
      READ (LINE(1:LL),FORM,ERR=10) INTE
      CHRCI = .TRUE.
      RETURN

   10 CHRCI = .FALSE.
      RETURN

      END
      LOGICAL FUNCTION CHRCMP(KWD,PART1,PART2)
      CHARACTER*(*) KWD
      CHARACTER*(*) PART1
      CHARACTER*(*) PART2

      CALL CHRTRM(KWD,LK)
      CALL CHRTRM(PART1,LF)
      CALL CHRTRM(PART2,LV)
      IF (LK.LT.LF .OR. LK.GT.LF+LV) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (KWD(1:LF).NE.PART1(1:LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      IF (LK.EQ.LF) THEN
         CHRCMP = .TRUE.
         RETURN

      END IF

      IF (KWD(LF+1:LK).NE.PART2(1:LK-LF)) THEN
         CHRCMP = .FALSE.
         RETURN

      END IF

      CHRCMP = .TRUE.
      RETURN

      END
      LOGICAL FUNCTION CHRCR(LINE,REAL)
      CHARACTER LINE* (*),FORM*10,CL*2

      CALL CHRTRM(LINE,LL)
      CALL CHRIC(LL,CL,NL)
      FORM = '(f'//CL(1:NL)//'.0)'
      READ (LINE(1:LL),FORM,ERR=10) REAL
      CHRCR = .TRUE.
      RETURN

   10 CHRCR = .FALSE.
      RETURN

      END
      SUBROUTINE CHRDN(LINE1,LINE2)
      CHARACTER*(*) LINE1,LINE2

      L = MIN(LEN(LINE1),LEN(LINE2))
      I = 1
 2000 IF (.NOT. (I.LE.L)) GO TO 2020
      ICH = ICHAR(LINE1(I:I))
      IF (ICH.GE.65 .AND. ICH.LE.90) THEN
         ICH = ICH + 32
      END IF

      LINE2(I:I) = CHAR(ICH)
 2010 I = I + 1
      GO TO 2000

 2020 CONTINUE
      RETURN

      END
      SUBROUTINE CHRIC(JIN,DS,ND)
      CHARACTER*(*) DS
      CHARACTER SDG*16

      DS = ' '
      J = ABS(JIN)
      ND = 0
 2030 CONTINUE
      ND = ND + 1
      KD = MOD(J,10)
      J = J/10
      SDG(ND:ND) = CHAR(ICHAR('0')+KD)
 2040 IF (.NOT. (J.EQ.0)) GO TO 2030
      IJ = 1
      IF (JIN.LT.0) THEN
         IJ = 2
         DS(1:1) = '-'
         ND = ND + 1
      END IF

      I = IJ
 2060 IF (.NOT. (I.LE.ND)) GO TO 2080
      DS(I:I) = SDG(ND-I+1:ND-I+1)
 2070 I = I + 1
      GO TO 2060

 2080 CONTINUE
      RETURN

      END
      INTEGER FUNCTION CHRLEN(S)
      CHARACTER*(*) S

      L = LEN(S)
      J = INDEX(S,CHAR(0))
      IF (J.EQ.0) THEN
         I = L

      ELSE
         I = J - 1
      END IF

 2090 IF (.NOT. (I.GT.0.AND.S(I:I).EQ.' ')) GO TO 2100
      I = I - 1
      GO TO 2090

 2100 CONTINUE
      CHRLEN = I
      RETURN

      END
      SUBROUTINE CHRRVC(BUFF,TXT,L)
      CHARACTER*16 LOCTXT
      CHARACTER*(*) TXT

      LT = LEN(TXT)
      WRITE (LOCTXT,'(1pg13.6)') BUFF
      DO 2110 I = 1,LT
         IF (LOCTXT(I:I).NE.' ') THEN
            L1 = I
            GO TO 2120

         END IF

 2110 CONTINUE
 2120 CONTINUE
      DO 2130 I = L1,LT
         IF (LOCTXT(I:I).EQ.' ') THEN
            L2 = I
            GO TO 2140

         END IF

 2130 CONTINUE
 2140 CONTINUE
      TXT = LOCTXT(L1:L2)
      L = L2 - L1
      RETURN

      END
      SUBROUTINE CHRSTR(LINE1,LINE2,L)
      CHARACTER*(*) LINE1,LINE2
      CHARACTER CH

      K = LEN(LINE1)
      J = 1
 2150 IF (.NOT. (J.LT.K)) GO TO 2170
      CH = LINE1(J:J)
      IF (CH.EQ.' ') THEN
         GO TO 2160

      END IF

      IF (CH.EQ.CHAR(9)) THEN
         GO TO 2160

      END IF

      GO TO 2170

 2160 J = J + 1
      GO TO 2150

 2170 CONTINUE
      IF (J.GT.K) THEN
         L = 0
         LINE2 = ' '
         RETURN

      END IF

      LINE2(1:K-J+1) = LINE1(J:K)
      CALL CHRTRM(LINE2(1:K-J+1),L)
      RETURN

      END
      SUBROUTINE CHRTRM(LINE,L)
      CHARACTER*(*) LINE
      CHARACTER CH

      L1 = LEN(LINE)
 2180 IF (.NOT. (L1.GT.0)) GO TO 2200
      CH = LINE(L1:L1)
      IF (CH.EQ.' ') THEN
         GO TO 2190

      END IF

      IF (CH.EQ.CHAR(9)) THEN
         GO TO 2190

      END IF

      IF (CH.EQ.CHAR(0)) THEN
         GO TO 2190

      END IF

      GO TO 2200

 2190 L1 = L1 - 1
      GO TO 2180

 2200 CONTINUE
      L = L1
      RETURN

      END
      SUBROUTINE CHRUP(LINE1,LINE2)
      CHARACTER*(*) LINE1,LINE2

      L = MIN(LEN(LINE1),LEN(LINE2))
      I = 1
 2210 IF (.NOT. (I.LE.L)) GO TO 2230
      ICH = ICHAR(LINE1(I:I))
      IF (ICH.GE.97 .AND. ICH.LE.122) THEN
         ICH = ICH - 32
      END IF

      LINE2(I:I) = CHAR(ICH)
 2220 I = I + 1
      GO TO 2210

 2230 CONTINUE
      RETURN

      END
      SUBROUTINE CPUDAC(DS)
      CHARACTER*(*) DS

      DS = ' '
      RETURN

      END
      SUBROUTINE CPUDAT(IY,MO,ID,IH,IM,IS)

      RETURN

      END
      SUBROUTINE CPUERR(STR,DISP)
      CHARACTER*(*) STR

      CALL SIORPT('CPU',STR,DISP)
      RETURN

      END
      SUBROUTINE CPUCML(LINE,PROMPT,L)
      CHARACTER*(*) LINE,PROMPT
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      IF (FIRST) THEN
         FIRST = .FALSE.
         LINE = ' '
         RETURN

      ELSE
         IF (L.LT.0) THEN
            CALL CHRTRM(PROMPT,LP)

         ELSE
            LP = L
         END IF

         IF (LP.NE.0) THEN
            WRITE (6,'(1x,a)') PROMPT(1:LP)
         END IF

         READ (5,'(a)',ERR=10,END=10) LINE
         CALL CHRTRM(LINE,L)
      END IF

      RETURN

   10 L = -1
      RETURN

      END
      SUBROUTINE CPUMVU(A,B,L)
      IMPLICIT INTEGER (A-Z)
      DIMENSION A(*),B(*)

      DO 2240 I = 1,L
         B(I) = A(I)
 2240 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION CPUNAL(IU)
      CHARACTER*40 CPUHLB
      COMMON /CPHLBN/CPUHLB
      INTEGER CPUNIT(20)
      INTEGER CPUNIF
      COMMON /CPUN/CPUNIT,CPUNIF
      CHARACTER*20 TERM
      CHARACTER*80 LINBUF(24)
      COMMON /REBUF/LINBUF,TERM
      INTEGER LINLEN(24)
      COMMON /REBUF2/NL,LINLEN

      CPUNAL = .TRUE.
      IF (CPUNIF.NE.12345) THEN
         I = 1
 2260    IF (.NOT. (I.LE.20)) GO TO 2280
         CPUNIT(I) = 80 + I - 1
 2270    I = I + 1
         GO TO 2260

 2280    CONTINUE
         CPUNIF = 12345
      END IF

      I = 1
 2290 IF (.NOT. (I.LE.20)) GO TO 2310
      IF (CPUNIT(I).EQ.0) THEN
         GO TO 2300

      END IF

      IU = CPUNIT(I)
      CPUNIT(I) = 0
      RETURN

 2300 I = I + 1
      GO TO 2290

 2310 CONTINUE
      CPUNAL = .FALSE.
      CALL CPUERR('Cannot allocate logical unit.',2)
      RETURN

      END
      SUBROUTINE CPUNDE(IU)
      CHARACTER*40 CPUHLB
      COMMON /CPHLBN/CPUHLB
      INTEGER CPUNIT(20)
      INTEGER CPUNIF
      COMMON /CPUN/CPUNIT,CPUNIF
      CHARACTER*20 TERM
      CHARACTER*80 LINBUF(24)
      COMMON /REBUF/LINBUF,TERM
      INTEGER LINLEN(24)
      COMMON /REBUF2/NL,LINLEN

      J = IU - 80 + 1
      CPUNIT(J) = IU
      RETURN

      END
      SUBROUTINE CPUTBK(V)
      LOGICAL V

      RETURN

      END
      SUBROUTINE CPUQA(DATE1,TIME,USER,JOBID,DIV,ENUM,CASENO,CLASS)
      CHARACTER*(*) DATE1
      CHARACTER*(*) TIME
      CHARACTER*(*) USER
      CHARACTER*(*) JOBID
      CHARACTER*(*) DIV
      CHARACTER*(*) ENUM
      CHARACTER*(*) CASENO
      CHARACTER*(*) CLASS

      DATE1 = ' '
      TIME = ' '
      JOBID = ' '
      USER = ' '
      DIV = ' '
      ENUM = ' '
      CASENO = ' '
      CLASS = ' '
      RETURN

      END
      SUBROUTINE CPUWAT(ISEC)
      INTEGER ISEC

      RETURN

      END
      SUBROUTINE LXCLN
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504

      ELSE
         JLINE = MIN(JLINE+K,504)
      END IF

      RETURN

      END
      LOGICAL FUNCTION LXRNL(VAL,N,CH)
      REAL*8 VAL(*)
      CHARACTER CH
      LOGICAL LDUM,LXGTCH,LXGTWH,LXREAL
      REAL*8 XX,XY

      I = 1
 2320 IF (.NOT. (I.LE.N)) GO TO 2340
      VAL(I) = 0.
 2330 I = I + 1
      GO TO 2320

 2340 CONTINUE
      N = 0
      LXRNL = .TRUE.
 2350 CONTINUE
      LDUM = LXGTWH(CH)
      IF (LXREAL(VAL(N+1),CH)) THEN
         N = N + 1
         LDUM = LXGTWH(CH)
         IF (LXGTCH('#',CH)) THEN
            RETURN

         END IF

         IF (LXGTCH(',',CH)) THEN
            GO TO 2360

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            RETURN

         END IF

         IF (LXGTCH('*',CH)) THEN
            LDUM = LXGTWH(CH)
            XX = VAL(N)
            IF (.NOT.LXREAL(XY,CH)) THEN
               LXRNL = .FALSE.
               RETURN

            END IF

            M = XX + .1
            N0 = N
 2380       IF (.NOT. (N.LT.M+N0)) GO TO 2400
            VAL(N) = XY
 2390       N = N + 1
            GO TO 2380

 2400       CONTINUE
            LDUM = LXGTWH(CH)
            N = N0 + MAX(M-1,0)
            IF (LXGTCH(',',CH)) THEN
               GO TO 2360

            END IF

            IF (LXGTCH('#',CH)) THEN
               RETURN

            END IF

            IF (CH.EQ.CHAR(0)) THEN
               RETURN

            END IF

         END IF

      ELSE IF (LXGTCH(',',CH)) THEN
         VAL(N+1) = 0.
         N = N + 1

      ELSE IF (LXGTCH('#',CH)) THEN
         RETURN

      ELSE IF (CH.EQ.CHAR(0)) THEN
         IF (N.EQ.0) THEN
            RETURN

         END IF

         N = N + 1
         VAL(N) = 0.
         RETURN

      ELSE
         LXRNL = .FALSE.
         RETURN

      END IF

 2360 GO TO 2350

      END
      SUBROUTINE LXSETP(LINE)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER*255 LOCLIN
      CHARACTER*80 TMPLIN
      INTEGER CHRLEN

      IF (LXINIT.NE.12345) THEN
         CALL LXRST
      END IF

      LOCLIN = LINE
      L = CHRLEN(LOCLIN)
      K = JLINE - L - 1
      IF (K.LE.0) THEN
         TMPLIN = 'Buffer overflow in lxsetp: '//LINE(1:L)
         CALL LXERR(TMPLIN,3)
         RETURN

      END IF

      IF (L.GT.0) THEN
         ILINE(K:JLINE-1) = LINE(1:L)//CHAR(0)

      ELSE
         ILINE(K:JLINE-1) = CHAR(0)
      END IF

      JLINE = K
      RETURN

      END
      SUBROUTINE LXRST
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      LXINIT = 12345
      JLINE = 504
      ILINE(JLINE:JLINE) = CHAR(0)
      RETURN

      END
      SUBROUTINE LXSTK(LINE)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE

      L = LEN(LINE)
      JLINE = JLINE - L
      ILINE(JLINE:JLINE+L-1) = LINE(1:L)
      RETURN

      END
      SUBROUTINE LXREM(LINE,L)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER*80 TMPLIN

      K = INDEX(ILINE(JLINE:),CHAR(0))
      IF (K.EQ.0) THEN
         JLINE = 504
         L = 0
         LINE = ' '

      ELSE
         L = K - 1
         IF (L.GT.LEN(LINE)) THEN
            L = LEN(LINE)
            LINE = ILINE(JLINE:JLINE+L-1)
            TMPLIN = 'Remainder truncated:'//LINE
            CALL LXERR(TMPLIN,1)

         ELSE IF (L.GT.0) THEN
            LINE = ILINE(JLINE:JLINE+K-1)
         END IF

         JLINE = JLINE + K
      END IF

      RETURN

      END
      SUBROUTINE LXERR(MSG,DISP)
      CHARACTER*(*) MSG
      INTEGER DISP
      CHARACTER*80 LOCMSG

      LOCMSG = MSG
      IF (DISP.GE.2) THEN
         CALL LXCLN
      END IF

      CALL SIORPT('LEX',LOCMSG,DISP)
      RETURN

      END
      INTEGER FUNCTION LXSV()
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      LXSV = JLINE
      RETURN

      END
      SUBROUTINE LXRS(IP)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT

      IF (IP.GT.0 .AND. IP.LE.504) THEN
         JLINE = IP

      ELSE
         CALL LXERR('Illegal pointer restoration',3)
      END IF

      RETURN

      END
      LOGICAL FUNCTION LXANY(CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER CH

      CH = ILINE(JLINE:JLINE)
      LXANY = CH .NE. CHAR(0)
      IF (LXANY) THEN
         JLINE = JLINE + 1
      END IF

      RETURN

      END
      LOGICAL FUNCTION LXSCAN(DELIM,REM,L,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) REM
      INTEGER L
      CHARACTER CH

      REM = ' '
      L = 0
      LXSCAN = .FALSE.
 2410 CONTINUE
      CH = ILINE(JLINE:JLINE)
      IF (CH.EQ.CHAR(0)) THEN
         RETURN

      END IF

      JLINE = JLINE + 1
      LXSCAN = (INDEX(DELIM,CH).NE.0)
      IF (LXSCAN) THEN
         RETURN

      END IF

      L = L + 1
      IF (L.GT.LEN(REM)) THEN
         GO TO 2420

      END IF

      REM(L:L) = CH
 2420 GO TO 2410

      END
      LOGICAL FUNCTION LXGTCH(CH1,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) CH1,CH

      L = LEN(CH1)
      IF (CH1.EQ.ILINE(JLINE:JLINE+L-1)) THEN
         JLINE = JLINE + L
         CH = ILINE(JLINE:JLINE)
         LXGTCH = .TRUE.
         RETURN

      END IF

      CH = ILINE(JLINE:JLINE)
      LXGTCH = .FALSE.
      RETURN

      END
      LOGICAL FUNCTION LXGTBP(STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) STR
      CHARACTER CH
      LOGICAL QFLAG
      CHARACTER QCH

      NS = 0
      CH = ILINE(JLINE:JLINE)
      LXGTBP = (CH.EQ.'(') .OR. (CH.EQ.'[') .OR. (CH.EQ.'{')
      IF (.NOT.LXGTBP) THEN
         RETURN

      END IF

      J = JLINE
      LP = 1
      QFLAG = .FALSE.
 2440 CONTINUE
      J = J + 1
      CH = ILINE(J:J)
      IF (.NOT.QFLAG) THEN
         IF (CH.EQ.'(' .OR. CH.EQ.'[' .OR. CH.EQ.'{') THEN
            LP = LP + 1

         ELSE IF (CH.EQ.')' .OR. CH.EQ.']' .OR. CH.EQ.'}') THEN
            LP = LP - 1

         ELSE IF (CH.EQ.'''' .OR. CH.EQ.'"' .OR. CH.EQ.CHAR(96)) THEN
            QFLAG = .TRUE.
            QCH = CH

         ELSE IF (CH.EQ.CHAR(0)) THEN
            LXGTBP = .FALSE.
            RETURN

         END IF

         IF (LP.EQ.0) THEN
            GO TO 2460

         END IF

         NS = NS + 1
         STR(NS:NS) = CH

      ELSE
         NS = NS + 1
         STR(NS:NS) = CH
         IF (CH.EQ.QCH) THEN
            IF (CH.EQ.ILINE(J+1:J+1)) THEN
               NS = NS + 1
               STR(NS:NS) = CH
               J = J + 1
               GO TO 2450

            ELSE
               QFLAG = .FALSE.
            END IF

         ELSE IF (CH.EQ.CHAR(0)) THEN
            LXGTBP = .FALSE.
            RETURN

         END IF

      END IF

 2450 GO TO 2440

 2460 CONTINUE
      JLINE = J + 1
      LXGTBP = .TRUE.
      RETURN

      END
      LOGICAL FUNCTION LXGTQT(STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) STR
      CHARACTER CH

      CH = ILINE(JLINE:JLINE)
      LXGTQT = (CH.EQ.CHAR(39)) .OR. (CH.EQ.CHAR(34))
      IF (.NOT.LXGTQT) THEN
         RETURN

      END IF

      NS = 0
      J = JLINE
      LP = 1
 2470 IF (.NOT. (LP.GT.0)) GO TO 2480
      J = J + 1
      CH = ILINE(J:J)
      IF (CH.EQ.CHAR(39)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(34)) THEN
         LP = LP - 1

      ELSE IF (CH.EQ.CHAR(0)) THEN
         LXGTQT = .FALSE.
         RETURN

      END IF

      NS = NS + 1
      STR(NS:NS) = CH
      GO TO 2470

 2480 CONTINUE
      NS = NS - 1
      JLINE = J + 1
      RETURN

      END
      LOGICAL FUNCTION LXSET(SET,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) SET,CH

      CH = ILINE(JLINE:JLINE)
      LXSET = INDEX(SET,CH) .NE. 0 .AND. CH .NE. CHAR(0)
      IF (LXSET) THEN
         JLINE = JLINE + 1
      END IF

      RETURN

      END
      LOGICAL FUNCTION LXGTWH(CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER CH

      J0 = JLINE
      CH = ILINE(JLINE:JLINE)
 2490 IF (.NOT. (CH.EQ.' '.OR.CH.EQ.CHAR(9))) GO TO 2510
      JLINE = JLINE + 1
 2500 CH = ILINE(JLINE:JLINE)
      GO TO 2490

 2510 CONTINUE
      LXGTWH = JLINE .NE. J0
      RETURN

      END
      LOGICAL FUNCTION LXNBS(LINE,L)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) LINE
      CHARACTER CH
      LOGICAL LDUM
      LOGICAL LXGTWH

      LDUM = LXGTWH(CH)
      CH = ILINE(JLINE:JLINE)
      J0 = JLINE
 2520 IF (.NOT. (INDEX(' '//CHAR(9)//CHAR(0),CH).EQ.0)) GO TO 2530
      JLINE = JLINE + 1
      CH = ILINE(JLINE:JLINE)
      GO TO 2520

 2530 CONTINUE
      L = JLINE - J0
      IF (J0.EQ.JLINE) THEN
         LXNBS = .FALSE.
         RETURN

      END IF

      LINE = ILINE(J0:JLINE-1)
      LXNBS = .TRUE.
      RETURN

      END
      LOGICAL FUNCTION LXSYMB(SYM,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) SYM,CH
      CHARACTER*26 ALPHA,BETA
      DATA ALPHA/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA BETA/'abcdefghijklmnopqrstuvwxyz'/

      SYM = ' '
      LS = LEN(SYM)
      CH = ILINE(JLINE:JLINE)
      IF (INDEX(ALPHA//BETA,CH).EQ.0) THEN
         LXSYMB = .FALSE.
         RETURN

      END IF

      LXSYMB = .TRUE.
      J0 = JLINE
 2540 IF (.NOT. (INDEX(ALPHA//BETA//'0123456789_$',CH).NE.0)) GO TO 2550
      JLINE = JLINE + 1
      CH = ILINE(JLINE:JLINE)
      GO TO 2540

 2550 CONTINUE
      SYM = ILINE(J0:JLINE-1)
      NS = MIN(LS,JLINE-J0)
      RETURN

      END
      LOGICAL FUNCTION LXSYM2(SYM,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) SYM,CH
      CHARACTER*65 S

      S =
     *'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_$.
     *'
      SYM = ' '
      CH = ILINE(JLINE:JLINE)
      KS = 0
 2560 IF (.NOT. (INDEX(S,CH).NE.0)) GO TO 2570
      CH = ILINE(JLINE:JLINE)
      JLINE = JLINE + 1
      KS = KS + 1
      SYM(KS:KS) = CH
      GO TO 2560

 2570 CONTINUE
      NS = KS
      LXSYM2 = KS .GT. 0
      RETURN

      END
      LOGICAL FUNCTION LXNUMB(N,ND,CH)
      CHARACTER*(*) CH
      REAL*8 N
      LOGICAL LXSET

      N = 0.
      ND = 0
 2580 IF (.NOT. (LXSET('0123456789',CH))) GO TO 2590
      N = N*10. + FLOAT(ICHAR(CH)-ICHAR('0'))
      ND = ND + 1
      GO TO 2580

 2590 CONTINUE
      IF (ND.EQ.0) THEN
         LXNUMB = .FALSE.
         RETURN

      END IF

      LXNUMB = .TRUE.
      RETURN

      END
      LOGICAL FUNCTION LXREAL(VALUE,CH)
      LOGICAL LXSET,LXNUMB
      LOGICAL PREF,POSF
      REAL*8 VALUE,PREDIG,ESIGN
      REAL*8 FN
      CHARACTER CH
      INTEGER ND

      LXREAL = .FALSE.
      ESIGN = 1.
      PREDIG = 0.
      ISIGN = 1
      ISAVE = LXSV()
      IF (LXSET('+-',CH)) THEN
         IF (CH.EQ.'-') THEN
            ISIGN = -1
         END IF

      END IF

      PREF = LXNUMB(PREDIG,ND,CH)
      POSF = .FALSE.
      IF (LXSET('.',CH)) THEN
         IF (LXNUMB(FN,ND,CH)) THEN
            POSF = .TRUE.
            PREDIG = PREDIG + FN*10.** (FLOAT(-ND))
         END IF

      END IF

      PREDIG = PREDIG*ISIGN
      IF (.NOT. (PREF.OR.POSF)) THEN
         CALL LXRS(ISAVE)
         RETURN

      END IF

      IF (LXSET('EeDdQq',CH)) THEN
         IF (LXSET('+-',CH)) THEN
            IF (CH.EQ.'-') THEN
               ESIGN = -1.
            END IF

         END IF

         IF (LXNUMB(FN,ND,CH)) THEN
            PREDIG = PREDIG*10.** (ESIGN*FN)

         ELSE
            CALL LXRS(ISAVE)
            RETURN

         END IF

      END IF

      VALUE = PREDIG
      LXREAL = .TRUE.
      RETURN

      END
c$$$      SUBROUTINE LXSWC
c$$$      CHARACTER CH
c$$$      LOGICAL LDUM
c$$$      LOGICAL LXGTCH
c$$$      LOGICAL LXGTWH
c$$$
c$$$      LDUM = LXGTWH(CH)
c$$$      LDUM = LXGTCH(',',CH)
c$$$      LDUM = LXGTWH(CH)
c$$$      RETURN
c$$$
c$$$      END
      LOGICAL FUNCTION LXSCNP(DELIM,STR,NS,CH)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*504 ILINE
      COMMON /LXCOM1/ILINE
      COMMON /LXCOM2/JLINE,LXINIT
      CHARACTER*(*) DELIM
      CHARACTER*(*) STR
      INTEGER NS
      CHARACTER CH
      LOGICAL QFLAG
      CHARACTER QCH

      NS = 0
      J = JLINE
      LP = 0
      QFLAG = .FALSE.
 2600 CONTINUE
      CH = ILINE(J:J)
      IF (LP.EQ.0) THEN
         ID = INDEX(DELIM,CH)
         IF (ID.GT.0) THEN
            GO TO 2620

         END IF

      END IF

      IF (.NOT.QFLAG) THEN
         IF (LP.EQ.0) THEN
            ID = INDEX(DELIM,CH)
            IF (ID.GT.0) THEN
               GO TO 2620

            END IF

         END IF

         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.'(' .OR. CH.EQ.'[' .OR. CH.EQ.'{') THEN
            LP = LP + 1

         ELSE IF (CH.EQ.')' .OR. CH.EQ.']' .OR. CH.EQ.'}') THEN
            LP = LP - 1

         ELSE IF (CH.EQ.'''' .OR. CH.EQ.'"' .OR. CH.EQ.CHAR(96)) THEN
            QFLAG = .TRUE.
            QCH = CH
         END IF

         IF (LP.LT.0) THEN
            GO TO 2620

         END IF

         NS = NS + 1
         STR(NS:NS) = CH

      ELSE
         NS = NS + 1
         STR(NS:NS) = CH
         IF (CH.EQ.CHAR(0)) THEN
            GO TO 2620

         ELSE IF (CH.EQ.QCH) THEN
            IF (CH.EQ.ILINE(J+1:J+1)) THEN
               NS = NS + 1
               STR(NS:NS) = CH
               J = J + 2
               GO TO 2610

            ELSE
               QFLAG = .FALSE.
            END IF

         END IF

      END IF

      J = J + 1
 2610 GO TO 2600

 2620 CONTINUE
      JLINE = J
      LXSCNP = ((LP.EQ.0) .AND. .NOT.QFLAG)
      RETURN

      END
      SUBROUTINE MEMINI(MEMRY,LENGTH)
      INTEGER MEMRY(*)
      INTEGER LENGTH

      DO 2630 I = 1,LENGTH
         MEMRY(1) = LENGTH
 2630 CONTINUE
      RETURN

      END
      INTEGER FUNCTION MEMALL(LENGTH,MEMRY)
      INTEGER LENGTH
      INTEGER MEMRY(*)
      INTEGER LR
      INTEGER LMX
      INTEGER IP
      INTEGER IPN
      INTEGER LAV
      INTEGER LT
      INTEGER IB

      IF (LENGTH.LE.0) THEN
         WRITE (6,*) ' Cannot allocate a segment of length zero.'
         CALL CPUTBK(.FALSE.)
      END IF

      LR = LENGTH + MOD(LENGTH,2)
      LMX = MEMRY(1)
      IB = 0
      IP = 3
 2650 IF (.NOT. (IB.EQ.0)) GO TO 2660
      IPN = MEMRY(IP)
      IF (IPN.LE.0) THEN
         LT = IP + LR + 1
         IF (LT.GE.LMX) THEN
            WRITE (6,*) ' Cannot allocate space in memall.'
            WRITE (6,*) ' lt >= lmx '
            WRITE (6,*) ' lt,lr,lmx,length: ',LT,LR,LMX,LENGTH
            MEMALL = (-1)
            RETURN

         END IF

         IF (IPN.EQ.0) THEN
            LAV = LMX - IP - 1

         ELSE IF (IPN.LT.0) THEN
            LAV = -IPN - IP - 2
         END IF

         IF (LAV.GT.LR) THEN
            MEMRY(IP) = LT + 1
            MEMRY(LT+1) = IPN
            IF (IPN.EQ.0) THEN
               MEMRY(2) = LT + 1
            END IF

            IB = IP + 2
            DO 2670 J = 1,LR
               MEMRY(IP+1+J) = 0
 2670       CONTINUE

         ELSE IF (LAV.EQ.LR) THEN
            MEMRY(IP) = -IPN
            IB = IP + 2
            DO 2690 J = 1,LR
               MEMRY(IP+1+J) = 0
 2690       CONTINUE
         END IF

      END IF

      IP = ABS(IPN)
      GO TO 2650

 2660 CONTINUE
      MEMALL = (IB)
      RETURN

      END
      LOGICAL FUNCTION MEMFRE(IB,MEMRY)
      INTEGER IB
      INTEGER MEMRY(*)

      IPR = IB - 2
      MEMFRE = .FALSE.
      IP = 3
      IQ = 3
      IPN = MEMRY(IP)
 2710 IF (.NOT. (IPN.NE.0)) GO TO 2720
 2730 IF (.NOT. (IPN.LT.0)) GO TO 2740
      IPM = MEMRY(-IPN)
      IF (IPM.GT.0) THEN
         IQ = IP
         IP = -IPN
         IPN = IPM

      ELSE IF (IPM.EQ.0) THEN
         MEMRY(IP) = 0
         MEMRY(2) = IP
         RETURN

      ELSE IF (IPM.LT.0) THEN
         IPN = IPM
         MEMRY(IP) = IPN
      END IF

      GO TO 2730

 2740 CONTINUE
      IF (IP.EQ.IPR) THEN
         IPN = -IPN
         IF (MEMRY(IQ).LT.0) THEN
            IP = IQ
         END IF

         MEMRY(IP) = IPN
         MEMFRE = .TRUE.
         IB = 0

      ELSE
         IQ = IP
         IP = ABS(IPN)
         IPN = MEMRY(IP)
      END IF

      GO TO 2710

 2720 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION MEMCMP(IB,L,MEMRY)
      INTEGER IB
      INTEGER L
      INTEGER MEMRY(*)
      LOGICAL MEMFRE

      IF (L.LE.0) THEN
         MEMCMP = MEMFRE(IB,MEMRY)
         RETURN

      END IF

      IPX = ABS(IB) + L + MOD(L,2)
      IP = 3
 2750 CONTINUE
      IPN = MEMRY(ABS(IP))
      IF (ABS(IP)+2.EQ.ABS(IB)) THEN
         IF (IPN.GT.0) THEN
            IF (IPN.LT.IPX) THEN
               MEMCMP = (.FALSE.)
               RETURN

            END IF

            MEMRY(ABS(IP)) = IPX
            IF (IPN.GT.IPX) THEN
               IF (MEMRY(IPN).GT.0) THEN
                  MEMRY(IPX) = -IPN

               ELSE IF (MEMRY(IPN).EQ.0) THEN
                  MEMRY(IPX) = 0

               ELSE IF (MEMRY(IPN).LT.0) THEN
                  MEMRY(IPX) = MEMRY(IPN)
               END IF

            END IF

         END IF

         MEMCMP = (.TRUE.)
         RETURN

      END IF

      IP = IPN
 2760 IF (.NOT. (IPN.EQ.0)) GO TO 2750
      MEMCMP = (.FALSE.)
      RETURN

      END
      LOGICAL FUNCTION MEMEXT(IB,MEMRY,INCL)
      INTEGER IB
      INTEGER MEMRY(*)
      INTEGER DELTA
      INTEGER MEMALL
      LOGICAL MEMFRE
      INTEGER LA
      INTEGER IBT
C ... Guess by GDS -- DELTA not defined, set to 0 since that is 
C     what it would be on VMS if undefined.
      DELTA = 0
      LA = MEMRY(IB-2) - IB
      IBT = IB
      INCL = LA + DELTA
      IB = MEMALL(INCL,MEMRY)
      DO 2780 I = 1,LA
         MEMRY(IB+I-1) = MEMRY(IBT+I-1)
 2780 CONTINUE
      IF (.NOT.MEMFRE(IBT,MEMRY)) THEN
         WRITE (6,*) 'Error deallocating old segment.'
         CALL CPUTBK(.TRUE.)
         MEMEXT = .FALSE.
         RETURN

      END IF

      MEMEXT = .TRUE.
      RETURN

      END
      SUBROUTINE MEMTRC(MEMRY)
      INTEGER MEMRY(*)

      IKL = MEMRY(1)
      IKM = MEMRY(2)
      WRITE (6,*) IKL,IKM
      IP = 3
 2800 CONTINUE
      IPN = MEMRY(ABS(IP))
      PRINT *,IP,IPN
      IP = IPN
 2810 IF (.NOT. (IPN.EQ.0)) GO TO 2800
      RETURN

      END
      SUBROUTINE MXMULT(N,MAT1,MAT2,MATR)
      REAL MAT1(N,*),MAT2(N,*),MATR(N,*)

      DO 2830 J = 1,N
         DO 2850 I = 1,N
            MATR(I,J) = 0.
            DO 2870 K = 1,N
               MATR(I,J) = MATR(I,J) + MAT1(I,K)*MAT2(K,J)
 2870       CONTINUE
 2850    CONTINUE
 2830 CONTINUE
      RETURN

      END
      SUBROUTINE MXCOPY(N,MAT1,MAT2)
      REAL MAT1(N,*),MAT2(N,*)

      DO 2890 J = 1,N
         DO 2910 I = 1,N
            MAT2(I,J) = MAT1(I,J)
 2910    CONTINUE
 2890 CONTINUE
      RETURN

      END
      SUBROUTINE MXIDEN(N,MAT)
      REAL MAT(N,*)

      DO 2930 J = 1,N
         DO 2950 I = 1,N
            MAT(I,J) = 0.
 2950    CONTINUE
 2930 CONTINUE
      DO 2970 I = 1,N
         MAT(I,I) = 1.
 2970 CONTINUE
      RETURN

      END
      SUBROUTINE MXVECT(N,VEC,MAT,RES)
      REAL VEC(*),MAT(N,*),RES(*)

      DO 2990 J = 1,N
         RES(J) = 0.0
         DO 3010 I = 1,N
            RES(J) = RES(J) + MAT(I,J)*VEC(I)
 3010    CONTINUE
 2990 CONTINUE
      RETURN

      END
      SUBROUTINE MXZERO(N,MAT)
      REAL MAT(N,*)

      DO 3030 J = 1,N
         DO 3050 I = 1,N
            MAT(I,J) = 0.
 3050    CONTINUE
 3030 CONTINUE
      RETURN

      END
      SUBROUTINE SIORPT(MODULE,MESS,DISP)
      IMPLICIT INTEGER (A-Z)
      CHARACTER*(*) MODULE
      CHARACTER*(*) MESS
      CHARACTER*80 LOCLIN
      CHARACTER*20 ERRORT
      CHARACTER*120 ERRMSG
      CHARACTER*8 ERRMOD

      ERRMOD = MODULE
      ERRMSG = MESS
      CALL CPUDAC(ERRORT)
      CALL CHRTRM(ERRORT,LT)
      CALL CHRTRM(ERRMOD,L)
      IF (L.EQ.0) THEN
         L = 1
      END IF

      IF (DISP.EQ.1) THEN

      ELSE IF (DISP.EQ.2) THEN
         LOCLIN = '#PLT error (warning) in module '//ERRMOD(1:L)//
     *            ' at '//ERRORT(1:LT)
         WRITE (6,10) LOCLIN

   10    FORMAT (1X,A)

         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN

      ELSE IF (DISP.EQ.3) THEN
         LOCLIN = '#PLT error (traceback) in module '//ERRMOD(1:L)//
     *            ' at '//ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN
         CALL CPUTBK(.FALSE.)

      ELSE IF (DISP.EQ.4) THEN
         LOCLIN = '#PLT error (fatal) in module '//ERRMOD(1:L)//' at '//
     *            ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN
         CALL CPUTBK(.TRUE.)

      ELSE
         LOCLIN = '#PLT error (fatal) in module '//ERRMOD(1:L)//' at '//
     *            ERRORT(1:LT)
         WRITE (6,10) LOCLIN
         LOCLIN = '#Error message: '//ERRMSG
         WRITE (6,10) LOCLIN
         CALL CPUTBK(.TRUE.)
      END IF

      RETURN

      END
      SUBROUTINE VECRGP(ND,V,VMAX,VMIN)
      DIMENSION V(*)

      IF (ND.LT.1) THEN
         RETURN

      END IF

      VMAX = V(1)
      DO 3070 I = 1,ND
         VMAX = MAX(VMAX,V(I))
 3070 CONTINUE
      IF (VMAX.LT.0) THEN
         VMAX = 1.
         VMIN = .1
         RETURN

      END IF

      VMIN = VMAX
      DO 3090 I = 1,ND
         IF (V(I).GT.0.) THEN
            VMIN = MIN(VMIN,V(I))
         END IF

 3090 CONTINUE
      IF (VMIN.EQ.VMAX) THEN
         VMIN = .1*VMAX
      END IF

      RETURN

      END
      SUBROUTINE VECRGS(ND,V,VMAX,VMIN)
      DIMENSION V(*)

      IF (ND.LT.1) THEN
         RETURN

      END IF

      VMAX = V(1)
      VMIN = V(1)
      I = 1
 3110 IF (.NOT. (I.LE.ND)) GO TO 3130
      T = V(I)
      VMAX = MAX(VMAX,T)
      VMIN = MIN(VMIN,T)
 3120 I = I + 1
      GO TO 3110

 3130 CONTINUE
      RETURN

      END
      LOGICAL FUNCTION TTYIFC(RESET)
      LOGICAL RESET

      TTYIFC = .FALSE.
      RETURN

      END
