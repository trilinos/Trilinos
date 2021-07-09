C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
