C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
C      CALL VDSTLA(X,Y)
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
