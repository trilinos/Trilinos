C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
