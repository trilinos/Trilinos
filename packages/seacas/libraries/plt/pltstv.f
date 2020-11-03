C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
            IBUFFT = INT(BUFF(1)) - 1
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
