C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
