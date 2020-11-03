C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
