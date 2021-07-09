C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
