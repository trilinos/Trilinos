C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
C=======================================================================
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

C .. Same as PLTSTG
      LOGICAL FUNCTION PLTSTG1(INDX,BUFF)
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
      REAL BUFF
      CHARACTER*16 IERROR

      PLTSTG1 = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSG

      ELSE IF (INDX.EQ.1) THEN
         GRAPHP(1) = BUFF

      ELSE IF (INDX.EQ.2) THEN
         GRAPHP(2) = BUFF

      ELSE IF (INDX.EQ.3) THEN
         GRAPHP(3) = BUFF

      ELSE IF (INDX.EQ.4) THEN
         GRAPHP(4) = BUFF

      ELSE IF (INDX.EQ.5) THEN
         GRAPHP(5) = BUFF

      ELSE IF (INDX.EQ.6) THEN
         GRAPHP(38) = BUFF

      ELSE IF (INDX.EQ.7) THEN
         IF (BUFF.EQ.0.) THEN
            GRAPHP(6) = 0.

         ELSE
            GRAPHP(6) = 1.
            GRAPHP(47) = BUFF + 4.
         END IF

      ELSE IF (INDX.EQ.8) THEN
         IF (BUFF.LE.0.) THEN
            CALL PLTFLU
            CALL SIORPT('PLTSTG',
     *                  'Symbol increment must be greater than zero.',2)
            PLTSTG1 = .FALSE.
            RETURN

         END IF

         GRAPHP(23) = BUFF

      ELSE IF (INDX.EQ.9) THEN
         GRAPHP(21) = BUFF

      ELSE IF (INDX.EQ.10) THEN
         GRAPHP(37) = BUFF

      ELSE IF (INDX.EQ.12) THEN
         GRAPHP(32) = BUFF

      ELSE IF (INDX.EQ.13) THEN
         GRAPHP(91) = BUFF

      ELSE IF (INDX.EQ.14) THEN
         GRAPHP(90) = BUFF

      ELSE IF (INDX.EQ.15) THEN
         GRAPHP(35) = BUFF

      ELSE IF (INDX.EQ.16) THEN
         GRAPHP(36) = BUFF

      ELSE IF (INDX.EQ.17) THEN
         GRAPHP(39) = BUFF

      ELSE IF (INDX.EQ.18) THEN
         GRAPHP(40) = BUFF

      ELSE IF (INDX.EQ.19) THEN
         GRAPHP(41) = BUFF

      ELSE IF (INDX.EQ.20) THEN
         GRAPHP(42) = BUFF

      ELSE IF (INDX.EQ.21) THEN
         GRAPHP(92) = BUFF

      ELSE IF (INDX.EQ.22) THEN
         GRAPHP(44) = BUFF

      ELSE IF (INDX.EQ.23) THEN
         GRAPHP(45) = BUFF

      ELSE IF (INDX.EQ.47) THEN
         GRAPHP(88) = BUFF

      ELSE IF (INDX.EQ.48) THEN
         GRAPHP(89) = BUFF

      ELSE IF (INDX.EQ.24) THEN
         GRAPHP(46) = BUFF

      ELSE IF (INDX.EQ.25) THEN
         GRAPHP(48) = BUFF

      ELSE IF (INDX.EQ.26) THEN
         GRAPHP(49) = BUFF

      ELSE IF (INDX.EQ.28) THEN
         GRAPHP(62) = BUFF

      ELSE IF (INDX.EQ.29) THEN
         GRAPHP(63) = BUFF

      ELSE IF (INDX.EQ.30) THEN
         GRAPHP(64) = BUFF

      ELSE IF (INDX.EQ.31) THEN
         GRAPHP(65) = BUFF

      ELSE IF (INDX.EQ.32) THEN
         GRAPHP(66) = BUFF

      ELSE IF (INDX.EQ.33) THEN
         GRAPHP(67) = BUFF

      ELSE IF (INDX.EQ.34) THEN
         GRAPHP(68) = BUFF

      ELSE IF (INDX.EQ.35) THEN
         GRAPHP(69) = BUFF

      ELSE IF (INDX.EQ.36) THEN
         GRAPHP(70) = BUFF

      ELSE IF (INDX.EQ.37) THEN
         GRAPHP(71) = BUFF

      ELSE IF (INDX.EQ.38) THEN
         GRAPHP(72) = BUFF

      ELSE IF (INDX.EQ.39) THEN
         GRAPHP(73) = BUFF

      ELSE IF (INDX.EQ.43) THEN
         GRAPHP(74) = BUFF

      ELSE IF (INDX.EQ.44) THEN
         GRAPHP(75) = BUFF

      ELSE IF (INDX.EQ.45) THEN
         GRAPHP(76) = BUFF

      ELSE IF (INDX.EQ.46) THEN
         GRAPHP(77) = BUFF

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTG','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTG1 = .FALSE.
         RETURN

      END IF

      RETURN

      END

