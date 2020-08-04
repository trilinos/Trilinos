C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
