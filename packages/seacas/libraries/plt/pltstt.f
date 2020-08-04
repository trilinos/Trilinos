C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
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
