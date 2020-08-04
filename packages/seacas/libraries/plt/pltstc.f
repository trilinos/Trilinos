C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION PLTSTC(INDX,BUFF)
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

      PLTSTC = .TRUE.
      IF (INDX.EQ.0) THEN
         CALL PLTRSC

      ELSE IF (INDX.EQ.1) THEN
         BTEMP = COLP(1)
         COLP(1) = BUFF(1)
         IF (BUFF(1).EQ.-1.) THEN
            CALL PLTSPC(0.,.706,0.,.706,.1875,0.,0.,1.)
            CALL PLTSPC(.1875,0.,0.,1.,.3708,0.,1.,0.)
            CALL PLTSPC(.3708,0.,1.,0.,.6208,1.,1.,0.)
            CALL PLTSPC(.6208,1.,1.,0.,.8292,1.,.659,0.)
            CALL PLTSPC(.8292,1.,.659,0.,1.,1.,0.,0.)

         ELSE IF (BUFF(1).EQ.0.) THEN
            CALL PLTSPC(0.,1.,0.,0.,.1708,1.,.659,0.)
            CALL PLTSPC(.1708,1.,.659,0.,.3792,1.,1.,0.)
            CALL PLTSPC(.3792,1.,1.,0.,.6292,0.,1.,0.)
            CALL PLTSPC(.6292,0.,1.,0.,.8125,0.,0.,1.)
            CALL PLTSPC(.8125,0.,0.,1.,1.,.706,0.,.706)

         ELSE IF (BUFF(1).EQ.1.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,2),PALETT(2,2),
     *                  PALETT(3,2))

         ELSE IF (BUFF(1).EQ.2.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,3),PALETT(2,3),
     *                  PALETT(3,3))

         ELSE IF (BUFF(1).EQ.3.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,4),PALETT(2,4),
     *                  PALETT(3,4))

         ELSE IF (BUFF(1).EQ.4.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,5),PALETT(2,5),
     *                  PALETT(3,5))

         ELSE IF (BUFF(1).EQ.5.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,6),PALETT(2,6),
     *                  PALETT(3,6))

         ELSE IF (BUFF(1).EQ.6.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,7),PALETT(2,7),
     *                  PALETT(3,7))

         ELSE IF (BUFF(1).EQ.7.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,8),PALETT(2,8),
     *                  PALETT(3,8))

         ELSE IF (BUFF(1).EQ.8.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,9),PALETT(2,9),
     *                  PALETT(3,9))

         ELSE IF (BUFF(1).EQ.9.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,10),PALETT(2,10),
     *                  PALETT(3,10))

         ELSE IF (BUFF(1).EQ.10.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,11),PALETT(2,11),
     *                  PALETT(3,11))

         ELSE IF (BUFF(1).EQ.11.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,12),PALETT(2,12),
     *                  PALETT(3,12))

         ELSE IF (BUFF(1).EQ.12.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,13),PALETT(2,13),
     *                  PALETT(3,13))

         ELSE IF (BUFF(1).EQ.13.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,14),PALETT(2,14),
     *                  PALETT(3,14))

         ELSE IF (BUFF(1).EQ.14.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,15),PALETT(2,15),
     *                  PALETT(3,15))

         ELSE IF (BUFF(1).EQ.15.) THEN
            CALL PLTSPC(0.,0.,0.,0.,1.,PALETT(1,16),PALETT(2,16),
     *                  PALETT(3,16))

         ELSE
            COLP(1) = BTEMP
            CALL CHRRVC(BUFF(1),IERROR,L)
            CALL PLTFLU
            CALL SIORPT('PLTSTC','Illegal buffer '//IERROR(1:L)//
     *                  ' passed to PLTSTC.',2)
            PLTSTC = .FALSE.
         END IF

      ELSE
         CALL CHRIC(INDX,IERROR,L)
         CALL PLTFLU
         CALL SIORPT('PLTSTC','Illegal index '//IERROR(1:L)//'.',2)
         PLTSTC = .FALSE.
      END IF

      RETURN

      END
