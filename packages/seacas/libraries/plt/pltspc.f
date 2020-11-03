C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTSPC(S1,RED1,GREEN1,BLUE1,S2,RED2,GREEN2,BLUE2)
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
      DIMENSION DCOLOR(3,256),ICOL(256)

      IF (S1.LT.0. .OR. S1.GE.1.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *               'Starting value of color spectrum is out of range',
     *               2)
         RETURN

      END IF

      IF (S2.LE.0. .OR. S2.GT.1.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *               'Ending value of color spectrum is out of range',2)
         RETURN

      END IF

      IF (S2.LE.S1) THEN
         CALL PLTFLU
         CALL SIORPT('PLTSPC',
     *'Starting value of color spectrum must be less than the ending val
     *ue',2)
         RETURN

      END IF

      IF (COLP(3).EQ.0.) THEN
         RETURN

      END IF

      I1 = INT(COLP(2)) + NINT(S1*COLP(3))
      I2 = INT(COLP(2)) + NINT(S2*COLP(3)-1.)
      IDS = I2 - I1 + 1
      DS = IDS
      DR = (RED2-RED1)/DS
      DG = (GREEN2-GREEN1)/DS
      DB = (BLUE2-BLUE1)/DS
      ICOL(1) = I1
      DCOLOR(1,1) = RED1
      DCOLOR(2,1) = GREEN1
      DCOLOR(3,1) = BLUE1
      DO 2000 I = 2,IDS
         ICOL(I) = ICOL(I-1) + 1
         DCOLOR(1,I) = DCOLOR(1,I-1) + DR
         DCOLOR(1,I) = MAX(0.,DCOLOR(1,I))
         DCOLOR(1,I) = MIN(1.,DCOLOR(1,I))
         DCOLOR(2,I) = DCOLOR(2,I-1) + DG
         DCOLOR(2,I) = MAX(0.,DCOLOR(2,I))
         DCOLOR(2,I) = MIN(1.,DCOLOR(2,I))
         DCOLOR(3,I) = DCOLOR(3,I-1) + DB
         DCOLOR(3,I) = MAX(0.,DCOLOR(3,I))
         DCOLOR(3,I) = MIN(1.,DCOLOR(3,I))
 2000 CONTINUE
      CALL VDSTCO(IDS,ICOL,DCOLOR,0)
      RETURN

      END
