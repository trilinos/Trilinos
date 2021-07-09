C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTRSC
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

      COLP(2) = MIN(16.,DEVCAP(4))
      COLP(3) = DEVCAP(4) - COLP(2)
      PALETT(1,1) = 0.
      PALETT(2,1) = 0.
      PALETT(3,1) = 0.
      CALL VDSTCO(1,0,PALETT(1,1),0)
      PALETT(1,2) = 1.
      PALETT(2,2) = 0.
      PALETT(2,3) = 0.
      CALL VDSTCO(1,1,PALETT(1,2),0)
      PALETT(1,3) = 0.
      PALETT(2,3) = 1.
      PALETT(3,3) = 0.
      CALL VDSTCO(1,2,PALETT(1,3),0)
      PALETT(1,4) = 1.
      PALETT(2,4) = 1.
      PALETT(3,4) = 0.
      CALL VDSTCO(1,3,PALETT(1,4),0)
      PALETT(1,5) = 0.
      PALETT(2,5) = 0.
      PALETT(3,5) = 1.
      CALL VDSTCO(1,4,PALETT(1,5),0)
      PALETT(1,6) = 1.
      PALETT(2,6) = 0.
      PALETT(3,6) = 1.
      CALL VDSTCO(1,5,PALETT(1,6),0)
      PALETT(1,7) = 0.
      PALETT(2,7) = 1.
      PALETT(3,7) = 1.
      CALL VDSTCO(1,6,PALETT(1,7),0)
      PALETT(1,8) = 1.
      PALETT(2,8) = 1.
      PALETT(3,8) = 1.
      CALL VDSTCO(1,7,PALETT(1,8),0)
      PALETT(1,9) = .4
      PALETT(2,9) = .4
      PALETT(3,9) = .4
      CALL VDSTCO(1,8,PALETT(1,9),0)
      PALETT(1,10) = .7
      PALETT(2,10) = .7
      PALETT(3,10) = .7
      CALL VDSTCO(1,9,PALETT(1,10),0)
      PALETT(1,11) = .225
      PALETT(2,11) = .225
      PALETT(3,11) = .225
      CALL VDSTCO(1,10,PALETT(1,11),0)
      PALETT(1,12) = 1.
      PALETT(2,12) = .35
      PALETT(3,12) = .45
      CALL VDSTCO(1,11,PALETT(1,12),0)
      PALETT(1,13) = .5
      PALETT(2,13) = 1.
      PALETT(3,13) = .8
      CALL VDSTCO(1,12,PALETT(1,13),0)
      PALETT(1,14) = .4
      PALETT(2,14) = .7
      PALETT(3,14) = 1.
      CALL VDSTCO(1,13,PALETT(1,14),0)
      PALETT(1,15) = .706
      PALETT(2,15) = 0.
      PALETT(3,15) = .706
      CALL VDSTCO(1,14,PALETT(1,15),0)
      PALETT(1,16) = 1.
      PALETT(2,16) = .659
      PALETT(3,16) = 0.
      CALL VDSTCO(1,15,PALETT(1,16),0)
      CALL PLTSTC(1,0.)
      RETURN

      END
