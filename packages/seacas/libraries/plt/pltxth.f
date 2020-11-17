C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTXTH(X,Y,TEXT)
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
      DIMENSION INCHAR(133)
      CHARACTER*(*) TEXT

      CALL PLTRIM(TEXT,L)
      IF (L.LE.0) THEN
         RETURN

      END IF

      NC = 0
      DO 2000 I = 1,L
         IC = ICHAR(TEXT(I:I))
         IF (IC.NE.0) THEN
            NC = NC + 1
            INCHAR(NC) = IC
         END IF

 2000 CONTINUE
      IF (NC.EQ.0) THEN
         RETURN

      END IF

      CALL PLTMOV(X,Y)
      CALL VDTEXT(NC,INCHAR)
      CALL VDIQCP(TEXTP(4),TEXTP(5))
      CALL PLTMOV(X,Y)
      TEXTP(6) = X
      TEXTP(7) = Y
      RETURN

      END
