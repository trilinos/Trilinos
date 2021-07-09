C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTINT
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
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      IF (FIRST) THEN
         ASP = 10./7.5
         CALL VDINIT(ASP,5)
         CALL VDFRAM(0)
         CALL PLTIQD(DEVCAP)
         CALL VDSTCS(1./58.)
         CALL VDIQOS(DEFOUT)
         DEFOUT(1) = 7.
         FIRST = .FALSE.
      END IF

      CALL VDIQND(DEVP(4),DEVP(5))
      CALL PLTRST
      CALL PLTRSD
      CALL PLTRSV
      CALL PLTRSG
      CALL PLTRSM
      CALL PLTRSC
      CALL MPINIT
      CALL MEMINI(MEMORY,1000)
      CALL PLTFLU
      RETURN

      END
