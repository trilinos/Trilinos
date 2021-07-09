C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CGCAL2(CRD,IX,MAT,MASS,VOL,DENS,VOLM,CG,ZITOT,XXX,
     *    XG,XI,XINI,AJ,NNODES,NDIM,NQUAD,VOLMN,IELM,NELBLK,
     *    AXI,NUMNP)

      DIMENSION CRD(NUMNP,*), IX(NNODES,*), MAT(6,*), MASS(*),
     *    DENS(*), VOLM(*), CG(*), ZITOT(*),VOLMN(4,*),IELM(4,*),
     *    XXX(NDIM+1,NQUAD,*),XG(NDIM,*), XI(NDIM,*), XINI(*),
     *    AJ(2,*)
      DIMENSION ZI(6), ZMOM(3)
      DIMENSION CCC(2,4)

      LOGICAL AXI
      REAL MASS, MASSE
      PI = ATAN2(0.0, -1.0)

C ... VOLMN(1,*) = MINIMUM VOLUME (AREAS FOR 2-D)
C     VOLMN(2,*) = MAXIMUM VOLUME
C     VOLMN(3,*) = TOTAL VOLUME

      DO 10 I=1, NELBLK
          VOLMN(1,I) = 1.0E30
          VOLMN(2,I) = 0.0
          VOLMN(3,I) = 0.0
          IELM (3,I) = 0
          MASS(I)    = 0.0
          VOLM(I)    = 0.0
   10 CONTINUE
      DO 20 I=1,3
          ZITOT(I)   = 0.0
          ZITOT(I+3) = 0.0
          ZMOM(I)    = 0.0
   20 CONTINUE
      ZMAS = 0.0
      VOL  = 0.0

C ... GET QUADRATURE POINT LOCATIONS, EVALUATE SHAPE FUNCTIONS

      CALL QUAD(XXX, XI, XG, NDIM, NNODES, NQUAD, WT)

      DO 110 IBLK = 1, NELBLK
          IF (MAT(5,IBLK) .NE. 1) GOTO 110
          IELBEG = MAT(3,IBLK)
          IELEND = MAT(4,IBLK)
          MIEL   = IBLK
          DO 100 IEL = IELBEG, IELEND

C ... CALCULATE AREA, VOLUME, AND MOMENTS OF INERTIA OF ELEMENT

              DO 30 I=1,3
                  ZI(I)   = 0.0
                  ZI(I+3) = 0.0
                  CG(I)   = 0.0
   30         CONTINUE
              VOLUME = 0.0

              DO 40 I=1,4
                  CCC(1,I) = CRD(IX(I,IEL),1)
                  CCC(2,I) = CRD(IX(I,IEL),2)
   40         CONTINUE
              DO 80 NG=1,NQUAD
                  DET = 0.0
                  DO 60 J=1,2
                      XINI(J) = 0.0
                      DO 50 K=1,2
                          AJ(K,J) = 0.0
   50                 CONTINUE
   60             CONTINUE

                  DO 70 I=1,4
                      XINI(1) = XINI(1)+XXX(1,I,NG) * CCC(1,I)
                      AJ(1,1) = AJ(1,1)+XXX(2,I,NG) * CCC(1,I)
                      AJ(2,1) = AJ(2,1)+XXX(3,I,NG) * CCC(1,I)

                      XINI(2) = XINI(2)+XXX(1,I,NG) * CCC(2,I)
                      AJ(1,2) = AJ(1,2)+XXX(2,I,NG) * CCC(2,I)
                      AJ(2,2) = AJ(2,2)+XXX(3,I,NG) * CCC(2,I)
   70             CONTINUE

                  DET = ( AJ(1,1)*AJ(2,2) - AJ(2,1)*AJ(1,2) ) * WT
                  DETW = DET * DENS(MIEL)

                  IF (AXI) THEN

C ... CG(1) IS THE VOLUME FOR AXI 2-D, VOLUME IS ACTUALLY C/S AREA

                      CG(1) = CG(1) + DET * XINI(1)
                      CG(2) = CG(2) + DETW * XINI(1) * XINI(2)
                      ZI(2) = ZI(2) + DETW * XINI(1)**3
                      ZI(1) = ZI(1) + DETW * XINI(1)*XINI(2)**2
                      VOLUME  = VOLUME  + DET
                  ELSE
                      CG(1) = CG(1) + DETW * XINI(1)
                      CG(2) = CG(2) + DETW * XINI(2)
                      ZI(1) = ZI(1) + DETW * XINI(2)**2
                      ZI(2) = ZI(2) + DETW * XINI(1)**2
                      ZI(3) = ZI(3) + DETW * XINI(1)*XINI(2)
                      VOLUME = VOLUME + DET
                  END IF

   80         CONTINUE

C ... DETERMINE MIN/MAX ELEMENT VOLUMES FOR EACH MATERIAL AND
C        COUNT NUMBER OF ELEMENTS FOR EACH MATERIAL

              IELM(3,MIEL)      = IELM(3,MIEL) + 1
              VOLMN(3,MIEL)     = VOLMN(3,MIEL) + VOLUME
              IF (VOLUME .LT. VOLMN(1,MIEL)) THEN
                  VOLMN(1,MIEL) = VOLUME
                  IELM(1,MIEL)  = IEL
              ELSE IF (VOLUME .GT. VOLMN(2,MIEL)) THEN
                  VOLMN(2,MIEL) = VOLUME
                  IELM(2,MIEL)  = IEL
              ENDIF

              IF (AXI) THEN
                  VOLUME = 2. * PI * CG(1)
                  ZI(2)  = ZI(2) * 2. * PI
                  ZI(1)  = ZI(1) * 2. * PI + ZI(2) / 2.0
                  ZI(3)  = ZI(1)
              END IF
              DO 90 I=1,3
                  ZITOT(I) = ZITOT(I) + ZI(I)
   90         CONTINUE
              MASSE = VOLUME * DENS(MIEL)
              MASS(MIEL)= MASS(MIEL) + MASSE
              VOLM(MIEL)= VOLM(MIEL) + VOLUME
              VOL  = VOL  + VOLUME
              ZMAS = ZMAS + MASSE
              IF (AXI) THEN
                  CG(1)   = 0.0
                  ZMOM(2) = ZMOM(2) + CG(2) * 2. * PI
              ELSE
                  ZMOM(1) = ZMOM(1) + CG(1)
                  ZMOM(2) = ZMOM(2) + CG(2)
              END IF
  100     CONTINUE
  110 CONTINUE
      FIX = SIGN(0.5, ZMAS) + SIGN(0.5, -ZMAS)
      DO 120 I=1,3
          CG(I) = ZMOM(I) / (ZMAS + FIX)
  120 CONTINUE
      IF (AXI) THEN
          ZITOT(1) = ZITOT(1) - ZMAS * CG(2)**2
          ZITOT(3) = ZITOT(3) - ZMAS * CG(2)**2
      ELSE
          ZITOT(1) = ZITOT(1) - ZMAS * CG(2)**2
          ZITOT(2) = ZITOT(2) - ZMAS * CG(1)**2
          ZITOT(4) = ZITOT(3) - ZMAS * CG(1) * CG(2)
          ZITOT(3) = ZITOT(1) + ZITOT(2)
      END IF

      RETURN
      END
