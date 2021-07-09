C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CGCAL3(CRD,IX,MAT,MASS,VOL,DENS,VOLM,CG,ZITOT,XXX,
     *    XG,XI,XINI,AJ,NNODES,NDIM,NQUAD,VOLMN,IELM,NELBLK,NUMNP)

C ... CALCULATE PROPERTIES FOR THREE-DIMENSIONAL MESH --- BRICKS ONLY

      DIMENSION CRD(NUMNP,*), IX(NNODES,*), MAT(6,*), MASS(*),
     *    DENS(*), VOLM(*), CG(*), ZITOT(*),VOLMN(4,*),IELM(4,*),
     *    XXX(3+1,NQUAD,*),XG(3,*), XI(3,*), XINI(*),
     *    AJ(3,*)
      DIMENSION ZI(6), ZMOM(3)
      DIMENSION CCC(3,8)

      REAL MASS, MASSE

C ... VOLMN(1,*) = MINIMUM VOLUME (AREAS FOR 2-D)
C     VOLMN(2,*) = MAXIMUM VOLUME
C     VOLMN(3,*) = TOTAL VOLUME
C     VOLMN(4,*) = MAXIMUM SUM OF SQUARES OF INVERSE LENGTHS

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

      DO 80 IBLK = 1, NELBLK
          IF (MAT(5,IBLK) .NE. 1) GOTO 80
          IELBEG = MAT(3,IBLK)
          IELEND = MAT(4,IBLK)
          MIEL   = IBLK
          DO 70 IEL = IELBEG, IELEND

C ... CALCULATE ARA, VOLUME, AND MOMENTS OF INERTIA OF ELEMENT

              DO 30 I=1,3
                  ZI(I)   = 0.0
                  ZI(I+3) = 0.0
                  CG(I)   = 0.0
   30         CONTINUE
              VOLUME = 0.0

              DO 40 I=1,8
                  CCC(1,I) = CRD(IX(I,IEL),1)
                  CCC(2,I) = CRD(IX(I,IEL),2)
                  CCC(3,I) = CRD(IX(I,IEL),3)
   40         CONTINUE

              DO 60 NG=1,NQUAD
                  DET = 0.0

                  XINI(1) = XXX(1,1,NG) * CCC(1,1) +
     *                XXX(1,2,NG) * CCC(1,2) + XXX(1,3,NG) * CCC(1,3) +
     *                XXX(1,4,NG) * CCC(1,4) + XXX(1,5,NG) * CCC(1,5) +
     *                XXX(1,6,NG) * CCC(1,6) + XXX(1,7,NG) * CCC(1,7) +
     *                XXX(1,8,NG) * CCC(1,8)
                  XINI(2) = XXX(1,1,NG) * CCC(2,1) +
     *                XXX(1,2,NG) * CCC(2,2) + XXX(1,3,NG) * CCC(2,3) +
     *                XXX(1,4,NG) * CCC(2,4) + XXX(1,5,NG) * CCC(2,5) +
     *                XXX(1,6,NG) * CCC(2,6) + XXX(1,7,NG) * CCC(2,7) +
     *                XXX(1,8,NG) * CCC(2,8)
                  XINI(3) = XXX(1,1,NG) * CCC(3,1) +
     *                XXX(1,2,NG) * CCC(3,2) + XXX(1,3,NG) * CCC(3,3) +
     *                XXX(1,4,NG) * CCC(3,4) + XXX(1,5,NG) * CCC(3,5) +
     *                XXX(1,6,NG) * CCC(3,6) + XXX(1,7,NG) * CCC(3,7) +
     *                XXX(1,8,NG) * CCC(3,8)

                  AJ(1,1) = XXX(2,1,NG) * CCC(1,1)
     *                + XXX(2,2,NG) * CCC(1,2) + XXX(2,3,NG) * CCC(1,3)
     *                + XXX(2,4,NG) * CCC(1,4) + XXX(2,5,NG) * CCC(1,5)
     *                + XXX(2,6,NG) * CCC(1,6) + XXX(2,7,NG) * CCC(1,7)
     *                + XXX(2,8,NG) * CCC(1,8)
                  AJ(2,1) = XXX(3,1,NG) * CCC(1,1)
     *                + XXX(3,2,NG) * CCC(1,2) + XXX(3,3,NG) * CCC(1,3)
     *                + XXX(3,4,NG) * CCC(1,4) + XXX(3,5,NG) * CCC(1,5)
     *                + XXX(3,6,NG) * CCC(1,6) + XXX(3,7,NG) * CCC(1,7)
     *                + XXX(3,8,NG) * CCC(1,8)
                  AJ(3,1) = XXX(4,1,NG) * CCC(1,1)
     *                + XXX(4,2,NG) * CCC(1,2) + XXX(4,3,NG) * CCC(1,3)
     *                + XXX(4,4,NG) * CCC(1,4) + XXX(4,5,NG) * CCC(1,5)
     *                + XXX(4,6,NG) * CCC(1,6) + XXX(4,7,NG) * CCC(1,7)
     *                + XXX(4,8,NG) * CCC(1,8)

                  AJ(1,2) = XXX(2,1,NG) * CCC(2,1)
     *                + XXX(2,2,NG) * CCC(2,2) + XXX(2,3,NG) * CCC(2,3)
     *                + XXX(2,4,NG) * CCC(2,4) + XXX(2,5,NG) * CCC(2,5)
     *                + XXX(2,6,NG) * CCC(2,6) + XXX(2,7,NG) * CCC(2,7)
     *                + XXX(2,8,NG) * CCC(2,8)
                  AJ(2,2) = XXX(3,1,NG) * CCC(2,1)
     *                + XXX(3,2,NG) * CCC(2,2) + XXX(3,3,NG) * CCC(2,3)
     *                + XXX(3,4,NG) * CCC(2,4) + XXX(3,5,NG) * CCC(2,5)
     *                + XXX(3,6,NG) * CCC(2,6) + XXX(3,7,NG) * CCC(2,7)
     *                + XXX(3,8,NG) * CCC(2,8)
                  AJ(3,2) = XXX(4,1,NG) * CCC(2,1)
     *                + XXX(4,2,NG) * CCC(2,2) + XXX(4,3,NG) * CCC(2,3)
     *                + XXX(4,4,NG) * CCC(2,4) + XXX(4,5,NG) * CCC(2,5)
     *                + XXX(4,6,NG) * CCC(2,6) + XXX(4,7,NG) * CCC(2,7)
     *                + XXX(4,8,NG) * CCC(2,8)

                  AJ(1,3) = XXX(2,1,NG) * CCC(3,1)
     *                + XXX(2,2,NG) * CCC(3,2) + XXX(2,3,NG) * CCC(3,3)
     *                + XXX(2,4,NG) * CCC(3,4) + XXX(2,5,NG) * CCC(3,5)
     *                + XXX(2,6,NG) * CCC(3,6) + XXX(2,7,NG) * CCC(3,7)
     *                + XXX(2,8,NG) * CCC(3,8)
                  AJ(2,3) = XXX(3,1,NG) * CCC(3,1)
     *                + XXX(3,2,NG) * CCC(3,2) + XXX(3,3,NG) * CCC(3,3)
     *                + XXX(3,4,NG) * CCC(3,4) + XXX(3,5,NG) * CCC(3,5)
     *                + XXX(3,6,NG) * CCC(3,6) + XXX(3,7,NG) * CCC(3,7)
     *                + XXX(3,8,NG) * CCC(3,8)
                  AJ(3,3) = XXX(4,1,NG) * CCC(3,1)
     *                + XXX(4,2,NG) * CCC(3,2) + XXX(4,3,NG) * CCC(3,3)
     *                + XXX(4,4,NG) * CCC(3,4) + XXX(4,5,NG) * CCC(3,5)
     *                + XXX(4,6,NG) * CCC(3,6) + XXX(4,7,NG) * CCC(3,7)
     *                + XXX(4,8,NG) * CCC(3,8)

                  DET = AJ(1,1) * (AJ(2,2)*AJ(3,3) - AJ(2,3)*AJ(3,2))
     *                -AJ(1,2) * (AJ(2,1)*AJ(3,3) - AJ(2,3)*AJ(3,1))
     *                +AJ(1,3) * (AJ(2,1)*AJ(3,2) - AJ(2,2)*AJ(3,1))
                  DET = DET * WT
                  DETW = DET * DENS(MIEL)

                  CG(1) = CG(1) + DETW * XINI(1)
                  CG(2) = CG(2) + DETW * XINI(2)
                  CG(3) = CG(3) + DETW * XINI(3)

                  ZI(1) = ZI(1) + DETW * XINI(1)**2
                  ZI(2) = ZI(2) + DETW * XINI(2)**2
                  ZI(3) = ZI(3) + DETW * XINI(3)**2

                  ZI(4) = ZI(4) + DETW * XINI(1) * XINI(2)
                  ZI(5) = ZI(5) + DETW * XINI(1) * XINI(3)
                  ZI(6) = ZI(6) + DETW * XINI(2) * XINI(3)

                  VOLUME = VOLUME + DET

   60         CONTINUE

C ... DETERMINE MIN/MAX ELEMENT VOLUMES FOR EACH MATERIAL AND
C        COUNT NUMBER OF ELEMENTS FOR EACH MATERIAL

              IELM(3,MIEL)      = IELM(3,MIEL) + 1
              VOLMN(3,MIEL)     = VOLMN(3,MIEL) + VOLUME
              IF (VOLUME .LT. VOLMN(1,MIEL)) THEN
                 VOLMN(1,MIEL) = VOLUME
                 IELM(1,MIEL)  = IEL
              END IF
C ... Changed from else if to if so 1 element and equal size blocks get correct volume
              IF (VOLUME .GT. VOLMN(2,MIEL)) THEN
                 VOLMN(2,MIEL) = VOLUME
                 IELM(2,MIEL)  = IEL
              END IF

              if (volume .le. 0.0) then
                 write (*,*) 'Zero or negative volume at element',
     &                iel
              end if
              ZITOT(1) = ZITOT(1) + ZI(2) + ZI(3)
              ZITOT(2) = ZITOT(2) + ZI(1) + ZI(3)
              ZITOT(3) = ZITOT(3) + ZI(1) + ZI(2)
              ZITOT(4) = ZITOT(4) + ZI(4)
              ZITOT(5) = ZITOT(5) + ZI(5)
              ZITOT(6) = ZITOT(6) + ZI(6)
              MASSE = VOLUME * DENS(MIEL)
              MASS(MIEL)= MASS(MIEL) + MASSE
              VOLM(MIEL)= VOLM(MIEL) + VOLUME
              VOL  = VOL  + VOLUME
              ZMAS = ZMAS + MASSE
              ZMOM(1) = ZMOM(1) + CG(1)
              ZMOM(2) = ZMOM(2) + CG(2)
              ZMOM(3) = ZMOM(3) + CG(3)
   70     CONTINUE
   80 CONTINUE
      FIX = SIGN(0.5, ZMAS) + SIGN(0.5, -ZMAS)
      DO 90 I=1,3
          CG(I) = ZMOM(I) / (ZMAS + FIX)
   90 CONTINUE
      ZITOT(1) = ZITOT(1) - ZMAS * (CG(2)**2 + CG(3)**2)
      ZITOT(2) = ZITOT(2) - ZMAS * (CG(1)**2 + CG(3)**2)
      ZITOT(3) = ZITOT(3) - ZMAS * (CG(1)**2 + CG(2)**2)
      ZITOT(4) = ZITOT(4) - ZMAS *  CG(1)    * CG(2)
      ZITOT(5) = ZITOT(5) - ZMAS *  CG(1)    * CG(3)
      ZITOT(6) = ZITOT(6) - ZMAS *  CG(2)    * CG(3)

      RETURN
      END
