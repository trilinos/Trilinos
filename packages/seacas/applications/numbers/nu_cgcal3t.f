C    Copyright(C) 2025 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CGCAL3T(CRD,IX,MAT,MASS,VOL,DENS,VOLM,CG,ZITOT,XXX,
     *    XG,XI,XINI,AJ,NNODES,NDIM,NQUAD,VOLMN,IELM,NELBLK,NUMNP)

C ... CALCULATE PROPERTIES FOR THREE-DIMENSIONAL TET MESH

      DIMENSION CRD(NUMNP,*), IX(4,*), MAT(7,*), MASS(*),
     *    DENS(*), VOLM(*), CG(*), ZITOT(*),VOLMN(4,*),IELM(4,*),
     *    XXX(3+1,NQUAD,*),XG(3,*), XI(3,*), XINI(*),
     *    AJ(3,*)
      DIMENSION ZI(6), ZMOM(3)
      DIMENSION CCC(3,8)

      REAL MASS, MASSE
      REAL M(3,3)

C ... VOLMN(1,*) = MINIMUM VOLUME (AREAS FOR 2-D)
C     VOLMN(2,*) = MAXIMUM VOLUME
C     VOLMN(3,*) = TOTAL VOLUME
C     VOLMN(4,*) = MAXIMUM SUM OF SQUARES OF INVERSE LENGTHS

      DO I=1, NELBLK
         VOLMN(1,I) = 1.0E30
         VOLMN(2,I) = 0.0
         VOLMN(3,I) = 0.0
         IELM (3,I) = 0
         MASS(I)    = 0.0
         VOLM(I)    = 0.0
      END DO

      DO I=1,3
         ZITOT(I)   = 0.0
         ZITOT(I+3) = 0.0
         ZMOM(I)    = 0.0
      END DO

      ZMAS = 0.0
      VOL  = 0.0

C ... GET QUADRATURE POINT LOCATIONS, EVALUATE SHAPE FUNCTIONS

      DO IBLK = 1, NELBLK
          IF (MAT(5,IBLK) .NE. 1) CONTINUE
          IELBEG = MAT(3,IBLK)
          IELEND = MAT(4,IBLK)
          MIEL   = IBLK
          DO IEL = IELBEG, IELEND

C ... CALCULATE AREA, VOLUME, AND MOMENTS OF INERTIA OF ELEMENT

             DO I=1,3
                ZI(I)   = 0.0
                ZI(I+3) = 0.0
                CG(I)   = 0.0
             END DO
             VOLUME = 0.0

             X1 = CRD(IX(1,IEL),1)
             X2 = CRD(IX(2,IEL),1)
             X3 = CRD(IX(3,IEL),1)
             X4 = CRD(IX(4,IEL),1)

             Y1 = CRD(IX(1,IEL),2)
             Y2 = CRD(IX(2,IEL),2)
             Y3 = CRD(IX(3,IEL),2)
             Y4 = CRD(IX(4,IEL),2)

             Z1 = CRD(IX(1,IEL),3)
             Z2 = CRD(IX(2,IEL),3)
             Z3 = CRD(IX(3,IEL),3)
             Z4 = CRD(IX(4,IEL),3)

             xini(1) = (x1 + x2 + x3 + x4) / 4.0
             xini(2) = (y1 + y2 + y3 + y4) / 4.0
             xini(3) = (z1 + z2 + z3 + z4) / 4.0

             m(1,1) = x1 - x2
             m(1,2) = y1 - y2
             m(1,3) = z1 - z2

             m(2,1) = x2 - x3
             m(2,2) = y2 - y3
             m(2,3) = z2 - z3

             m(3,1) = x3 - x4
             m(3,2) = y3 - y4
             m(3,3) = z3 - z4

             det = (m(1,1) * (m(2,2) * m(3,3) - m(2,3) * m(3,2)) -
     $            m(2,1) * (m(1,2) * m(3,3) - m(1,3) * m(3,2)) +
     $            m(3,1) * (m(1,2) * m(2,3) - m(1,3) * m(2,2)))

             det = -det / 6.0

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
           END DO
        END DO
      FIX = SIGN(0.5, ZMAS) + SIGN(0.5, -ZMAS)
      DO I=1,3
          CG(I) = ZMOM(I) / (ZMAS + FIX)
       END DO
      ZITOT(1) = ZITOT(1) - ZMAS * (CG(2)**2 + CG(3)**2)
      ZITOT(2) = ZITOT(2) - ZMAS * (CG(1)**2 + CG(3)**2)
      ZITOT(3) = ZITOT(3) - ZMAS * (CG(1)**2 + CG(2)**2)
      ZITOT(4) = ZITOT(4) - ZMAS *  CG(1)    * CG(2)
      ZITOT(5) = ZITOT(5) - ZMAS *  CG(1)    * CG(3)
      ZITOT(6) = ZITOT(6) - ZMAS *  CG(2)    * CG(3)

      RETURN
      END
