C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MP2VC(N,X0,Y0,X1,Y1,PX,PY,QX,QY,MASK)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),X0(*),Y0(*),X1(*),Y1(*),PX(*),PY(*),QX(*),
     *          QY(*),MASK(*)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      KM = 0
      J = 0
 2140 IF (.NOT. (J.LT.N)) GO TO 2150
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(KM) = -1
      IF (NOCLIP) THEN
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MVP,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MVP,TARR5,TARR6,TARR7,TARR8)

      ELSE
         CALL MPMUL2(JN,X0(1+J1),Y0(1+J1),MODEL,TARR1,TARR2,TARR3,TARR4)
         CALL MPMUL2(JN,X1(1+J1),Y1(1+J1),MODEL,TARR5,TARR6,TARR7,TARR8)
         DO 2160 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            CALL PLTCV2(JN,MASK(KM),TARR1,TARR2,TARR5,TARR6,TARR1,TARR2,
     *                  TARR5,TARR6,C1,C2)
 2160    CONTINUE
         CALL MPMUL4(JN,MASK(KM),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK(KM),TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      DO 2180 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
 2180 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK(KM),TARR1,TARR2,TARR5,TARR6)
      DO 2200 I = 1,JN
         PX(J1+I) = TARR1(I)
         PY(J1+I) = TARR2(I)
         QX(J1+I) = TARR5(I)
         QY(J1+I) = TARR6(I)
 2200 CONTINUE
      GO TO 2140

 2150 CONTINUE
      RETURN

      END
