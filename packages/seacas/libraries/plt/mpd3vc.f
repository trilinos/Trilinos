C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MPD3VC(N,X0,Y0,Z0,X1,Y1,Z1)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),X0(*),Y0(*),Z0(*),X1(*),Y1(*),Z1(*)
      INTEGER MASK(1)

      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      J = 0
      KM = 0
 2840 IF (.NOT. (J.LT.N)) GO TO 2850
      JN = MIN(N-J,32)
      J1 = J
      KM = KM + 1
      J = J + JN
      MASK(1) = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MVP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MVP,TARR5,TARR6,
     *               TARR7,TARR8)

      ELSE
         CALL MPMUL3(JN,X0(1+J1),Y0(1+J1),Z0(1+J1),MODEL,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL3(JN,X1(1+J1),Y1(1+J1),Z1(1+J1),MODEL,TARR5,TARR6,
     *               TARR7,TARR8)
         MASK(1) = -1
         DO 2860 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCV3(JN,MASK,TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,
     *                  TARR1,TARR2,TARR3,TARR5,TARR6,TARR7,C1,C2)
 2860    CONTINUE
         CALL MPMUL4(JN,MASK(1),TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
         CALL MPMUL4(JN,MASK(1),TARR5,TARR6,TARR7,TARR8,VP,TARR5,TARR6,
     *               TARR7,TARR8)
      END IF

      CALL PLTZCV(CPNEAR,CPFAR,JN,MASK,TARR1,TARR2,TARR4,TARR5,TARR6,
     *            TARR8)
      DO 2880 I = 1,JN
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
         TARR5(I) = (TARR5(I)/TARR8(I))*VSX + VCX
         TARR6(I) = (TARR6(I)/TARR8(I))*VSY + VCY
 2880 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      CALL PLTVWV(C1,C2,JN,MASK,TARR1,TARR2,TARR5,TARR6)
      CALL PLTVCM(JN,MASK,TARR1,TARR2,TARR5,TARR6)
      GO TO 2840

 2850 CONTINUE
      RETURN

      END
