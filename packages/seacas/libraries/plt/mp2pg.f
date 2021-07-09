C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MP2PG(N,XV,YV,NO,XVO,YVO)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),XV(*),YV(*),XVO(*),YVO(*)

      IF (N.GT.32) THEN
         CALL PLTFLU
         CALL SIORPT('MP2PG','Too many vertices specified',2)
         RETURN

      END IF

      NOSAVE = NO
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      IF (NOCLIP) THEN
         CALL MPMUL2(N,XV,YV,MVP,TARR1,TARR2,TARR3,TARR4)
         NT = N

      ELSE
         CALL MPMUL2(N,XV,YV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2000 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            NT = 32
            CALL PLTCG2(N,TARR1,TARR2,NT,TARR5,TARR6,C1,C2)
 2000    CONTINUE
         IF (NT.GT.N) THEN
            DO 2020 J = N + 1,NT
               TARR4(J) = 1.
 2020       CONTINUE
         END IF

         CALL MPMUL4(NT,-1,TARR5,TARR6,TARR3,TARR4,VP,TARR1,TARR2,TARR3,
     *               TARR4)
      END IF

      DO 2040 I = 1,NT
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2040 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NO = NOSAVE
      CALL PLTVWG(C1,C2,NT,TARR1,TARR2,TARR3,NO,XVO,YVO,TARR3)
      RETURN

      END
