C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MPD2PG(N,XV,YV,MODE)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(2),C2(2),XV(*),YV(*)
      CHARACTER MODE* (*),TMODE*1

      IF (N.GT.32) THEN
         CALL PLTFLU
         CALL SIORPT('MPD2PG','Too many vertices specified',2)
         RETURN

      END IF

      TMODE = MODE
      CALL CHRDN(TMODE,TMODE)
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLIN.EQ.0)
      IF (NOCLIP) THEN
         CALL MPMUL2(N,XV,YV,MVP,TARR1,TARR2,TARR3,TARR4)
         NO = N

      ELSE
         CALL MPMUL2(N,XV,YV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2460 K = 1,NCPLIN
            C1(1) = CPLINE(1,1,K)
            C1(2) = CPLINE(1,2,K)
            C2(1) = CPLINE(2,1,K)
            C2(2) = CPLINE(2,2,K)
            NO = 32
            CALL PLTCG2(N,TARR1,TARR2,NO,TARR5,TARR6,C1,C2)
 2460    CONTINUE
         IF (NO.GT.N) THEN
            DO 2480 J = N + 1,NO
               TARR4(J) = 1.
 2480       CONTINUE
         END IF

         CALL MPMUL4(NO,-1,TARR5,TARR6,TARR3,TARR4,VP,TARR1,TARR2,TARR3,
     *               TARR4)
      END IF

      DO 2500 I = 1,NO
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2500 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NVO = 32
      CALL PLTVWG(C1,C2,NO,TARR1,TARR2,TARR3,NVO,TARR5,TARR6,TARR7)
      IF (TMODE.EQ.'s') THEN
         CALL PLTPLY(NVO,TARR5,TARR6)

      ELSE IF (TMODE.EQ.'o') THEN
         CALL PLTMOV(TARR5(1),TARR6(1))
         DO 2520 J = 2,NVO
            CALL PLTDRW(TARR5(J),TARR6(J))
 2520    CONTINUE
         CALL PLTDRW(TARR5(1),TARR6(1))

      ELSE
         CALL PLTFLU
         CALL SIORPT('MPD2PG','Unrecognized drawing mode: '//TMODE,2)
      END IF

      RETURN

      END
