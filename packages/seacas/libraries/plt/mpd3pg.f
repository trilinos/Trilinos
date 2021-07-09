C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MPD3PG(NV,XV,YV,ZV,MODE)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      LOGICAL NOCLIP
      DIMENSION C1(3),C2(3),XV(NV),YV(NV),ZV(NV),XVO(50),YVO(50)
      CHARACTER MODE* (*),TMODE*1

      TMODE = MODE
      CALL CHRDN(TMODE,TMODE)
      VCX = (VWPORT(1)+VWPORT(2))/2.
      VSX = (VWPORT(2)-VWPORT(1))/2.
      VCY = (VWPORT(3)+VWPORT(4))/2.
      VSY = (VWPORT(4)-VWPORT(3))/2.
      NOCLIP = (NCPLAN.EQ.0)
      MASK = -1
      IF (NOCLIP) THEN
         CALL MPMUL3(NV,XV,YV,ZV,MVP,TARR1,TARR2,TARR3,TARR4)

      ELSE
         CALL MPMUL3(NV,XV,YV,ZV,MODEL,TARR1,TARR2,TARR3,TARR4)
         DO 2720 I = 1,NCPLAN
            C1(1) = CPPLAN(1,1,I)
            C1(2) = CPPLAN(1,2,I)
            C1(3) = CPPLAN(1,3,I)
            C2(1) = CPPLAN(2,1,I)
            C2(2) = CPPLAN(2,2,I)
            C2(3) = CPPLAN(2,3,I)
            CALL PLTCP3(NV,MASK,TARR1,TARR2,TARR3,C1,C2)
            IF (MASK.NE.-1) THEN
               GO TO 2720

            END IF

 2720    CONTINUE
         CALL MPMUL4(NV,MASK,TARR1,TARR2,TARR3,TARR4,VP,TARR1,TARR2,
     *               TARR3,TARR4)
      END IF

      CALL PLTZCP(CPNEAR,CPFAR,NV,MASK,TARR4)
      IF (MASK.NE.-1) THEN
         RETURN

      END IF

      DO 2740 I = 1,NV
         TARR1(I) = (TARR1(I)/TARR4(I))*VSX + VCX
         TARR2(I) = (TARR2(I)/TARR4(I))*VSY + VCY
 2740 CONTINUE
      C1(1) = VWPORT(1)
      C1(2) = VWPORT(3)
      C2(1) = VWPORT(2)
      C2(2) = VWPORT(4)
      NO = 50
      CALL PLTVWG(C1,C2,NV,TARR1,TARR2,TARR3,NO,XVO,YVO,TARR4)
      IF (TMODE.EQ.'s') THEN
         CALL PLTPLY(NO,XVO,YVO)

      ELSE IF (TMODE.EQ.'o') THEN
         CALL PLTMOV(XVO(1),YVO(1))
         DO 2760 J = 2,NO
            CALL PLTDRW(XVO(J),YVO(J))
 2760    CONTINUE
         CALL PLTDRW(XVO(1),YVO(1))

      ELSE
         CALL PLTFLU
         CALL SIORPT('MPD3PG','Unrecognized drawing mode: '//TMODE,2)
         RETURN

      END IF

      RETURN

      END
