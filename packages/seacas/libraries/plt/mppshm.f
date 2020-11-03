C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MPPSHM()
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPPSHM')
      DIMENSION EQMAP(193)
      EQUIVALENCE (EQMAP(1),MODEL(1,1))
      COMMON /MPSTCK/SVMAP(195,10),MAPDEP
      EXTERNAL PLTBLK

      MPPSHM = .FALSE.
      IF (MAPDEP.EQ.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Stack is full; depth is ten',2)
         RETURN

      END IF

      MPPSHM = .TRUE.
      MAPDEP = MAPDEP + 1
      DO 2960 I = 1,193
         SVMAP(I,MAPDEP) = EQMAP(I)
 2960 CONTINUE
      SVMAP(194,MAPDEP) = DBLE(NCPLIN)
      SVMAP(195,MAPDEP) = DBLE(NCPLAN)
      RETURN

      END
