C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MPWIND(LEFT,RIGHT,BOTTOM,TOP,NEAR,FAR)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      REAL LEFT,NEAR
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPWIND')

      MPWIND = .FALSE.
      IF (NEAR.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *         'You cannot specify the near clipping plane less than 0.'
     *               ,2)
         RETURN

      END IF

      IF (NEAR.GE.FAR) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the near clipping plane >= to the far clipping
     * plane',2)
         RETURN

      END IF

      IF (RIGHT.EQ.LEFT) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the right and left edges of the viewing frustr
     *um as equal',2)
         RETURN

      END IF

      IF (TOP.EQ.BOTTOM) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the top and bottom edges of the viewing frustr
     *um as equal',2)
         RETURN

      END IF

      MPWIND = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(1,1) = 2.*NEAR/ (RIGHT-LEFT)
      PROJ(2,2) = 2.*NEAR/ (TOP-BOTTOM)
      PROJ(3,1) = (RIGHT+LEFT)/ (RIGHT-LEFT)
      PROJ(3,2) = (TOP+BOTTOM)/ (TOP-BOTTOM)
      PROJ(3,3) = - (FAR+NEAR)/ (FAR-NEAR)
      PROJ(3,4) = -1.
      PROJ(4,3) = - (2.*FAR*NEAR)/ (FAR-NEAR)
      CPNEAR = NEAR
      CPFAR = FAR
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
