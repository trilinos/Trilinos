C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MPROTA(ANGLE,AXIS)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPROTA')
      CHARACTER*1 AXIS,TAXIS
      PARAMETER (DPR=57.2958)

      MPROTA = .FALSE.
      CALL CHRUP(AXIS,TAXIS)
      IF (TAXIS.NE.'X' .AND. TAXIS.NE.'Y' .AND. TAXIS.NE.'Z') THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,'Illegal axis '//AXIS//' specified',2)
         RETURN

      END IF

      MPROTA = .TRUE.
      CALL LDROTA(AXIS,COS(ANGLE/DPR),SIN(ANGLE/DPR),TMAT1)
      CALL MXCOPY(4,MODEL,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,MODEL)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
