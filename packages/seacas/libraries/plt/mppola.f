C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MPPOLA(DIST,AZIM,INC,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      REAL INC
      PARAMETER (DPR=57.2958)

      ANG = -AZIM/DPR
      CALL LDROTA('y',COS(ANG),SIN(ANG),TMAT1)
      ANG = -INC/DPR
      CALL LDROTA('x',COS(ANG),SIN(ANG),TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,TMAT3)
      ANG = -TWIST/DPR
      CALL LDROTA('z',COS(ANG),SIN(ANG),TMAT1)
      CALL MXMULT(4,TMAT3,TMAT1,TMAT2)
      CALL LDTRAN(0.,0.,-DIST,TMAT1)
      CALL MXMULT(4,TMAT2,TMAT1,VIEW)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RINC = INC/DPR
      RAZIM = AZIM/DPR
      PEYE(1) = DIST*COS(RINC)*SIN(RAZIM)
      PEYE(2) = -DIST*SIN(RINC)
      PEYE(3) = DIST*COS(RINC)*COS(RAZIM)
      PLOOK(1) = 0.
      PLOOK(2) = 0.
      PLOOK(3) = 0.
      ETWIST = TWIST
      RETURN

      END
