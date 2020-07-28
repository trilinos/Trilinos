C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MPGETM(TMODEL,TVIEW,TPROJ,TVWPT)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      DIMENSION TMODEL(4,4),TVIEW(4,4),TPROJ(4,4),TVWPT(4)

      DO 2900 I = 1,4
         DO 2920 J = 1,4
            TMODEL(J,I) = MODEL(J,I)
            TPROJ(J,I) = PROJ(J,I)
            TVIEW(J,I) = VIEW(J,I)
 2920    CONTINUE
         TVWPT(I) = VWPORT(I)
 2900 CONTINUE
      RETURN

      END
