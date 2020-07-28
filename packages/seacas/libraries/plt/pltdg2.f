C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTDG2(MAP,N,XV,YV)
      REAL MAP(*)
      INTEGER N
      REAL XV(*),YV(*)
      REAL XWORK(50),YWORK(50),XWORK1(50),YWORK1(50)
      INTEGER NWORK

      AXX = MAP(1)
      AYY = MAP(4)
      AXY = MAP(3)
      AYX = MAP(2)
      BX = MAP(5)
      BY = MAP(6)
      DO 2570 I = 1,N
         XWORK1(I) = AXX*XV(I) + AXY*YV(I) + BX
         YWORK1(I) = AYX*XV(I) + AYY*YV(I) + BY
 2570 CONTINUE
      NWORK = 50
      CALL PLTCG2(N,XWORK1,YWORK1,NWORK,XWORK,YWORK,MAP(7),MAP(9))
      NO = 50
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XWORK1,YWORK1,MAP(9),MAP(11))
      NWORK = 50
      CALL PLTCG2(NO,XWORK1,YWORK1,NWORK,XWORK,YWORK,MAP(11),MAP(13))
      NO = 50
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XWORK1,YWORK1,MAP(13),MAP(7))
      CALL PLTPLY(NO,XWORK1,YWORK1)
      RETURN

      END
