C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMG2(MAP,N,XV,YV,NO,XVO,YVO)
      REAL MAP(*)
      INTEGER N
      REAL XV(*),YV(*)
      INTEGER NO
      REAL XVO(*),YVO(*)
      REAL XWORK(50),YWORK(50)
      INTEGER NWORK

      NOSAVE = NO
      AXX = MAP(1)
      AYY = MAP(4)
      AXY = MAP(3)
      AYX = MAP(2)
      BX = MAP(5)
      BY = MAP(6)
      DO 2220 I = 1,N
         XVO(I) = AXX*XV(I) + AXY*YV(I) + BX
         YVO(I) = AYX*XV(I) + AYY*YV(I) + BY
 2220 CONTINUE
      NWORK = 50
      CALL PLTCG2(N,XVO,YVO,NWORK,XWORK,YWORK,MAP(7),MAP(9))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(9),MAP(11))
      NWORK = 50
      CALL PLTCG2(NO,XVO,YVO,NWORK,XWORK,YWORK,MAP(11),MAP(13))
      NO = NOSAVE
      CALL PLTCG2(NWORK,XWORK,YWORK,NO,XVO,YVO,MAP(13),MAP(7))
      RETURN

      END
