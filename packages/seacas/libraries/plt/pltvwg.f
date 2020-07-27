C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTVWG(PLL,PUR,N,XV,YV,ZV,NO,XVO,YVO,ZVO)
      DIMENSION XV(*),YV(*),ZV(*),XVO(*),YVO(*),ZVO(*),PLL(*),PUR(*)
      DIMENSION XVT(50),YVT(50),QUR(2),QLL(2)

      QUR(1) = PUR(1) + .0001
      QUR(2) = PUR(2) + .0001
      QLL(1) = PLL(1) - .0001
      QLL(2) = PLL(2) - .0001
      NOSAVE = NO
      NT = 50
      CALL PLTLI1(QLL,QUR,N,XV,YV,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI2(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      NT = 50
      CALL PLTLI3(QLL,QUR,NO,XVO,YVO,NT,XVT,YVT)
      NO = NOSAVE
      CALL PLTLI4(QLL,QUR,NT,XVT,YVT,NO,XVO,YVO)
      DO 2180 I = 1,NO
         ZVO(I) = PLTPGZ(NO,XV,YV,ZV,XVO(I),YVO(I))
 2180 CONTINUE
      RETURN

      END
