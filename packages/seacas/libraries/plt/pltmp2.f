C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTMP2(UMP,N,MASK,PX,PY,QX,QY)
      DIMENSION UMP(*),MASK(*),PX(*),PY(*),QX(*),QY(*)

      AXX = UMP(1)
      AYY = UMP(4)
      AXY = UMP(3)
      AYX = UMP(2)
      BX = UMP(5)
      BY = UMP(6)
      DO 2040 I = 1,N
         QX(I) = AXX*PX(I) + AXY*PY(I) + BX
         QY(I) = AYX*PX(I) + AYY*PY(I) + BY
 2040 CONTINUE
      CALL PLTCP2(N,MASK,QX,QY,UMP(7),UMP(9))
      CALL PLTCP2(N,MASK,QX,QY,UMP(9),UMP(11))
      CALL PLTCP2(N,MASK,QX,QY,UMP(11),UMP(13))
      CALL PLTCP2(N,MASK,QX,QY,UMP(13),UMP(7))
      RETURN

      END
