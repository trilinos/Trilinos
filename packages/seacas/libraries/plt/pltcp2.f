C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltcp2.f,v 1.3 1993/07/19 17:06:33 gdsjaar Exp $
C $Log: pltcp2.f,v $
C Revision 1.3  1993/07/19 17:06:33  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:10  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:50  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTCP2(N,MASK,PX,PY,C1,C2)
      DIMENSION MASK(*),PX(*),PY(*),C1(*),C2(*)
      include 'izbit.inc'

      CX = C1(1)
      CY = C1(2)
      DX = C2(1) - CX
      DY = C2(2) - CY
      J = 0
      KM = 0
 2020 IF (.NOT. (J.LT.N)) GO TO 2030
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2020

      END IF

      DO 2040 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).NE.0) THEN
            FP = (PY(J1+K)-CY)*DX - (PX(J1+K)-CX)*DY
            IF (FP.LT.0.) THEN
               M = IAND(M,NOT(JB))
            END IF

         END IF

 2040 CONTINUE
      MASK(KM) = M
      GO TO 2020

 2030 CONTINUE
      RETURN

      END
