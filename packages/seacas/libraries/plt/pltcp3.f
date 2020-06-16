C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltcp3.f,v 1.3 1993/07/19 17:06:34 gdsjaar Exp $
C $Log: pltcp3.f,v $
C Revision 1.3  1993/07/19 17:06:34  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:11  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:51  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTCP3(N,MASK,PX,PY,PZ,V,Q)
      DIMENSION MASK(*),PX(*),PY(*),PZ(*),V(*),Q(*)
      include 'izbit.inc'

      J = 0
      KM = 0
 2060 IF (.NOT. (J.LT.N)) GO TO 2070
      JN = MIN(N-J,32)
      KM = 1 + KM
      J1 = J
      J = J + JN
      M = MASK(KM)
      IF (M.EQ.0) THEN
         GO TO 2060

      END IF

      DO 2080 K = 1,JN
         JB = IZBIT(K)
         IF (IAND(M,JB).NE.0) THEN
            FP = (PX(J1+K)-V(1))*Q(1) + (PY(J1+K)-V(2))*Q(2) +
     *           (PZ(J1+K)-V(3))*Q(3)
            IF (FP.LT.0.) THEN
               M = IAND(M,NOT(JB))
            END IF

         END IF

 2080 CONTINUE
      MASK(KM) = M
      GO TO 2060

 2070 CONTINUE
      RETURN

      END
