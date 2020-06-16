C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltdp2.f,v 1.2 2000/10/25 18:55:02 gdsjaar Exp $
C $Log: pltdp2.f,v $
C Revision 1.2  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.1  1993/07/16 16:48:01  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTDP2(MAP,N,PX,PY)
      REAL MAP(*),PX(*),PY(*)
      REAL XWORK(32),YWORK(32)
      INTEGER MASK(1)

      J = 0
 2360 IF (.NOT. (J.LT.N)) GO TO 2370
      JN = MIN(N-J,32)
      J1 = J
      J = J + JN

      do 2400 i=1, jn
        XWORK(I) = MAP(1)*PX(J1+I) + MAP(3)*PY(J1+I) + MAP(5)
        YWORK(I) = MAP(2)*PX(J1+I) + MAP(4)*PY(J1+I) + MAP(6)
 2400 CONTINUE

      MASK(1) = -1
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(7),MAP(9))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(9),MAP(11))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(11),MAP(13))
      CALL PLTCP2(JN,MASK,XWORK,YWORK,MAP(13),MAP(7))
      CALL PLTPTM(JN,MASK,XWORK,YWORK)
      GO TO 2360

 2370 CONTINUE
      RETURN

      END
