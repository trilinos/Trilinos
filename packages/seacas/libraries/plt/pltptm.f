C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltptm.f,v 1.3 1993/07/19 17:06:42 gdsjaar Exp $
C $Log: pltptm.f,v $
C Revision 1.3  1993/07/19 17:06:42  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:20  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:49:09  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTPTM(N,MASK,X,Y)
      DIMENSION X(*),Y(*),MASK(*)
      include 'izbit.inc'

      J = 0
      J1 = J
      KM = 0
 2120 IF (.NOT. (J.LT.N)) GO TO 2130
      JN = MIN(N-J,32)
      KM = KM + 1
      M = MASK(KM)
      J = J + JN
      IF (M.EQ.0) THEN
         GO TO 2120

      END IF

      DO 2140 K = 1,JN
         IF (IAND(M,IZBIT(K)).NE.0) THEN
            CALL PLTPNT(1,X(J1+K),Y(J1+K))
         END IF

 2140 CONTINUE
      GO TO 2120

 2130 CONTINUE
      RETURN

      END
