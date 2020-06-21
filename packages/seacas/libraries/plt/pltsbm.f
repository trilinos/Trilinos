C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltsbm.f,v 1.3 1993/07/19 17:06:43 gdsjaar Exp $
C $Log: pltsbm.f,v $
C Revision 1.3  1993/07/19 17:06:43  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:21  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:49:28  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTSBM(N,MASK,X,Y,SYMB)
      DIMENSION X(*),Y(*),MASK(*)
      CHARACTER*(*) SYMB
      include 'izbit.inc'

      J = 0
      KM = 0
 2170 IF (.NOT. (J.LT.N)) GO TO 2180
      JN = MIN(N-J,32)
      KM = KM + 1
      M = MASK(KM)
      J1 = J
      J = J + JN
      IF (M.EQ.0) THEN
         GO TO 2170

      END IF

      DO 2190 K = 1,JN
         IF (IAND(M,IZBIT(K)).NE.0) THEN
            CALL PLTXTS(X(J1+K),Y(J1+K),SYMB)
         END IF

 2190 CONTINUE
      GO TO 2170

 2180 CONTINUE
      RETURN

      END
