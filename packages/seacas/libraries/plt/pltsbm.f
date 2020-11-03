C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
