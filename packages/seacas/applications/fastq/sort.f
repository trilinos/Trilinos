C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SORT (N, IX, IY)
C***********************************************************************

C  SUBROUTINE SORT = SORT THE ARRAY IX,  CARRYING ALONG IY

C***********************************************************************

      DIMENSION IX (N), IY (N)
      NN = N
      M = NN
  100 CONTINUE
      M =  (9 * M) / 16
      IF (M .LE. 0) RETURN
      M1 = M + 1
      DO 120 J = M1, NN
         L = J
         I = J - M
  110    CONTINUE
         IF (IX (L) .LT. IX (I)) THEN
            KEEPX = IX (I)
            KEEPY = IY (I)
            IX (I) = IX (L)
            IY (I) = IY (L)
            IX (L) = KEEPX
            IY (L) = KEEPY
            L = I
            I = I - M
            IF (I .GE. 1)GOTO 110
         ENDIF
  120 CONTINUE

      GOTO 100

      END
