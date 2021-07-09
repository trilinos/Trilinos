C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE BUBBLE (X, KARRY, NORD, N)
C***********************************************************************

C  SUBROUTINE BUBBLE=SORTS ALL VALUES X(I), KARRY(I) INTO DECREASING
C                      ORDER, ASSUMING THAT VALUES 1 TO NORD ARE SORTED

C***********************************************************************

      DIMENSION X (N), KARRY (N)

      IF (N .LE. 1) RETURN

      ISTART = MAX0 (NORD + 1, 2)
      IF (ISTART .GT. N) RETURN
      DO 120 J = ISTART, N
         XVAL = X (J)
         KVAL = KARRY (J)
         JM1 = J - 1
         I = J
         DO 100 II = 1, JM1
            IF  (XVAL .LE. X (I - 1)) GO TO 110
            X (I) = X (I - 1)
            KARRY (I) = KARRY (I - 1)
            I = I - 1
  100    CONTINUE
  110    CONTINUE
         X (I) = XVAL
         KARRY (I) = KVAL
  120 CONTINUE

      RETURN

      END
