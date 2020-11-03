C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE REVERS (X, N)
C***********************************************************************

C  SUBROUTINE REVERS = REVERS THE REAL ARRAY OF X (I), I=1, N

C***********************************************************************

      DIMENSION X (N)

      IF (N .LE. 1) RETURN

      NUP = N + 1
      M = N / 2
      DO 100 I = 1, M
         NUP = NUP - 1
         XK = X (I)
         X (I) = X (NUP)
         X (NUP) = XK
  100 CONTINUE

      RETURN

      END
