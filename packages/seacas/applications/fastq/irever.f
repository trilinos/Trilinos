C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE IREVER (L, N)
C***********************************************************************

C  SUBROUTINE IREVER = REVERS THE INTEGER ARRAY L (I), I=1, N

C***********************************************************************

      DIMENSION L (N)

      IF (N .LE. 1) RETURN
      NUP = N + 1
      M = N / 2
      DO 100 I = 1, M
         NUP = NUP - 1
         ITEMP = L (I)
         L (I) = L (NUP)
         L (NUP) = ITEMP
  100 CONTINUE

      RETURN

      END
