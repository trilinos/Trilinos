C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FQ_ROTATE (N, X, Y, NID, NEWF)
C***********************************************************************

C  SUBROUTINE ROTATE = CIRCULARLY SHIFTS THE DATA IN X,  Y,  AND NID

C***********************************************************************

      DIMENSION X (N), Y (N), NID (N)

      IF ((NEWF .LE. 1) .OR. (NEWF .GT. N)) RETURN

C  BUBBLE UP THROUGH THE ARRAYS AS MANY TIMES AS NEEDED

      DO 110 I = 1, NEWF - 1
         XLAST = X (1)
         YLAST = Y (1)
         NLAST = NID (1)
         DO 100 J = 1, N - 1
            X(J) = X (J + 1)
            Y(J) = Y (J + 1)
            NID(J) = NID (J + 1)
  100    CONTINUE
         X(N)   = XLAST
         Y(N)   = YLAST
         NID(N) = NLAST
  110 CONTINUE

      RETURN

      END
