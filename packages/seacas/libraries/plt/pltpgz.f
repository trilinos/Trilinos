C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      FUNCTION PLTPGZ(N,X,Y,Z,XQ,YQ)
      DIMENSION X(3),Y(3),Z(3),D(2)

      X21 = X(2) - X(1)
      X31 = X(3) - X(1)
      Y21 = Y(2) - Y(1)
      Y31 = Y(3) - Y(1)
      Z21 = Z(2) - Z(1)
      Z31 = Z(3) - Z(1)
      A = Y31*Z21 - Z31*Y21
      B = Z31*X21 - X31*Z21
      C = X31*Y21 - Y31*X21
      IF (C.EQ.0.) THEN
         DO 2200 I = 1,N - 2
            D(1) = X(I+1) - X(I)
            D(2) = Y(I+1) - Y(I)
            DDD = D(1)*D(1) + D(2)*D(2)
            IF (DDD.NE.0.) THEN
               GO TO 2210

            END IF

 2200    CONTINUE
 2210    CONTINUE
         IF (DDD.EQ.0.) THEN
            PLTPGZ = 0.0
            RETURN

         END IF

         DDP = D(1)*XQ + D(2)*YQ
         DDP1 = D(1)*X(1) + D(2)*Y(1)
         ALPHA = (DDP-DDP1)/DDD
         PLTPGZ = Z(1) + ALPHA* (Z(2)-Z(1))
         RETURN

      END IF

      PLTPGZ = ((A* (X(1)-XQ)+B* (Y(1)-YQ))/C) + Z(1)
      RETURN

      END
