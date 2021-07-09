C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PENDIS (SCORD, FCORD, DIST, NDIM, NNOD)
      DIMENSION SCORD(*), FCORD(NDIM, NNOD)

      IF (NDIM .EQ. 3) THEN

C -- DETERMINE PLANE EQUATION

         XI = FCORD(1, 1)
         YI = FCORD(2, 1)
         ZI = FCORD(3, 1)

         XJ = FCORD(1, 2)
         YJ = FCORD(2, 2)
         ZJ = FCORD(3, 2)

         XK = FCORD(1, 3)
         YK = FCORD(2, 3)
         ZK = FCORD(3, 3)

         XL = FCORD(1, 4)
         YL = FCORD(2, 4)
         ZL = FCORD(3, 4)

         A =  (YK - YI) * (ZL - ZJ) - (ZK - ZI) * (YL - YJ)
         B =  (ZK - ZI) * (XL - XJ) - (XK - XI) * (ZL - ZJ)
         C =  (XK - XI) * (YL - YJ) - (YK - YI) * (XL - XJ)
         RMAG = SQRT (A**2 + B**2 + C**2)

         A = A / RMAG
         B = B / RMAG
         C = C / RMAG
         D = A * FCORD(1,1) + B * FCORD(2,1) + C * FCORD(3,1)

         DIST = ABS(A * SCORD(1) + B * SCORD(2) + C * SCORD(3) - D) /
     *      SQRT(A**2 + B**2 + C**2)

      ELSE IF (NDIM .EQ. 2) THEN
         A  = FCORD(1,2) - FCORD(1,1)
         B  = FCORD(2,2) - FCORD(2,1)

         X1 = FCORD(1,1)
         Y1 = FCORD(2,1)

         X0 = SCORD(1)
         Y0 = SCORD(2)
         T  = -1. * (A * (X1 - X0) + B * (Y1 - Y0)) / (A**2 + B**2)

         X = X1 + A * T
         Y = Y1 + B * T

         DIST = SQRT((X - X0)**2 + (Y - Y0)**2)
      END IF
      RETURN
      END
