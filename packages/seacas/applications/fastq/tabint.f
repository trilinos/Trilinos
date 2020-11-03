C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE TABINT (X1, X2, Y1, Y2, CT, ST, SCALE, XX1, YY1, XX2,
     &   YY2, DRWTAB)
C***********************************************************************

C     SUBROUTINE TABINT = INITIALIZES THE TABLET TO THE PLOT LIMITS

C***********************************************************************

      LOGICAL DRWTAB

      IF (DRWTAB) THEN
         THETA = ATAN2 (YY2 - YY1, XX2 - XX1) -
     &      ATAN2 (Y2 - Y1, X2 - X1)
         CT = COS (THETA)
         ST = SIN (THETA)
         SCALE = SQRT (((X2 - X1) ** 2 + (Y2 - Y1) ** 2 ) /
     &      ((XX2 - XX1) ** 2 + (YY2 - YY1) ** 2 ))
      ELSE
         CT = 1.
         ST = 0.
         XX1 = 2000
         XX2 = 15000
         YY1 = 2000
         YY2 = 10000
         SCALEX =  (X2 - X1) / (XX2 - XX1)
         SCALEY =  (Y2 - Y1) / (YY2 - YY1)
         IF (SCALEX .GT. SCALEY) THEN
            SCALE = SCALEX
            YY1 =  (YY2 - YY1) -  ( (Y2 - Y1) / SCALE)
         ELSE
            SCALE = SCALEY
            XX1 =  (XX2 - XX1) -  ( (X2 - X1) / SCALE)
         ENDIF
      ENDIF

      RETURN

      END
