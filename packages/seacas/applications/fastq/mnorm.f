C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MNORM (MXND, XN, YN, NXL, LLL, STDLEN)
C***********************************************************************

C  SUBROUTINE MNORM = FINDS THE AVERAGE LENGTH OF THOSE LINES NOT MUCH
C                     LONGER THAN THE AVERAGE

C***********************************************************************

      DIMENSION NXL (2, 3 * MXND), XN (MXND), YN (MXND)

      STDLEN = 0.
      NUML = 0
      S = 0.0
      DO 100 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            S = S + D
            NUML = NUML + 1
         ENDIF
  100 CONTINUE

      IF (NUML .LE. 0) RETURN
      TOL = 1.25 * S / DBLE(NUML)
      NUML = 0
      S = 0.0
      DO 110 L = 1, LLL
         N1 = NXL (1, L)
         N2 = NXL (2, L)
         IF (N1 .GT. 0) THEN
            D = SQRT ((XN (N1) - XN (N2)) **2 + (YN (N1) - YN (N2)) **2)
            IF (D .LT. TOL) THEN
               S = S + D
               NUML = NUML + 1
            ENDIF
         ENDIF
  110 CONTINUE
      STDLEN = S / DBLE(NUML)

      RETURN

      END
