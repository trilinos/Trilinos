C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      REAL FUNCTION ELIPSE (A7, A8, A2, ANG)
C***********************************************************************

C  FUNCTION ELIPSE = CALCULATES THE ANGULAR EQUATION ERROR WHEN FINDING
C                    AN ELIPSE PATTERN

C***********************************************************************

      PI = ATAN2(0.0, -1.0)
      A4 = A8 - ANG
      A5 = A7 - ANG
      A3 = A2 - A4
      A6 = PI - A5 - A2
      ELIPSE = SIN(A4) * SIN(A6) - SIN(A5) * SIN (A3)
      RETURN

      END
