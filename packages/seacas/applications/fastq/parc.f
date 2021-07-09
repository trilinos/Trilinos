C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      FUNCTION PARC (AL, TCOEF)
C***********************************************************************

C    SUBROUTINE PARC = CALCULATES PARABOLIC ARC LOCATIONS

C***********************************************************************

      PARC = 0.5 * (SQRT (1.0 + AL **2) * AL +
     &   ALOG (SQRT (1.0 + AL **2) + AL)) / TCOEF

      END
