C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE MINMAX_FQ (NDIM, N, X, Y, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************

C  SUBROUTINE MINMAX_FQ = COMPARES AND RECORDS X AND Y EXTREMES

C***********************************************************************

      DIMENSION X (NDIM), Y (NDIM)
      XMIN = X (1)
      XMAX = X (1)
      YMIN = Y (1)
      YMAX = Y (1)

      DO 100 I = 1, N
         XMIN = AMIN1 (X (I), XMIN)
         XMAX = AMAX1 (X (I), XMAX)
         YMIN = AMIN1 (Y (I), YMIN)
         YMAX = AMAX1 (Y (I), YMAX)
  100 CONTINUE

      RETURN

      END
