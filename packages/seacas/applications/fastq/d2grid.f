C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE D2GRID (X1, Y1, X2, Y2)
C***********************************************************************

C  SUBROUTINE D2GRID = DRAWS A LINE BETWEEN TWO GRIDS

C***********************************************************************

      DIMENSION X (2), Y (2)

      X (1) = X1
      X (2) = X2
      Y (1) = Y1
      Y (2) = Y2
      CALL MPD2VC (1, X (1), Y (1), X (2), Y (2))
      RETURN

      END
