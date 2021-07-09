C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE D2NODE (MXND, XN, YN, NODE1, NODE2)
C***********************************************************************

C  SUBROUTINE D2NODE = DRAWS A LINE BETWEEN TWO NODES

C***********************************************************************

      DIMENSION X (2), Y (2), XN (MXND), YN (MXND)

      X (1) = XN (NODE1)
      X (2) = XN (NODE2)
      Y (1) = YN (NODE1)
      Y (2) = YN (NODE2)
      CALL MPD2VC (1, X (1), Y (1), X (2), Y (2))
      RETURN

      END
