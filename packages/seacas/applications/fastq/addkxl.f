C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDKXL (MXND, KXL, K, L)
C***********************************************************************

C  SUBROUTINE ADDKXL = ADDS TO THE LIST OF ELEMENTS FOR THIS LINE

C***********************************************************************

      DIMENSION KXL (2, 3*MXND)

      IF ( KXL(1, L) .EQ. 0) THEN
         KXL(1, L) = K
      ELSE
         KXL(2, L) = K
      ENDIF
      RETURN

      END
