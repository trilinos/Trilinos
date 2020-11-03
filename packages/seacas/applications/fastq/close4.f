C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CLOSE4 (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, N1, N2, N3, KKK, ERR)
C***********************************************************************

C  SUBROUTINE CLOSE4 = CLOSES THE AREA AROUND A FOUR NODE ELEMENT

C***********************************************************************

      DIMENSION LXK (4, MXND), KXL (2, 3*MXND), LNODES (MLN, MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)

      LOGICAL ERR

C  SET ALL THE LOOP NODES TO BE INTERIOR

      LNODES (4, N0) = - 2
      LNODES (4, N1) = - 2
      LNODES (4, N2) = - 2
      LNODES (4, N3) = - 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N1, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N2, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N3, ERR)
      IF (ERR) GOTO 100

C  MAKE THE LXK AND THE KXL ARRAY

      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N0)
      LXK (2, KKK) = LNODES (5, N1)
      LXK (3, KKK) = LNODES (5, N2)
      LXK (4, KKK) = LNODES (5, N3)

      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N2))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N3))

  100 CONTINUE
      RETURN

      END
