C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SNAPXY (MP, MSNAP, N, IPOINT, COOR, LINKP, SNAPDX,
     &   NSNAP)
C***********************************************************************

C  SUBROUTINE SNAPXY = ADDS SNAP GRIDS AT EVERY X, Y LOCATION

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)
      DIMENSION SNAPDX (2, MSNAP), NSNAP (2)

      LOGICAL ERR, ADDLNK

      ADDLNK = .FALSE.

      DO 100 I = 1, N
         CALL LTSORT (MP, LINKP, IABS (IPOINT (I)), II, ADDLNK)
         IF (II.GT.0)THEN
            INDEX = 1
            CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, COOR (1, II), ERR)
            IF  (ERR) RETURN
            INDEX = 2
            CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, COOR (2, II), ERR)
            IF  (ERR) RETURN
         ENDIF
  100 CONTINUE

      RETURN

      END
