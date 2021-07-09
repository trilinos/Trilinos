C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE INGRID (MSNAP, SNAPDX, NSNAP, II, RIN, IFOUND, ERR)
C***********************************************************************

C  SUBROUTINE INGRID = INPUTS A X OR Y GRID CARD

C***********************************************************************

      DIMENSION SNAPDX(2, MSNAP), NSNAP(2), RIN(IFOUND)
      LOGICAL ERR

      ERR = .FALSE.
      IF (II .LT. 1 .OR. II .GT. 2) THEN
         ERR = .TRUE.
         WRITE (*, 10000) II
         GO TO 110
      END IF

      DO 100 I = 1, IFOUND
         CALL ADDSNP (MSNAP, SNAPDX, NSNAP, II, RIN(I), ERR)
         IF (ERR) GO TO 110
  100 CONTINUE

  110 CONTINUE
      RETURN

10000 FORMAT (' GRID INDEX OUT-OF-RANGE [1,2] IN INGRID: ', I3)
      END
