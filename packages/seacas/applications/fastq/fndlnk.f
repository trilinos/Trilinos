C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FNDLNK (MXND, LXK, NXL, K, N1, N2, L, ERR)
C***********************************************************************

C  SUBROUTINE FNDLNK = FIND THE LINE IN ELEMENT K WITH NODES N1 AND N2

C***********************************************************************

      DIMENSION LXK (4, MXND), NXL (2, 3 * MXND)

      LOGICAL ERR

      ERR = .FALSE.
      DO 100 I = 1, 4
         LL = LXK (I, K)
         M1 = NXL (1, LL)
         M2 = NXL (2, LL)
         IF ( ( (M1 .EQ. N1) .AND. (M2 .EQ. N2)) .OR.
     &      ( (M2 .EQ. N1) .AND. (M1 .EQ. N2) ) ) THEN
            L = LL
            RETURN
         ENDIF
  100 CONTINUE
      L = 0
      ERR = .TRUE.
      WRITE ( * , 10000) K, N1, N2
10000 FORMAT (' IN FNDLNK, NO LINE CAN BE FOUND FOR K, N1, N2: ', 3I5)
      RETURN
      END
