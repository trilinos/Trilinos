C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE CHKKXL (MXND, LXK, KXL, LLL, ERR)
C***********************************************************************

C  SUBROUTINE CHKKXL = CHECKS TO SEE IF THE KXL COMPARES CORRECTLY TO
C                      THE LXK ARRAY

C***********************************************************************

      DIMENSION LXK (4, MXND), KXL (2, 3 * MXND)

      LOGICAL ERR

      ERR = .TRUE.

      DO 130 L = 1, LLL
         DO 120 IK = 1, 2
            K = KXL (IK, L)
            IF (K .NE. 0) THEN
               DO 100 I = 1, 4
                  IF (LXK (I, K) .EQ. L)GOTO 110
  100          CONTINUE
               WRITE ( * , 10000)IK, L, K
               RETURN
            ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
      ERR = .FALSE.

      RETURN

10000 FORMAT ('KXL(', I4, ',', I4,') = ', I4,
     &   ' IS NOT IN LXK ARRAY  -  CHKKXL')

      END
