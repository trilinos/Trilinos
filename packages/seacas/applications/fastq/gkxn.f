C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GKXN (MXND, KXL, LXN, N, KS, KLIST, ERR)
C***********************************************************************

C  SUBROUTINE GKXN = GENERATES THE LIST OF ELEMENTS ASSOCIATED WITH
C                    NODE N

C***********************************************************************

      DIMENSION KLIST (1), KL (20), LINES (20)
      DIMENSION KXL (2, 3 * MXND), LXN (4, MXND)

      LOGICAL ERR

      ERR = .FALSE.
      KS = 0
      IF  (LXN (1, N) .LE. 0) RETURN
      CALL GETLXN (MXND, LXN, N, LINES, NL, ERR)
      IF (ERR) RETURN

C  LOOP THROUGH ALL LINES CONNECTED TO THIS NODE

      KOUNT = 0
      DO 140 IL = 1, NL
         L = LINES (IL)

C  LOOK AT ELEMENTS ON BOTH SIDES OF THIS LINE

         DO 130 IK = 1, 2
            K = KXL (IK, L)
            IF (K .GT. 0) THEN
               IF (KOUNT .GT. 0) THEN
                  DO 100 I = 1, KOUNT
                     IF  (K .EQ. KL (I)) GOTO 120
  100             CONTINUE
                  IF (KOUNT .GE. 20) THEN
                     ERR = .TRUE.
                     DO 110 I = 1, KOUNT
                        KLIST (I) = KL (I)
  110                CONTINUE
                     KS = KOUNT
                     RETURN
                  ENDIF
               ENDIF
               KOUNT = KOUNT + 1
               KL (KOUNT) = K
            ENDIF
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE

C  RETURN RESULTS

      DO 150 I = 1, KOUNT
         KLIST (I) = KL (I)
  150 CONTINUE
      KS = KOUNT

      RETURN

      END
