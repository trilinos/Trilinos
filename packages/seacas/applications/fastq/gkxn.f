C $Id: gkxn.f,v 1.1 1990/11/30 11:08:55 gdsjaar Exp $
C $Log: gkxn.f,v $
C Revision 1.1  1990/11/30 11:08:55  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]GKXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GKXN (MXND, KXL, LXN, N, KS, KLIST, ERR)
C***********************************************************************
C
C  SUBROUTINE GKXN = GENERATES THE LIST OF ELEMENTS ASSOCIATED WITH
C                    NODE N
C
C***********************************************************************
C
      DIMENSION KLIST (1), KL (20), LINES (20)
      DIMENSION KXL (2, 3 * MXND), LXN (4, MXND)
C
      LOGICAL ERR
C
      ERR = .FALSE.
      KS = 0
      IF  (LXN (1, N) .LE. 0) RETURN
      CALL GETLXN (MXND, LXN, N, LINES, NL, ERR)
      IF (ERR) RETURN
C
C  LOOP THRU ALL LINES CONNECTED TO THIS NODE
C
      KOUNT = 0
      DO 140 IL = 1, NL
         L = LINES (IL)
C
C  LOOK AT ELEMENTS ON BOTH SIDES OF THIS LINE
C
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
C
C  RETURN RESULTS
C
      DO 150 I = 1, KOUNT
         KLIST (I) = KL (I)
  150 CONTINUE
      KS = KOUNT
C
      RETURN
C
      END
