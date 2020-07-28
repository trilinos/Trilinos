C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DELLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NODE, LINE,
     &   NNN, ERR, NOROOM)
C***********************************************************************

C  SUBROUTINE DELLXN = DELETE LINE FROM THE LIST OF LINES FOR THIS NODE

C***********************************************************************

      DIMENSION LINES (20), LXN (4, MXND), NUID (MXND)

      LOGICAL ERR, NOROOM

      CALL GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
      IF (NL.LT.1) THEN
         WRITE (*, 10000)NODE
         GOTO 110
      ENDIF
      IF (ERR) GOTO 110

      K = 0
      DO 100 I = 1, NL
         IF (LINES (I) .NE. LINE) THEN
            K = K + 1
            LINES (K) = LINES (I)
         ENDIF
  100 CONTINUE

      IF (K .NE. NL - 1) THEN
         WRITE (*, 10010) NODE, (LINES (I), I = 1, NL)
         ERR = .TRUE.
         GOTO 110
      ENDIF
      NL = NL-1
      CALL PUTLXN (MXND, NL, LXN, NUID, NODE, LINES, NAVAIL, IAVAIL,
     &   NNN, ERR, NOROOM)

  110 CONTINUE
      RETURN

10000 FORMAT (' ERROR IN DELLXN - NODE', I5, ' HAS NO LINES')
10010 FORMAT (' ERROR IN DELLXN - NODE:', I5, /, ' LINES:', 20I5)

      END
