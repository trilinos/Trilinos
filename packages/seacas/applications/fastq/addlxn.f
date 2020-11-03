C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NODE, LINE,
     &   NNN, ERR, NOROOM)
C***********************************************************************

C  SUBROUTINE ADDLXN = ADDS LINE TO THE LIST OF LINES FOR THIS NODE

C***********************************************************************

      DIMENSION LXN (4, MXND), NUID (MXND)

      LOGICAL ERR, NOROOM

      ERR = .FALSE.
      NN = NODE
  100 CONTINUE

C     LINK TO CONTINUATION

      IF (LXN (4, NN) .LT. 0) THEN
         NN = IABS (LXN (4, NN))
         GOTO 100

C  THERE IS ROOM FOR THE NEW LINE WITHOUT CONTINUING

      ELSEIF (LXN (4, NN) .EQ. 0) THEN
         DO 110 I = 1, 4
            IF  (LXN (I, NN) .EQ. 0) THEN
               LXN (I, NN) = LINE
               RETURN
            ENDIF
  110    CONTINUE

C  THIS CAN'T HAPPEN

         CALL MESAGE ('ERROR IN ADDLXN')
         ERR = .TRUE.
      ELSE

C  CREATE A CONTINUATION ENTRY

         IF (NAVAIL .LT. 1) THEN
            WRITE ( * , 10000)NODE
            ERR = .TRUE.
            NOROOM = .TRUE.
            RETURN
         ENDIF

         NEW = IAVAIL
         IF (NEW .GT. NNN)NNN = NEW
         IAVAIL = LXN (4, IAVAIL)
         NAVAIL = NAVAIL - 1
         LXN (1, NEW) =  - LXN (4, NN)
         LXN (2, NEW) = LINE
         IF (LXN (2, NN) .LT. 0)LXN (2, NEW) =  - LINE
         LXN (3, NEW) = 0
         LXN (4, NEW) = 0
         LXN (4, NN) =  - NEW
         NUID (NEW) = 0
      ENDIF
      RETURN

10000 FORMAT (' NODE TABLE OVERFLOW IN ADDLXN AT NODE', I5)

      END
