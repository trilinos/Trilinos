C $Id: addlxn.f,v 1.1 1990/11/30 11:02:57 gdsjaar Exp $
C $Log: addlxn.f,v $
C Revision 1.1  1990/11/30 11:02:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]ADDLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NODE, LINE,
     &   NNN, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE ADDLXN = ADDS LINE TO THE LIST OF LINES FOR THIS NODE
C
C***********************************************************************
C
      DIMENSION LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, NOROOM
C
      ERR = .FALSE.
      NN = NODE
  100 CONTINUE
C
C     LINK TO CONTINUATION
C
      IF (LXN (4, NN) .LT. 0) THEN
         NN = IABS (LXN (4, NN))
         GOTO 100
C
C  THERE IS ROOM FOR THE NEW LINE WITHOUT CONTINUING
C
      ELSEIF (LXN (4, NN) .EQ. 0) THEN
         DO 110 I = 1, 4
            IF  (LXN (I, NN) .EQ. 0) THEN
               LXN (I, NN) = LINE
               RETURN
            ENDIF
  110    CONTINUE
C
C  THIS CAN'T HAPPEN
C
         CALL MESAGE ('ERROR IN ADDLXN')
         ERR = .TRUE.
      ELSE
C
C  CREATE A CONTINUATION ENTRY
C
         IF (NAVAIL .LT. 1) THEN
            WRITE ( * , 10000)NODE
            ERR = .TRUE.
            NOROOM = .TRUE.
            RETURN
         ENDIF
C
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
C
10000 FORMAT (' NODE TABLE OVERFLOW IN ADDLXN AT NODE', I5)
C
      END
