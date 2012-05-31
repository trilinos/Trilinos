C $Id: getlxn.f,v 1.1 1990/11/30 11:08:19 gdsjaar Exp $
C $Log: getlxn.f,v $
C Revision 1.1  1990/11/30 11:08:19  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]GETLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
C***********************************************************************
C
C  SUBROUTINE GETLXN = GET THE FULL LIST OF LINES CONNECTED TO NODE
C
C***********************************************************************
C
      DIMENSION LINES (20), LXN (4, MXND)
C
      LOGICAL ERR
C
      ERR = .FALSE.
      NN = NODE
      NUM = 0
      IF (LXN (1, NN) .LE. 0) THEN
         NL = 0
         ERR = .TRUE.
         RETURN
      ENDIF
  100 LINES (NUM + 1) = IABS (LXN (1, NN))
      NUM = NUM + 2
      LINES (NUM) = IABS (LXN (2, NN))
      L = LXN (3, NN)
      IF (L.EQ.0) THEN
         NL = NUM
         RETURN
      ENDIF
      NUM = NUM + 1
      LINES (NUM) = IABS (L)
      L = LXN (4, NN)
      IF (L.LT.0) THEN
         NN = -L
         IF (NUM .LT. 18) THEN
            GOTO 100
         ELSE
            WRITE (*, 10000)NODE
            ERR = .TRUE.
            RETURN
         ENDIF
      ELSEIF (L .EQ. 0) THEN
         NL = NUM
         RETURN
      ELSE
         NUM = NUM + 1
         LINES (NUM) = L
         NL = NUM
         RETURN
      ENDIF
C
10000 FORMAT (' IN GETLXN, TOO MANY NODES CONNECTED TO NODE', I5)
C
      END
