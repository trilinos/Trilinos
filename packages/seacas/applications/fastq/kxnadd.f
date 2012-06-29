C $Id: kxnadd.f,v 1.1 1990/11/30 11:10:51 gdsjaar Exp $
C $Log: kxnadd.f,v $
C Revision 1.1  1990/11/30 11:10:51  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]KXNADD.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE KXNADD (MAXKXN, NNXK, KXN, NUMKXN, K, NODE, ERR)
C************************************************************************
C
C  SUBROUTINE KXNADD = ADDS K AS AN ELEMENT OF NODE
C
C***********************************************************************
C
C  NOTE:
C     IT IS ASSUMED K IS NOT ALREADY AN ELEMENT OF NODE
C
C***********************************************************************
C
      DIMENSION KXN (NNXK, MAXKXN)
C
      LOGICAL ERR
C
      ERR = .FALSE.
      NN = NODE
  100 CONTINUE
C
C  LINE CONTINUES  -  FIND NEW CONTINUATION LINE
C
      IF (KXN (4, NN) .LT. 0) THEN
         NN = IABS (KXN (4, NN))
         GOTO 100
C
C  ADD THE ELEMENT TO NODE
C
      ELSEIF (KXN (4, NN) .EQ. 0) THEN
         DO 110 I = 1, 4
            IF (KXN (I, NN) .EQ. 0) THEN
               KXN (I, NN) = K
               RETURN
            ENDIF
  110    CONTINUE
         CALL MESAGE ('IMPOSSIBLE SITUATION IN KXNADD')
         WRITE ( * , 10000)K, NODE
         ERR = .TRUE.
         RETURN
C
C  ADD A CONTINUATION LINE,  AND ADD THE ELEMENT TO NODE
C
      ELSE
         IF (NUMKXN .GE. MAXKXN) THEN
            CALL MESAGE ('NO ROOM FOR KXN TABLE IN KXNADD')
            ERR = .TRUE.
            RETURN
         ENDIF
         NUMKXN = NUMKXN + 1
         KXN (1, NUMKXN) = KXN (4, NN)
         KXN (2, NUMKXN) = K
         KXN (3, NUMKXN) = 0
         KXN (4, NUMKXN) = 0
         KXN (4, NN) =  - NUMKXN
         RETURN
      ENDIF
C
10000 FORMAT ('FOR ELEMENT', I5, ',  AND NODE', I5)
C
      END
