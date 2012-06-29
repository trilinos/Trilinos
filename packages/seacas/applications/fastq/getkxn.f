C $Id: getkxn.f,v 1.1 1990/11/30 11:08:16 gdsjaar Exp $
C $Log: getkxn.f,v $
C Revision 1.1  1990/11/30 11:08:16  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]GETKXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETKXN (NPNODE, MAXKXN, NNXK, KXN, NUID, NODE, KLIST,
     &   NUMK, ERR)
C***********************************************************************
C
C  SUBROUTINE GETKXN = GET THE LIST OF ELEMENTS RELATED TO THIS NODE
C
C***********************************************************************
C
      DIMENSION KLIST (20), NUID (NPNODE), KXN (NNXK, MAXKXN)
C
      LOGICAL ERR
C
      ERR = .FALSE.
      NUM = 0
      NN = NODE
C
C  ADD IN THE FIRST THREE NODES LISTED
C
  100 CONTINUE
      DO 110 I = 1, 3
         IF (KXN (I, NN) .EQ. 0) THEN
            NUMK = NUM
            IF (NUMK .GE. 1) THEN
               RETURN
            ELSE
               WRITE (*, 10000)NODE, NUID (NODE)
               ERR = .TRUE.
               RETURN
            ENDIF
         ENDIF
         NUM = NUM + 1
         KLIST (NUM) = KXN (I, NN)
  110 CONTINUE
C
C  CHECK THE FOURTH NODE FOR CONTINUATION
C
      IF (KXN (4, NN) .LT. 0) THEN
         NN = IABS (KXN (4, NN))
         IF (NUM .LT. 18) THEN
            GOTO 100
         ELSE
            WRITE (*, 10010)NODE, NUID (NODE)
            ERR = .TRUE.
            RETURN
         ENDIF
C
C  ADD IN THE LAST NODE IF IT IS NONZERO
C
      ELSE
         IF (KXN (4, NN) .NE. 0) THEN
            NUM = NUM + 1
            KLIST (NUM) = KXN (4, NN)
         ENDIF
         NUMK = NUM
         IF (NUMK .GE. 1) THEN
            RETURN
         ELSE
            WRITE (*, 10000)NODE, NUID (NODE)
            ERR = .TRUE.
            RETURN
         ENDIF
      ENDIF
C
10000 FORMAT (' NO ELEMENTS CONNECTED TO NODE', I5, ', NUID  = ', I10)
10010 FORMAT (' TOO MANY ELEMENTS CONNECTED TO NODE', I5, ', NUID  = ',
     &   I10)
C
      END
