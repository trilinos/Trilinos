C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETKXN (NPNODE, MAXKXN, NNXK, KXN, NUID, NODE, KLIST,
     &   NUMK, ERR)
C***********************************************************************

C  SUBROUTINE GETKXN = GET THE LIST OF ELEMENTS RELATED TO THIS NODE

C***********************************************************************

      DIMENSION KLIST (20), NUID (NPNODE), KXN (NNXK, MAXKXN)

      LOGICAL ERR

      ERR = .FALSE.
      NUM = 0
      NN = NODE

C  ADD IN THE FIRST THREE NODES LISTED

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

C  CHECK THE FOURTH NODE FOR CONTINUATION

      IF (KXN (4, NN) .LT. 0) THEN
         NN = IABS (KXN (4, NN))
         IF (NUM .LT. 18) THEN
            GOTO 100
         ELSE
            WRITE (*, 10010)NODE, NUID (NODE)
            ERR = .TRUE.
            RETURN
         ENDIF

C  ADD IN THE LAST NODE IF IT IS NONZERO

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

10000 FORMAT (' NO ELEMENTS CONNECTED TO NODE', I5, ', NUID  = ', I10)
10010 FORMAT (' TOO MANY ELEMENTS CONNECTED TO NODE', I5, ', NUID  = ',
     &   I10)

      END
