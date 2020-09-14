C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE KXNADD (MAXKXN, NNXK, KXN, NUMKXN, K, NODE, ERR)
C************************************************************************

C  SUBROUTINE KXNADD = ADDS K AS AN ELEMENT OF NODE

C***********************************************************************

C  NOTE:
C     IT IS ASSUMED K IS NOT ALREADY AN ELEMENT OF NODE

C***********************************************************************

      DIMENSION KXN (NNXK, MAXKXN)

      LOGICAL ERR

      ERR = .FALSE.
      NN = NODE
  100 CONTINUE

C  LINE CONTINUES  -  FIND NEW CONTINUATION LINE

      IF (KXN (4, NN) .LT. 0) THEN
         NN = IABS (KXN (4, NN))
         GOTO 100

C  ADD THE ELEMENT TO NODE

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

C  ADD A CONTINUATION LINE,  AND ADD THE ELEMENT TO NODE

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

10000 FORMAT ('FOR ELEMENT', I5, ',  AND NODE', I5)

      END
