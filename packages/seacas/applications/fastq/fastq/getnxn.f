C $Id: getnxn.f,v 1.1 1990/11/30 11:08:31 gdsjaar Exp $
C $Log: getnxn.f,v $
C Revision 1.1  1990/11/30 11:08:31  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]GETNXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETNXN (NPNODE, NPELEM, MAXKXN, NNXK, KXN, NXK, NUID,
     &   NODE, NLIST, NUMN, ALL, ERR)
C***********************************************************************
C
C  SUBROUTINE GETNXN = GETS THE LIST OF NODES CONNECTED TO NODE
C
C***********************************************************************
C
C  NOTE:
C     NODES FOR WHICH NUID (NODE) IS NEGATIVE WILL NOT BE INCLUDED.
C     IF ALL=.FALSE.,  ONLY DIRECTLY CONNECTED NODES WILL BE LISTED.
C     IF ALL=.TRUE.,  INDIRECTLY CONNECTED NODES WILL ALSO BE LISTED.
C
C***********************************************************************
C
      DIMENSION NLIST (20), KLIST (20), M (3)
      DIMENSION KXN (NNXK, MAXKXN), NUID (NPNODE)
      DIMENSION NXK (NNXK, NPELEM)
C
      LOGICAL ALL, ERR
C
      ERR = .FALSE.
      CALL GETKXN (NPNODE, MAXKXN, NNXK, KXN, NUID, NODE, KLIST, NUMK,
     &   ERR)
      IF (ERR) RETURN
      IF (ALL) THEN
         NDO = 3
      ELSE
         NDO = 2
      ENDIF
      NUM = 0
      NOD = NODE
C
      DO 130 IK = 1, NUMK
         K = KLIST (IK)
         IF (NXK (1, K) .EQ. NOD) THEN
            M (1) = 4
            M (2) = 2
            M (3) = 3
         ELSEIF (NXK (2, K) .EQ. NOD) THEN
            M (1) = 1
            M (2) = 3
            M (3) = 4
         ELSEIF (NXK (3, K) .EQ. NOD) THEN
            M (1) = 2
            M (2) = 4
            M (3) = 1
         ELSEIF (NXK (4, K) .EQ. NOD) THEN
            M (1) = 3
            M (2) = 1
            M (3) = 2
         ELSE
            CALL MESAGE ('IMPOSSIBLE SITUATION IN GETNXN, LOOP 50')
            ERR = .TRUE.
            RETURN
         ENDIF
C
         NLK = NUM
         DO 120 IDO = 1, NDO
            MIDO = M (IDO)
            N = NXK (MIDO, K)
            IF ( (N .GT. 0) .AND. (NUID (N) .GT. 0)) THEN
               IF (NLK .LE. 0) THEN
                  NUM = NUM + 1
                  NLIST (NUM) = N
               ELSE
                  DO 100 I = 1, NLK
                     IF (NLIST (I) .EQ. N)GOTO 110
  100             CONTINUE
                  NUM = NUM + 1
                  NLIST (NUM) = N
               ENDIF
            ENDIF
  110       CONTINUE
  120    CONTINUE
  130 CONTINUE
C
      NUMN = NUM
      IF (NUMN .GT. 20) THEN
         WRITE (*, 10000)NODE, NUID (NODE)
         ERR = .TRUE.
      ENDIF
C
      RETURN
C
10000 FORMAT  (' TOO MANY NODES CONNECTED TO NODE', I5,
     &   ', NUID  = ', I10)
C
      END
