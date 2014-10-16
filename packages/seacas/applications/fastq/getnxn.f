C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

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
