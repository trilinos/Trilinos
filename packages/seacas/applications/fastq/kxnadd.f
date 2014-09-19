C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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
