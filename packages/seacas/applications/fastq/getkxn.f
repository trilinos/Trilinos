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
