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

C $Id: marksm.f,v 1.1 1990/11/30 11:11:51 gdsjaar Exp $
C $Log: marksm.f,v $
C Revision 1.1  1990/11/30 11:11:51  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]MARKSM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES, NODE,
     &   ERR)
C***********************************************************************
C
C  SUBROUTINE MARKSM = MARKS NODES WITHIN 2 LINE CONNECTIONS FROM NODE
C                      FOR SMOOTHING
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION LNODES (MLN, MXND)
      DIMENSION L1LIST(20), L2LIST(20)
C
      LOGICAL ERR
C
      IF (LXN (1, NODE) .LE. 0) GOTO 120
      CALL GETLXN (MXND, LXN, NODE, L1LIST, NL1, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** PROBLEMS IN MARKSM FINDING LXN **')
         GOTO 120
      ENDIF
C
      LNODES (4, NODE) = - IABS (LNODES (4, NODE))
      DO 110 I = 1, NL1
         NODE2 = NXL (1, L1LIST (I)) + NXL (2, L1LIST (I)) - NODE
         CALL GETLXN (MXND, LXN, NODE2, L2LIST, NL2, ERR)
         IF (ERR) THEN
            CALL MESAGE ('** PROBLEMS IN MARKSM FINDING LXN **')
            GOTO 120
         ENDIF
         LNODES (4, NODE2) = - IABS (LNODES (4, NODE2))
         DO 100 J = 1, NL2
            NODE1 = NXL (1, L2LIST (J)) + NXL (2, L2LIST (J)) - NODE2
            LNODES (4, NODE1) = - IABS (LNODES (4, NODE1))
  100    CONTINUE
  110 CONTINUE
C
  120 CONTINUE
      RETURN
C
      END
