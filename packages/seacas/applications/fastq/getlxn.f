C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
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
