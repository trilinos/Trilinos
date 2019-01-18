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

CC* FILE: [.QMESH]FNDLIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FNDLIN_FQ (MXND, LXN, NODE1, NODE2, LINE, ERR)
C
C***********************************************************************
C
C  SUBROUTINE FNDLIN =  FINDS THE LINE WITH ENDS NODE1 & NODE2
C
C***********************************************************************
C
      DIMENSION LXN(4, MXND)
      DIMENSION LINES1(20), LINES2(20)
      LOGICAL ERR
C
      ERR = .FALSE.
C
      CALL GETLXN (MXND, LXN, NODE1, LINES1, NL1, ERR)
      CALL GETLXN (MXND, LXN, NODE2, LINES2, NL2, ERR)
C
      IF (.NOT.ERR) THEN
         ERR = .TRUE.
         DO 110 I = 1, NL1
            DO 100 J = 1, NL2
               IF (LINES1(I) .EQ. LINES2(J)) THEN
                  LINE = LINES1(I)
                  ERR = .FALSE.
                  RETURN
               END IF
  100       CONTINUE
  110    CONTINUE
      END IF
C
      RETURN
C
      END
