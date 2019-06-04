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

C $Id: addlxn.f,v 1.1 1990/11/30 11:02:57 gdsjaar Exp $
C $Log: addlxn.f,v $
C Revision 1.1  1990/11/30 11:02:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]ADDLXN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDLXN (MXND, LXN, NUID, NAVAIL, IAVAIL, NODE, LINE,
     &   NNN, ERR, NOROOM)
C***********************************************************************
C
C  SUBROUTINE ADDLXN = ADDS LINE TO THE LIST OF LINES FOR THIS NODE
C
C***********************************************************************
C
      DIMENSION LXN (4, MXND), NUID (MXND)
C
      LOGICAL ERR, NOROOM
C
      ERR = .FALSE.
      NN = NODE
  100 CONTINUE
C
C     LINK TO CONTINUATION
C
      IF (LXN (4, NN) .LT. 0) THEN
         NN = IABS (LXN (4, NN))
         GOTO 100
C
C  THERE IS ROOM FOR THE NEW LINE WITHOUT CONTINUING
C
      ELSEIF (LXN (4, NN) .EQ. 0) THEN
         DO 110 I = 1, 4
            IF  (LXN (I, NN) .EQ. 0) THEN
               LXN (I, NN) = LINE
               RETURN
            ENDIF
  110    CONTINUE
C
C  THIS CAN'T HAPPEN
C
         CALL MESAGE ('ERROR IN ADDLXN')
         ERR = .TRUE.
      ELSE
C
C  CREATE A CONTINUATION ENTRY
C
         IF (NAVAIL .LT. 1) THEN
            WRITE ( * , 10000)NODE
            ERR = .TRUE.
            NOROOM = .TRUE.
            RETURN
         ENDIF
C
         NEW = IAVAIL
         IF (NEW .GT. NNN)NNN = NEW
         IAVAIL = LXN (4, IAVAIL)
         NAVAIL = NAVAIL - 1
         LXN (1, NEW) =  - LXN (4, NN)
         LXN (2, NEW) = LINE
         IF (LXN (2, NN) .LT. 0)LXN (2, NEW) =  - LINE
         LXN (3, NEW) = 0
         LXN (4, NEW) = 0
         LXN (4, NN) =  - NEW
         NUID (NEW) = 0
      ENDIF
      RETURN
C
10000 FORMAT (' NODE TABLE OVERFLOW IN ADDLXN AT NODE', I5)
C
      END
