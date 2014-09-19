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

C $Id: nxkbdy.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: nxkbdy.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:12:58  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:12:57  gdsjaar
c Initial revision
c 
C
CC* FILE: [.RENUM]NXKBDY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NXKBDY (MDIM, NNXK, NPELEM, NXK, KKK, LIST, THREE,
     &   EIGHT, NINE)
C***********************************************************************
C
C  SUBROUTINE NXKBDY = FLAGS ALL SIDES OF ELEMENTS ONLY ONCE BY USE OF
C                      A HASH SCHEME
C
C***********************************************************************
C
C     NXK    = THE NODES PER ELEMENT ARRAY  (CONNECTIVITY)
C               (A NEGATIVE VALUE WILL INDICATE UNIQUENESS OF FOLLOWING
C              SIDE)
C***********************************************************************
C
      LOGICAL THREE, EIGHT, NINE, ITSOK
C
      DIMENSION NXK (NNXK, NPELEM), LIST (MDIM)
C
      DO 100 I = 1, MDIM
         LIST (I)  =  0
  100 CONTINUE
C
      DO 140 K = 1, KKK
         IF ((NXK (3, K) .EQ. 0) .AND. (THREE)) THEN
            ITSOK = .TRUE.
            NEND = 1
         ELSEIF  ((NXK (3, K) .NE. 0) .AND. ((EIGHT) .OR. (NINE))) THEN
            ITSOK = .TRUE.
            NEND = 4
         ELSE
            ITSOK = .FALSE.
         ENDIF
         IF (ITSOK) THEN
            DO 130 N = 1, NEND
C
C  COMPUTE HASH CODE FOR LINE
C
               N2 = N + 1
               IF (N .GE. 4)N2 = 1
               NODE1 = NXK (N, K)
               NODE2 = IABS (NXK (N2, K))
               IF (NODE2 .GT. NODE1) THEN
                  LVAL = NODE1 * 100000 + NODE2
               ELSE
                  LVAL = NODE2 * 100000 + NODE1
               ENDIF
C
C  CALCULATE THE BEGINNING HASH VALUE
C
               HOLD = FLOAT (NODE1 + NODE2) * 3.1830989
               LHASH = INT((HOLD-IFIX (HOLD)) * FLOAT (MDIM) + 1)
               LKEEP = LHASH
C
C  FIRST-TIMERS CLAIM THE NODE
C
  110          CONTINUE
               IF (LIST (LHASH).NE.0) THEN
                  IF (LIST (LHASH) .EQ. LVAL) GOTO 120
                  LHASH = LHASH + 1
                  IF (LHASH .EQ. MDIM) LHASH = 1
                  IF (LHASH.NE.LKEEP) GOTO 110
                  CALL MESAGE ('HASH SCHEME OVERFLOW IN KXNBDY')
                  STOP
               ENDIF
C
               LIST (LHASH) = LVAL
               NXK (N, K) = ISIGN (NXK (N, K), -1)
  120          CONTINUE
  130       CONTINUE
         ENDIF
  140 CONTINUE
      RETURN
      END
