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

C $Id: setn02.f,v 1.1 1990/11/30 11:15:30 gdsjaar Exp $
C $Log: setn02.f,v $
C Revision 1.1  1990/11/30 11:15:30  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SETN02.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SETN02 (MXND, NXL, LXK, KXL, LINE, NEND, NODE, N0, N2)
C***********************************************************************
C
C  SUBROUTINE SETN02 = PICKS THE NEXT LINE AROUND THE ELEMENTS ATTACHED
C                      TO LINE WITH ONE END AT NEND, AND THE OTHER END
C                      NOT AT NODE, AND FROM THE CONNECTIVITY OF THE
C                      ELEMENTS DETERMINES THE BOUNDING ANGULAR LINES
C                      AND NODES.
C
C***********************************************************************
C
      DIMENSION NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
C
      K1 = KXL (1, LINE)
      K2 = KXL (2, LINE)
C
C  FIND THE NEXT LINE IN K1
C
      DO 100 I = 1, 4
         IL = LXK (I, K1)
         IF ((NXL (1, IL) .EQ. NEND) .AND.
     &      (NXL (2, IL) .NE. NODE)) THEN
            L1 = IL
            NNEW1 = NXL (2, IL)
            GOTO 110
         ELSEIF ((NXL (2, IL) .EQ. NEND) .AND.
     &      (NXL (1, IL) .NE. NODE)) THEN
            L1 = IL
            NNEW1 = NXL (1, IL)
            GOTO 110
         ENDIF
  100 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING NNEW1 **')
      RETURN
C
  110 CONTINUE
C
C  FIND THE NEXT LINE IN K2
C
      DO 120 I = 1, 4
         IL = LXK (I, K2)
         IF ((NXL (1, IL) .EQ. NEND) .AND.
     &      (NXL (2, IL) .NE. NODE)) THEN
            NNEW2 = NXL (2, IL)
            GOTO 130
         ELSEIF ((NXL (2, IL) .EQ. NEND) .AND.
     &      (NXL (1, IL) .NE. NODE)) THEN
            NNEW2 = NXL (1, IL)
            GOTO 130
         ENDIF
  120 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING NNEW2 **')
      RETURN
C
  130 CONTINUE
C
C  NOW DETERMINE WHICH OF THESE NODES IS N0 AND WHICH IS N2 BASED
C  ON THE FACT THAT THE CONNECTIVITY OF THE ELEMENTS LINES IS ALWAYS IN
C  COUNTER-CLOCKWISE ORDER
C
      DO 140 I = 1, 4
         IF (LXK (I, K1) .EQ. LINE) THEN
            I0 = I - 1
            I2 = I + 1
            IF (I .EQ. 1) THEN
               I0 = 4
            ELSEIF (I .EQ. 4) THEN
               I2 = 1
            ENDIF
            L0 = LXK (I0, K1)
            L2 = LXK (I2, K1)
            IF (L0 .EQ. L1) THEN
               N0 = NNEW1
               N2 = NNEW2
            ELSEIF (L2 .EQ. L1) THEN
               N0 = NNEW2
               N2 = NNEW1
            ELSE
               CALL MESAGE ('** PROBLEMS IN SETN02 FINDING N0 '//
     &            'AND N2 **')
            ENDIF
            GOTO 150
         ENDIF
  140 CONTINUE
      CALL MESAGE ('** PROBLEMS IN SETN02 FINDING LINE AGAIN **')
C
  150 CONTINUE
C
      RETURN
C
      END
