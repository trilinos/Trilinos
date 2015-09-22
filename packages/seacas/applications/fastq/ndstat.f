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

C $Id: ndstat.f,v 1.1 1990/11/30 11:12:34 gdsjaar Exp $
C $Log: ndstat.f,v $
C Revision 1.1  1990/11/30 11:12:34  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]NDSTAT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NDSTAT (NODE, LXN, ANGLE, JSTAT)
C***********************************************************************
C
C  SUBROUTINE NDSTAT = DETERMINES THE MOST APPROPRIATE STATUS OF A
C                      GIVEN NODE
C
C***********************************************************************
C
      DIMENSION LXN(4)
C
C  AN UNDECIDED NODE HAS BEEN FOUND - TEST ANGLE AND CONNECTIVITY
C
C  THE NODE IS ON THE BOUNDARY
      IF (LXN (2) .LT. 0) THEN
C
C IF THE NODE HAS LESS THAN FOUR LINES ATTACHED
C  CUTOFFS ARE:    0 TO 135 DEGREES = ROW END
C                135 TO 225 DEGREES = ROW SIDE
C                225 TO 290 DEGREES = ROW CORNER
C                OVER 290 DEGREES   = ROW REVERSAL
C
         IF (LXN (4) .LE. 0) THEN
            IF (ANGLE .LT. 2.3561945) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.0614548) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF
C
C IF THE NODE HAS FOUR LINES ATTACHED
C  CUTOFFS ARE:    0 TO 110 DEGREES = ROW END
C                110 TO 225 DEGREES = ROW SIDE
C                OVER 225 DEGREES   = ROW CORNER (NEARLY IMPOSSIBLE)
C
         ELSE
            IF (ANGLE .LT. 1.9198622) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE
               JSTAT = 5
            ENDIF
         ENDIF
C
C  THE NODE IS NOT ON THE BOUNDARY - CUTOFFS ARE ADJUSTED BASED
C  ON THE CONNECTIVITY AND THE ANGLE
C
      ELSE
C
C  ONLY TWO LINES ARE ATTACHED - LEAN TOWARDS A ROW CORNER NODE
C  OR A ROW END NODE
C
         IF (LXN(3) .EQ. 0) THEN
C
C  CUTOFFS ARE:    0 TO 135 DEGREES = ROW END
C                135 TO 210 DEGREES = ROW SIDE
C                210 TO 320 DEGREES = ROW CORNER
C                OVER 320 DEGREES   = ROW REVERSAL
C
            IF (ANGLE .LT. 2.3561945) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.6651914) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.5850536) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF
C
C  THREE LINES ARE ATTACHED - LEAN TOWARDS A ROW SIDE
C
         ELSEIF (LXN(4) .EQ. 0) THEN
C
C  CUTOFFS ARE:    0 TO 110 DEGREES = ROW END
C                110 TO 240 DEGREES = ROW SIDE
C                240 TO 320 DEGREES = ROW CORNER
C                OVER 320 DEGREES   = ROW REVERSAL (REALLY IMPOSSIBLE)
C
            IF (ANGLE .LT. 1.9198622) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 4.1887902) THEN
               JSTAT = 3
            ELSE IF (ANGLE .LT. 5.5850536) THEN
               JSTAT = 5
            ELSE
               JSTAT = 7
            ENDIF
C
C  FOUR LINES ARE ATTACHED - LEAN TOWARDS A ROW END NODE
C
         ELSE
C
C  CUTOFFS ARE:    0 TO 145 DEGREES = ROW END
C                145 TO 225 DEGREES = ROW SIDE
C                OVER 225 DEGREES   = ROW CORNER (REALLY IMPOSSIBLE)
C
            IF (ANGLE .LT. 2.5307274) THEN
               JSTAT = 1
            ELSE IF (ANGLE .LT. 3.9269908) THEN
               JSTAT = 3
            ELSE
               JSTAT = 5
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
C
      END
