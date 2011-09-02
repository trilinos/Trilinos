C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      SUBROUTINE PACKT (TITLE,LENGTH)
C
C ... REMOVE MULTIPLE BLANKS FROM A TITLE OR LABEL
C
      CHARACTER*(*) TITLE
      CHARACTER*1 BLANK
      DATA BLANK/' '/
      I=1
      L=1
C
C ... SKIP LEADING BLANKS
C
   10 CONTINUE
      IF (TITLE(I:I) .NE. BLANK) GO TO 20
      I=I+1
      GO TO 10
   20 CONTINUE
      TITLE(L:L)=TITLE(I:I)
   30 CONTINUE
      IF (I .GE. LENGTH) GO TO 60
      I=I+1
      L=L+1
      TITLE(L:L)=TITLE(I:I)
      IF (TITLE(I:I) .EQ. BLANK) THEN
   40    CONTINUE
         IF (TITLE(I:I) .NE. BLANK .OR. I .GE. LENGTH) GO TO 50
         I=I+1
         GO TO 40
   50    CONTINUE
         L=L+1
         TITLE(L:L)=TITLE(I:I)
       END IF
      GO TO 30
   60 CONTINUE
      DO 70 I=L+1,LENGTH
         TITLE(I:I)=BLANK
   70    CONTINUE
      END
