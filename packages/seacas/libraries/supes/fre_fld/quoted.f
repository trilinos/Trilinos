C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
      SUBROUTINE QUOTED ( LINE, ILEFT, IRIGHT )
C
      CHARACTER*(*) LINE
C
      CALL STRIPB ( LINE, ILEFT, IRIGHT )
C
C     The first character is required to be a quote, so remove it.
C
      LINE(1:1) = ' '
      ILEFT = 2
      IBEGIN = 2
C
C     Begin loop looking for more quotes.  There should be at least 1 more.
C
 100  CONTINUE
         IQUOTE = INDEX ( LINE(IBEGIN:IRIGHT), '''' )
C
C        Has the quote ended within this record?
C
         IF ( IQUOTE .EQ. 0 ) THEN
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF
C
         IQUOTE = IQUOTE + IBEGIN - 1
         IF ( IQUOTE .EQ. IRIGHT ) THEN
C
C           The quote is at the end of the record.
C
            LINE(IRIGHT:IRIGHT) = ' '
            IRIGHT = IQUOTE - 1
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF
C
C        The quote is internal -- check for double quote.
C
         IF ( LINE(IQUOTE+1:IQUOTE+1) .NE. '''' ) THEN
C
C           The quote is single, thus ending the quoted string.
C
            LINE(IQUOTE:IQUOTE) = ' '
            IRIGHT = IQUOTE - 1
            IF ( ILEFT .GT. IRIGHT ) IRIGHT = 0
            RETURN
         END IF
C
C        The quote is a double quote.  Remove the repeated quote and loop.
C
         DO 10 I = IQUOTE, ILEFT, -1
            LINE(I+1:I+1) = LINE(I:I)
 10      CONTINUE
         LINE(ILEFT:ILEFT) = ' '
         ILEFT = ILEFT + 1
         IBEGIN = IQUOTE + 2
      GO TO 100
      END
