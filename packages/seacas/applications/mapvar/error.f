C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 

C=======================================================================
*DECK,ERROR
      SUBROUTINE ERROR (SUBNAM,MESSAG,LABEL1,I,LABEL2,J,LABEL3,WORD,    
     1ISTOP)
C
C     ******************************************************************
C
C     SUBROUTINE TO PRINT ERROR MESSAGE AND TERMINATE EXECUTION
C
C     Calls subroutine CLSFIL
C
C     Called by everything
C
C     ******************************************************************
C
      CHARACTER*(*) SUBNAM,MESSAG,LABEL1,LABEL2,LABEL3,WORD
C
      include 'tapes.blk'
C
C     ******************************************************************
C
      WRITE (NOUT, 60)
      WRITE (NTPOUT, 60)
      WRITE (NOUT, 10) SUBNAM
      WRITE (NTPOUT, 10) SUBNAM
      WRITE (NOUT, 20) MESSAG
      WRITE (NTPOUT, 20) MESSAG
      WRITE (NOUT, 30)
      WRITE (NTPOUT, 30)
      IF (LABEL1.NE.' ') THEN
        WRITE (NOUT, 40) LABEL1,I
        WRITE (NTPOUT, 40) LABEL1,I
      END IF
      IF (LABEL2.NE.' ') THEN
        WRITE (NOUT, 40) LABEL2,J
        WRITE (NTPOUT, 40) LABEL2,J
      END IF
      IF (LABEL3.NE.' ') THEN
        WRITE (NOUT, 50) LABEL3,WORD
        WRITE (NTPOUT, 50) LABEL3,WORD
      END IF
      WRITE (NOUT, 60)
      WRITE (NTPOUT, 60)
C
      IF (ISTOP.EQ.0) RETURN
C
      CALL CLSFIL
C
      STOP 'ERROR'
C
   10 FORMAT (/,10X,' ERROR FOUND IN - ' ,A)
   20 FORMAT (/,10X,' DESCRIPTION - ' ,A)
   30 FORMAT (/,10X,' RELEVANT PARAMETERS - ')
   40 FORMAT (/,15X, A, ' = ' ,I10)
   50 FORMAT (/,15X, A, ' = ', A)
   60 FORMAT (/,10X,'* * * * * * * * * * * * * * * * * * * * * * * * *',
     *' * * * * * ',/,10X,'* * * * * * * * * * * * * * * * * * * * * *',
     *' * * * * * * * *',/)
      END
