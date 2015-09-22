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
      PROGRAM FFRTEST
      PARAMETER (MFIELD=5)
      CHARACTER*16 CVALUE(MFIELD)
      INTEGER KVALUE(MFIELD),IVALUE(MFIELD)
      REAL RVALUE(MFIELD)
      CHARACTER*32 STRING

      CALL GSUPEV(STRING)
      WRITE (*,'(A, A)') 'SUPES VERSION ', STRING
   10 CALL FREFLD( 0,0,'AUTO',MFIELD,IOSTAT,NFIELD,KVALUE,
     *                   CVALUE,IVALUE,RVALUE )
      IF ( IOSTAT .EQ. 0 ) THEN
         PRINT 1000,NFIELD
         PRINT 2000,
     *         (I,KVALUE(I),CVALUE(I),RVALUE(I),IVALUE(I),I=1,MFIELD)
         GO TO 10
      END IF
 1000 FORMAT( ' NFIELD =',I5 /
     *        4X,'I',5X,'KV(I)',15X,'CV(I)',10X,'RV(I)',7X,'IV(I)' )
 2000 FORMAT( I5,I10,'  "',A,'"',1PG15.3,I12 )
      END
