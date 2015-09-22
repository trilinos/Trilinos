C Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

      PROGRAM TSTEXT
      CHARACTER*32 VERSN
      CHARACTER*8 DATE,TIME,NAME,LINE,HARD,SOFT
      character*1 b(5)
      DIMENSION A(5)
      CALL EXCPUS( CPU0 )

      CALL GSUPEV(VERSN)
      WRITE (*,'(A, A)') ' SUPES Version ', VERSN

      CALL EXREAD( 'TST: ',LINE,IOSTAT )
      CALL EXUPCS( LINE )
      PRINT *,'Input line = ',LINE

      CALL EXDATE( DATE )
      PRINT *,'Date = ',DATE

      CALL EXTIME( TIME )
      PRINT *,'Time = ',TIME

      CALL EXNAME( 1,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Unit 1 name = ',NAME(1:LN)

      CALL EXNAME( 10,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Unit 10 name = ',NAME(1:LN)

      CALL EXNAME( -1,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Symbol 1 = ',NAME(1:LN)

      CALL EXPARM( HARD,SOFT,MODE,KSCU,KNSU,IDAU )
      PRINT *,'Processor = ',HARD,'  System = ',SOFT,'  Mode = ',MODE
      PRINT *,'Character, Numeric, D/A Units: ',KSCU,KNSU,IDAU

      CALL EXMEMY( 10,LOCBLK,MEMRTN )
      PRINT *,'Memory block location and length: ',LOCBLK,MEMRTN

      MEMRTN = -998
      CALL EXMEMY( -10,LOCBLK,MEMRTN )
      PRINT *,'Memory block location and length: ',LOCBLK,MEMRTN

      NN = IXLNUM( A(5) ) - IXLNUM( A )
      PRINT *,'Numeric difference = ',NN

      NN = IXLCHR( B(5) ) - IXLCHR( B )
      PRINT *,'Character difference = ',NN

C ... Burn up some time so we can test the excpus function
      ra = 0.0
      do 100 i=1, 10000000
         ra = ra + sqrt(float(i)) 
 100  continue
      CALL EXCPUS( CPU1 )
      CPUS = CPU1 - CPU0
      PRINT *,'CPU time = ',CPUS
      END
