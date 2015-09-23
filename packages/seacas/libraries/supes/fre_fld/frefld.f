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
      SUBROUTINE FREFLD( KIN,KOUT,PROMPT,MFIELD,IOSTAT,NFIELD,KVALUE,
     *                   CVALUE,IVALUE,RVALUE )
      CHARACTER*(*) PROMPT,CVALUE(MFIELD)
      CHARACTER*132 LINE
      CHARACTER*132 PREFIX
      INTEGER KVALUE(MFIELD),IVALUE(MFIELD)
      REAL RVALUE(MFIELD)
C
************************************************************************
C
C     FREFLD INPUT SYSTEM - ANSI FORTRAN - USER INTERFACE ROUTINE
C
C     DESCRIPTION:
C     This routine is the main user interface to the SUPES Free Field
C     Input system. It obtains a record from the input stream, then
C     call FFISTR to parse the record into data fields.
C
C     FORMAL PARAMETERS:
C     KIN     INTEGER    Unit from which to read input.
C     KOUT    INTEGER    Unit to which to echo input.
C     PROMPT  CHARACTER  Prompt string.
C     MFIELD  INTEGER    Maximum number of data fields to be returned.
C     IOSTAT  INTEGER    ANSI FORTRAN I/O status.
C     NFIELD  INTEGER    Number of data fields found.
C     KVALUE  INTEGER    Translation states of the data fields:
C                           -1 = This is a null field.
C                            0 = This is a non-numeric field.
C                            1 = This is a REAL numeric field.
C                            2 = This is an INTEGER numeric field.
C     CVALUE  CHARACTER  Character values of the data fields.
C     RVALUE  REAL       Floating-point values of the data fields.
C     IVALUE  INTEGER    Integer values of the data fields.
C
C
C     ROUTINES CALLED:
C     GETINP            Get input line.
C     FFISTR            Parse input line.
C
************************************************************************
C
C PHASE 1: Initialize output arrays to their default values and zero
C          field counter. Set continuation flag to suppress futher
C          initialization by FFISTR.
C
      DO 300 I = 1 , MFIELD
         KVALUE(I) = -1
         CVALUE(I) = ' '
         RVALUE(I) = 0.
         IVALUE(I) = 0
  300 CONTINUE
      NFIELD = 0
      IDCONT = 1
C
C Initialize prompt to the caller's -
      PREFIX = PROMPT
      LPRE = LEN( PROMPT )
C
************************************************************************
C
C PHASE 2: Get the next input record via GETINP. Return to caller if an
C          end-of-file or error was detected by GETINP.
C          Re-enter here to process a continuation line.
C
  500 CONTINUE
C
C Get the next input line -
      CALL GETINP( KIN,KOUT,PREFIX(1:LPRE),LINE,IOSTAT )
C
C Return if I/O error or EOF detected -
      IF ( IOSTAT .NE. 0 ) RETURN
C
C Call FFISTR to parse input record -
      CALL FFISTR( LINE,MFIELD,IDCONT,NFIELD,KVALUE,CVALUE,IVALUE,
     *             RVALUE )
C
C If the continuation flag is set, define a continuation prompt and
C re-enter at phase 2. Otherwise, return to the caller -
      IF ( IDCONT .NE. 0 ) THEN
         IF ( PROMPT .EQ. 'AUTO' ) THEN
            PREFIX = '   *: '
            LPRE = 6
         ELSE
            IF ( LPRE .GT. 3 ) PREFIX(1:LPRE-3) = ' '
            IF ( LPRE .GE. 3 ) PREFIX(LPRE-2:LPRE-2) = '*'
         END IF
         GO TO 500
      END IF
C
      END
