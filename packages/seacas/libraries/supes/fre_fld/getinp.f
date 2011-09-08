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
      SUBROUTINE GETINP( KIN,KOUT,PROMPT,LINE,IOSTAT )
      CHARACTER*(*) PROMPT,LINE
      CHARACTER*132 PREFIX
      DATA KOUNT /0/
C
************************************************************************
C
C     FREFLD INPUT SYSTEM - ANSI FORTRAN - UTILITY ROUTINE
C
C     DESCRIPTION:
C     This routine is performs all I/O for the SUPES Free Field Input
C     system. Its operation depends on the input and output units
C     specified by the caller as follows.
C
C     KIN     KOUT    Source           Echo
C     ------------------------------------------------------------------
C     0       0       Standard Input   Standard Ouput
C     0       M       Standard Input   Standard Ouput and File (M)
C     N       M       File (N)         File (N)
C     N       0       File (N)         none
C
C     If the prompt string is 'AUTO' this routine will generate a prompt
C     of the form '   n: ', where  "n" is the current input line number.
C     Only lines read under the AUTO feature are counted.
C
C     This routine does not restrict the length of the input string, but
C     no more than 132 characters, including the prompt, will be echoed.
C
C     FORMAL PARAMETERS:
C     KIN     INTEGER    Unit from which to read input.
C     KOUT    INTEGER    Unit to which to echo input.
C     PROMPT  CHARACTER  Prompt string.
C     LINE    CHARACTER  Input record.
C     IOSTAT  INTEGER    ANSI FORTRAN I/O status.
C
C     ROUTINES CALLED:
C     STRIPB            Strip leading/trailing blanks from a string.
C     EXREAD            Prompt, read, and echo and input record.
C
************************************************************************
C
C Generate a prompt if autoincrement mode is requested; in any case
C PREFIX(1:LPRE) will contain the prompt -
      IF( PROMPT .EQ. 'AUTO' ) THEN
         KOUNT = KOUNT + 1
C Wrap-around counter, if it exceeds 4 digits -
         IF ( KOUNT .EQ. 10000 ) KOUNT = 0
         WRITE( PREFIX,1000 ) KOUNT
         LPRE = 6
      ELSE
         PREFIX = PROMPT
         LPRE = LEN( PROMPT )
      END IF
C
C Read the input line -
      IF ( KIN .EQ. 0 ) THEN
         CALL EXREAD( PREFIX(1:LPRE),LINE,IOSTAT )
      ELSE
         READ( KIN,2000,IOSTAT=IOSTAT ) LINE
      END IF
C
C Return if I/O error or EOF detected -
      IF ( IOSTAT .NE. 0 ) RETURN
C
C Find the last non-blank character -
      CALL STRIPB( LINE,ILEFT,IRIGHT )
C
C Truncate the string for echo, if necessary -
      IRIGHT = MIN( IRIGHT,132-LPRE )
C
C Echo the input line, if requested -
      IF ( KOUT .GT. 0 ) THEN
         IF ( IRIGHT .EQ. 0 ) THEN
            WRITE( KOUT,3000 ) PREFIX(1:LPRE)
         ELSE
            WRITE( KOUT,3000 ) PREFIX(1:LPRE),LINE(1:IRIGHT)
         END IF
      END IF
C
 1000 FORMAT( I4,': ' )
 2000 FORMAT( A )
 3000 FORMAT( 1X,A,A )
C
      END
