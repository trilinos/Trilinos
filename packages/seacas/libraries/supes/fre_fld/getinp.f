C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE GETINP( KIN,KOUT,PROMPT,LINE,IOSTAT )
      CHARACTER*(*) PROMPT,LINE
      CHARACTER*132 PREFIX
      DATA KOUNT /0/

************************************************************************

C     FREFLD INPUT SYSTEM - ANSI FORTRAN - UTILITY ROUTINE

C     DESCRIPTION:
C     This routine is performs all I/O for the SUPES Free Field Input
C     system. Its operation depends on the input and output units
C     specified by the caller as follows.

C     KIN     KOUT    Source           Echo
C     ------------------------------------------------------------------
C     0       0       Standard Input   Standard Output
C     0       M       Standard Input   Standard Output and File (M)
C     N       M       File (N)         File (N)
C     N       0       File (N)         none

C     If the prompt string is 'AUTO' this routine will generate a prompt
C     of the form '   n: ', where  "n" is the current input line number.
C     Only lines read under the AUTO feature are counted.

C     This routine does not restrict the length of the input string, but
C     no more than 132 characters, including the prompt, will be echoed.

C     FORMAL PARAMETERS:
C     KIN     INTEGER    Unit from which to read input.
C     KOUT    INTEGER    Unit to which to echo input.
C     PROMPT  CHARACTER  Prompt string.
C     LINE    CHARACTER  Input record.
C     IOSTAT  INTEGER    ANSI FORTRAN I/O status.

C     ROUTINES CALLED:
C     STRIPB            Strip leading/trailing blanks from a string.
C     EXREAD            Prompt, read, and echo and input record.

************************************************************************

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

C Read the input line -
      IF ( KIN .EQ. 0 ) THEN
         CALL EXREAD( PREFIX(1:LPRE),LINE,IOSTAT )
      ELSE
         READ( KIN,2000,IOSTAT=IOSTAT ) LINE
      END IF

C Return if I/O error or EOF detected -
      IF ( IOSTAT .NE. 0 ) RETURN

C Find the last non-blank character -
      CALL STRIPB( LINE,ILEFT,IRIGHT )

C Truncate the string for echo, if necessary -
      IRIGHT = MIN( IRIGHT,132-LPRE )

C Echo the input line, if requested -
      IF ( KOUT .GT. 0 ) THEN
         IF ( IRIGHT .EQ. 0 ) THEN
            WRITE( KOUT,3000 ) PREFIX(1:LPRE)
         ELSE
            WRITE( KOUT,3000 ) PREFIX(1:LPRE),LINE(1:IRIGHT)
         END IF
      END IF

 1000 FORMAT( I4,': ' )
 2000 FORMAT( A )
 3000 FORMAT( 1X,A,A )

      END
