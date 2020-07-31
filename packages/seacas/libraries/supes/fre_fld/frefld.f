C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE FREFLD( KIN,KOUT,PROMPT,MFIELD,IOSTAT,NFIELD,KVALUE,
     *                   CVALUE,IVALUE,RVALUE )
      CHARACTER*(*) PROMPT,CVALUE(MFIELD)
      CHARACTER*132 LINE
      CHARACTER*132 PREFIX
      INTEGER KVALUE(MFIELD),IVALUE(MFIELD)
      REAL RVALUE(MFIELD)

************************************************************************

C     FREFLD INPUT SYSTEM - ANSI FORTRAN - USER INTERFACE ROUTINE

C     DESCRIPTION:
C     This routine is the main user interface to the SUPES Free Field
C     Input system. It obtains a record from the input stream, then
C     call FFISTR to parse the record into data fields.

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

C     ROUTINES CALLED:
C     GETINP            Get input line.
C     FFISTR            Parse input line.

************************************************************************

C PHASE 1: Initialize output arrays to their default values and zero
C          field counter. Set continuation flag to suppress further
C          initialization by FFISTR.

      DO 300 I = 1 , MFIELD
         KVALUE(I) = -1
         CVALUE(I) = ' '
         RVALUE(I) = 0.
         IVALUE(I) = 0
  300 CONTINUE
      NFIELD = 0
      IDCONT = 1

C Initialize prompt to the caller's -
      PREFIX = PROMPT
      LPRE = LEN( PROMPT )

************************************************************************

C PHASE 2: Get the next input record via GETINP. Return to caller if an
C          end-of-file or error was detected by GETINP.
C          Re-enter here to process a continuation line.

  500 CONTINUE

C Get the next input line -
      CALL GETINP( KIN,KOUT,PREFIX(1:LPRE),LINE,IOSTAT )

C Return if I/O error or EOF detected -
      IF ( IOSTAT .NE. 0 ) RETURN

C Call FFISTR to parse input record -
      CALL FFISTR( LINE,MFIELD,IDCONT,NFIELD,KVALUE,CVALUE,IVALUE,
     *             RVALUE )

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

      END
