C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE STRIPB( STRING,ILEFT,IRIGHT )
      CHARACTER*(*) STRING

************************************************************************

C     FREFLD INPUT SYSTEM - ANSI FORTRAN - UTILITY ROUTINE

C     DESCRIPTION:
C     This routine strips leading and trailing blanks from a string. It
C     does not modify or copy the string, but simply returns the
C     location of the first and last non-blank characters.  If the
C     string is completely blank, ILEFT=LEN(STRING)+1 and IRIGHT=0 will
C     be returned.

C     FORMAL PARAMETERS:
C     STRING  CHARACTER  Any character string.
C     ILEFT   INTEGER    Position of first non-blank character.
C     IRIGHT  INTEGER    Position of last non-blank character.

************************************************************************

C ... Needed for 64-bit Solaris compile. Arg to CHAR must be
C     an integer of correct size.

      INTEGER TABC
      TABC = 9

C Get length of the string -
      LS = LEN( STRING )

C Find the first non-blank character -
      DO 10 N = 1 , LS
         IF ( STRING(N:N) .NE. ' '
     $        .AND. STRING(N:N) .NE. CHAR(TABC) ) GO TO 15
   10 CONTINUE
   15 CONTINUE
      ILEFT = N

C Find the last non-blank character -
      DO 20 N = LS , 1 , -1
         IF ( STRING(N:N) .NE. ' '
     $        .AND. STRING(N:N) .NE. CHAR(TABC) ) GO TO 25
   20 CONTINUE
   25 CONTINUE
      IRIGHT = N

      END
