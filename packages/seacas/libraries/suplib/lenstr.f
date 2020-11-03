C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      INTEGER FUNCTION LENSTR (STRING)
C=======================================================================

C   --*** LENSTR *** (STRLIB) Return string length
C   --
C   --LENSTR returns the length of a passed string.  There can be blanks
C   --embedded within the string.  To prevent problems, an empty string
C   --is returned with a length of 1.
C   --
C   --Parameters:
C   --   STRING - IN - the string to find the length of

      CHARACTER*(*) STRING

      LENSTR = LEN(STRING)
      IF (STRING(LENSTR:LENSTR) .EQ. ' ') THEN
         N = INDEX (STRING, ' ')
   10    CONTINUE
         IF (STRING(N:) .NE. ' ') THEN
            N = N + INDEX (STRING(N+1:), ' ')
            GOTO 10
         END IF
         LENSTR = MAX (1, N-1)
      END IF

      RETURN
      END
