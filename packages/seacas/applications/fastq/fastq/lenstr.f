C=======================================================================
      INTEGER FUNCTION LENSTR (STRING)
C=======================================================================
C$Id: lenstr.f,v 1.1 1999/01/27 15:17:48 gdsjaar Exp $
C$Log: lenstr.f,v $
CRevision 1.1  1999/01/27 15:17:48  gdsjaar
CAdded typical summary of mesh data on output.
C
CBetter filename handling
C
CCleaned up some character string handling
C
CRevision 1.1.1.1  1990/08/14 16:15:16  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:15  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:33  gdsjaar
c Initial revision
c 

C   --*** LENSTR *** (STRLIB) Return string length
C   --   Written by Amy Gilkey - revised 02/14/86
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
