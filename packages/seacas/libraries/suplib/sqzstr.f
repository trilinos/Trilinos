C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SQZSTR (STRING, LSTR)
C=======================================================================

C   --*** SQZSTR *** (STRLIB) Remove extra blanks from string
C   --   Written by Amy Gilkey - revised 06/02/87
C   --
C   --SQZSTR deletes leading and extra blanks within a string.
C   --To prevent problems, an empty string is returned with a length of 1.
C   --
C   --Parameters:
C   --   STRING - IN/OUT - the string to be compressed, returned, may be
C   --      up to 1024 characters long
C   --   LSTR - OUT - the new string length

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      CHARACTER*(*) STRING
      INTEGER LSTR

      PARAMETER (MXSTR = 1024)
      CHARACTER*(MXSTR) TMPSTR

      LSTR = 1
      IF (STRING .EQ. ' ') RETURN

      LSTR = LENSTR (STRING)
      if (lstr .gt. MXSTR) then
        call prterr ('PROGRAM', 'String is too long in SQZSTR')
        return
      end if

      IBLK = INDEX (STRING, '  ')
  100 CONTINUE
      IF ((IBLK .GT. 0) .AND. (IBLK .LE. LSTR)) THEN

         INOBLK = IBLK + 2
  110    CONTINUE
         IF (STRING(INOBLK:INOBLK) .EQ. ' ') THEN
            INOBLK = INOBLK + 1
            GOTO 110
         END IF

         TMPSTR = STRING(INOBLK:LSTR)
         STRING(IBLK+1:LSTR) = TMPSTR
         LSTR = LSTR - (INOBLK-IBLK-1)
         IBLK = INDEX (STRING, '  ')
         GOTO 100
      END IF

      IF (STRING(1:1) .EQ. ' ') THEN
         TMPSTR = STRING(2:LSTR)
         STRING(1:LSTR) = TMPSTR
         LSTR = LSTR - 1
      END IF

      RETURN
      END
