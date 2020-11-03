C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PCKSTR (NSTR, STR)
C=======================================================================

C   --*** PCKSTR *** (STRLIB) Remove all blanks from string
C   --   Written by Amy Gilkey - revised 03/21/88
C   --
C   --PCKSTR removes all embedded blanks (left-justified) from an array
C   --of strings.
C   --
C   --Parameters:
C   --   NSTR - IN - the number of strings to be packed
C   --   STR - IN/OUT - the array of strings, returned packed, may be up
C   --      to 1024 characters long

      INTEGER NSTR
      CHARACTER*(*) STR(*)

      PARAMETER (MXSTR = 1024)
      CHARACTER*(MXSTR) TMPSTR

      if (nstr .eq. 0) return

      LSTR = LENSTR (STR)
      if (lstr .gt. MXSTR) then
        call prterr ('PROGRAM', 'String is too long in SQZSTR')
        return
      end if

      DO 20 I = 1, NSTR
         LSTR = LENSTR (STR(I))

         IBLK = INDEX (STR(I), ' ')
   10    CONTINUE
         IF ((IBLK .GT. 0) .AND. (IBLK .LT. LSTR)) THEN
            TMPSTR = STR(I)(IBLK+1:LSTR)
            STR(I)(IBLK:LSTR) = TMPSTR
            LSTR = LSTR - 1
            IBLK = INDEX (STR(I), ' ')
            GOTO 10
         END IF

   20 CONTINUE

      RETURN
      END
