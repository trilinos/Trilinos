C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: grikey.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:28  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:41  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRIKEY (PROMPT, DX, DY, KEY, *)
C=======================================================================

C   --*** GRIKEY *** (GRPLIB) Position cursor and wait for a key (PLT)
C   --   Written by Amy Gilkey - revised 04/28/88
C   --
C   --GRIKEY positions a cursor and waits until a key is pressed.
C   --A message may be output before the cursor is positioned.  The
C   --cursor position and the key are returned.
C   --
C   --The first time this routine is called, a warning message is output.
C   --
C   --Parameters:
C   --   PROMPT - IN - prompt for key if non-blank
C   --   DX, DY - IN/OUT - the device coordinates of cursor
C   --   KEY - OUT - returned key pressed; upper-case
C   --   * - the alternate return if cancel requested or cursor input error

C   --Routines Called:
C   --   PLTCRS - (PLTLIB) Select cursor position
C   --   PLTFLU - (PLTLIB) Flush buffer
C   --   LENSTR - (STRLIB) Find string length

      CHARACTER*(*) PROMPT
      REAL DX, DY
      CHARACTER KEY

      LOGICAL FIRST
      SAVE FIRST

      DATA FIRST / .TRUE. /

      CALL PLTFLU

      IF (FIRST) THEN
         WRITE (*, 10000)
         FIRST = .FALSE.
      END IF

      IF (PROMPT .NE. ' ') THEN
         WRITE (*, '(1X, A)') PROMPT(:LENSTR(PROMPT))
      END IF

      CALL PLTCRS (DX, DY, KEY)

      IF ((DX .LT. 0.0) .OR. (DX .GT. 1.0)
     &   .OR. (DY .LT. 0.0) .OR. (DY .GT. 1.0)) THEN
         WRITE (*, 10010)
         RETURN 1
      END IF

      CALL EXUPCS (KEY)

      RETURN
10000  FORMAT (/,
     &   15X,'*** Press SPACE or a letter to pick a point ***', /)
10010  FORMAT (/,
     &   ' *** ERROR on cursor input ***', /)
      END
