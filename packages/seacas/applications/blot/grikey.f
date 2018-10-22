C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
