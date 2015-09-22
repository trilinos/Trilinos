C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C $Log: setfor.f,v $
C Revision 1.4  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  1998/06/12 15:53:38  gdsjaar
C 1. Problem with TIMES array. Blot accesses a dummy timestep even if
C there were no timesteps on the database. Array wasn't allocated, so
C writing off into never-never land.
C
C 2. Inconsistency among some common blocks. Some places weren't using
C the include but had the definition hardwired in. Removed those.
C
C 3. Added 'EXTERNAL BLKDAT' to all routines that used data values set
C in BLKDAT
C
C 4. Cleanup of some A vs. IA argument passing.
C
C Revision 1.2  1994/04/26 16:37:58  gdsjaar
C Fixed problem with setting background and foreground colors
C
c Revision 1.1  1994/04/07  20:11:33  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.2  1990/12/14  08:57:10  gdsjaar
c Added RCS Id and Log to all files
c
C============================================================================
      SUBROUTINE SETFOR (IFUNC, INLINE, IFLD, INTYP, CFIELD, *)
C============================================================================

C   --*** SETFOR ***  (BLOT) Process FOREGROUND command
C   --   Written by John Glick - 2/27/89
C   --
C   --Parameters:
C   --   IFUNC  - IN - = 1 if the call is to parse the FOREGROUND
C   --                     command.
C   --                 = 2 if the call is to reset the foreground color
C   --                     to the default color.
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   IFLD, INTYP, CFIELD, - IN/OUT - the free-field reader
C   --          index and charcter field.
C   --   * - return statement if command error; message is printed.

      include 'params.blk'
      INTEGER IFUNC
      CHARACTER*(*) INLINE(*)
      INTEGER IFLD, INTYP(*)
      CHARACTER*(*) CFIELD(*)

      include 'plcolr.blk'
      include 'plcol2.blk'

      LOGICAL FFEXST
      INTEGER IDCOL, LOCSTR

      CHARACTER*(mxstln) COLA, COLF
      include 'cmap-lst.blk'
      INTEGER DEFFOR
      SAVE DEFFOR

      DATA DEFFOR / 2 /

c      DATA FORGND /'WHITE   '/

C *****************************************************************

      IF (IFUNC .EQ. 1) THEN

         IF (FFEXST (IFLD, INTYP)) THEN

C              Check that next field has characters in it.

            IF (INTYP(IFLD) .GE. 0) THEN
               COLA = CFIELD(IFLD)
               IFLD = IFLD + 1
               CALL ABRSTR (COLF, COLA, COLLST)
               IF (COLF .EQ. ' ') THEN
                  WRITE (*, 10000) COLA
10000              FORMAT (1X, A, ' not a valid color name.')
                  GO TO 100
               ELSE
                  IDCOL = LOCSTR (COLF, NCOLOR, COLLST)
                  CALL FFADDC (COLF, INLINE)
                  IF (IDCOL .GT. 0) THEN
                     FORGND = COLF
                     IDFOR = IDCOL
                     IDFORT = IDCOL
                  ELSE
                     WRITE (*, 10000) COLA
                     GO TO 100
                  ENDIF
               ENDIF
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected color name following FOREGROUND command')
               GOTO 100
            ENDIF
         ELSE

            CALL PRTERR ('CMDERR',
     &         'Expected color name following FOREGROUND command')
            GOTO 100

         ENDIF

      ELSE IF (IFUNC .EQ. 2) THEN

         IDFOR = DEFFOR
         IDFORT = DEFFOR
         FORGND = 'WHITE   '

      ENDIF

      RETURN

  100 CONTINUE
      RETURN 1

      END
