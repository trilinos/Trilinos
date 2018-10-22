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

C $Log: cmdcon.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1996/06/21 16:07:02  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:55:54  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:12  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CMDCON (VERB, INLINE, IFLD, INTYP, IFIELD, RFIELD, *)
C=======================================================================

C   --*** CMDCON *** (DETOUR) Process contour commands
C   --   Written by Amy Gilkey - revised 04/04/88
C   --
C   --Parameters:
C   --   VERB   - I/O - the verb for the SHOW command
C   --   INLINE - I/O - the parsed input line for the log file
C   --   IFLD   - I/O - the free-field reader index
C   --   INTYP  - I/O - the free-field reader index
C   --   IFIELD - I/O - the free-field reader integer field
C   --   RFIELD - I/O - the free-field reader real field
C   --
C   --Common Variables:
C   --   Sets and uses CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC, CINTV of /CNTR/

      include 'cntr.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)

      LOGICAL FFEXST
      LOGICAL LDUM

      IF ((VERB .EQ. 'NCNTRS') .OR. (VERB .EQ. 'CRANGE')
     &   .OR. (VERB .EQ. 'CMIN') .OR. (VERB .EQ. 'CMAX')) THEN
         CALL FFADDC (VERB, INLINE)
         IF (VERB .EQ. 'NCNTRS') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD, 'number of contours',
     &                   6, NCNTR, *150)
            CALL FFADDI (NCNTR, INLINE)
            NCNTR = MAX (1, NCNTR)
            IF (NCNTR .GE. 256) THEN
               CALL PRTERR ('CMDWARN',
     &              'Number of contours reduced to 255')
               NCNTR = 255
            END IF

         ELSE IF ((VERB .EQ. 'CRANGE') .OR. (VERB .EQ. 'CMIN')) THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'contour minimum', CMIN, CMINX, *150)
            CALL FFADDR (CMINX, INLINE)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'contour maximum', CMAX, CMAX, *150)
            CALL FFADDR (CMAX, INLINE)
            CMIN = CMINX

         ELSE IF (VERB .EQ. 'CMAX') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'contour maximum', CMAX, CMAX, *150)
            CALL FFADDR (CMAX, INLINE)
         END IF

         IF (VERB .NE. 'NCNTR') CINTOK = .FALSE.

         CALL ADJCON (.FALSE.)

      ELSE IF (VERB .EQ. 'CSHIFT') THEN
         CALL FFADDC (VERB, INLINE)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'contour value', CMIN, CVAL, *150)
         CALL FFADDR (CVAL, INLINE)

         CINTOK = .FALSE.

         N = NINT ((CVAL - CMIN) / DELC)
         DIFF = CMIN + N * DELC - CVAL
         CMIN = CMIN - DIFF
         CMAX = CMAX - DIFF

      ELSE IF (VERB .EQ. 'DELCNTR') THEN
         CALL FFADDC (VERB, INLINE)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'contour interval', DELC, DELCX, *150)
         CALL FFADDR (DELCX, INLINE)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'contour minimum', CMIN, CMIN, *150)
         CALL FFADDR (CMIN, INLINE)
         DELC = DELCX

         CINTOK = .FALSE.

         CALL ADJCON (.TRUE.)

      ELSE IF (VERB .EQ. 'CINTV') THEN
         CALL FFADDC (VERB, INLINE)
         IF (NCNTR .GE. 256) THEN
            CALL PRTERR ('CMDWARN',
     &         'Number of contours reduced to 255')
            NCNTR = 255
         END IF

         IF (.NOT. CINTOK) THEN
            DO 100 I = 1, NCNTR+1
               CINTV(I) = CNTRI (I)
  100       CONTINUE
         END IF

         CINTOK = .TRUE.

         NC = 0
  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            NC = NC + 1
            IF (NC .GT. NCNTR+1) THEN
               CALL PRTERR ('CMDWARN',
     &            'Too many contour values given, ignored')
               GOTO 130
            END IF
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'contour value', CINTV(NC), CINTV(NC), *120)
            CALL FFADDR (CINTV(NC), INLINE)
  120       CONTINUE
            GOTO 110
         END IF
  130    CONTINUE

         DO 140 I = NCNTR+2, 256
            CINTV(I) = 0.0
  140    CONTINUE

         CALL ADJCON (.TRUE.)

         CALL CKCNTR (LDUM)
      END IF

      RETURN

  150 CONTINUE
      RETURN 1
      END
