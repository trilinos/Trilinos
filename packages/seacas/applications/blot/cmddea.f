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

C=======================================================================
      SUBROUTINE CMDDEA (VERB, INLINE, IFLD, INTYP, CFIELD, RFIELD,
     &                   NAMEEV, NALVAR, ALIVAL, *)
C=======================================================================

C   --*** CMDDEA *** (MESH) Process element birth/death commands
C   --   Written by Amy Gilkey - revised 03/11/88
C   --
C   --Parameters:
C   --   VERB   - I/O - the verbs for the SHOW command
C   --   INLINE - I/O - the parsed input line for the log file
C   --   IFLD   - I/O - the free-field reader index and fields
C   --   INTYP  - I/O - the free-field reader index and fields
C   --   CFIELD - I/O - the free-field reader index and fields
C   --   NAMEEV - IN  - the element variable names
C   --   NALVAR - I/O - the birth/death variable (0 if none)
C   --
C   --Common Variables:
C   --   Uses NUMEL, NVAREL, NSTEPW of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      REAL          RFIELD(*)
      CHARACTER*(*) NAMEEV(*)

      CHARACTER*(MXNAME) WORD
      LOGICAL ISON

      INTEGER NALOLD
      SAVE NALOLD
C      --NALOLD - the last birth/death variable (element variable index)

      DATA NALOLD / 0 /

      IF (VERB .EQ. 'DEATH') THEN
         CALL FFADDC (VERB, INLINE)
         IF (NSTEPW .LE. 0) THEN
            CALL PRTERR ('CMDERR',
     &         'No time steps with element variables are defined')
            GOTO 100
         END IF

         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *100)
         CALL FFADDO (ISON, INLINE)
         IF (ISON) THEN
            IF (NALOLD .GT. 0) THEN
               CALL FFCHAR (IFLD, INTYP, CFIELD, NAMEEV(NALOLD), WORD)
            ELSE
               CALL FFCHAR (IFLD, INTYP, CFIELD, 'DEATH', WORD)
            END IF
            CALL FFADDC (WORD, INLINE)
            IVAR = LOCSTR (WORD, NVAREL, NAMEEV)
            IF (IVAR .LE. 0) THEN
               CALL PRTERR ('CMDERR', 'Element variable "'
     &            // WORD(:LENSTR(WORD)) // '" does not exist')
               GOTO 100
            END IF

            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'alive value', 0.0, ALIVAL, *100)
            CALL FFADDR (ALIVAL, INLINE)

            CALL DBVIX_BL ('E', IVAR, NALVAR)
            NALOLD = IVAR
         ELSE
            NALVAR = 0
         END IF
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
