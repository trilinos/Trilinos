C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
