C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDMSH (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   IVIEW, JVIEW, NEWMOD,
     &   IDNPS, ISSNPS, IDESS, ISSESS, *)
C=======================================================================

C   --*** CMDMSH *** (MESH) Process display mode commands
C   --   Written by Amy Gilkey - revised 04/11/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   IVIEW - IN - the view number
C   --   JVIEW - IN - IVIEW (if non-zero) or a defined non-empty view number
C   --      (if any)
C   --   NEWMOD - OUT - the mode status of each view:
C   --     -1 = unchanged
C   --      0 = changed to default
C   --      n = changed to be like view n
C   --   IDNPS - IN - the node set ID for each set
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   IDESS - IN - the side set ID for each set
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Sets and uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET
C   --      of /MSHOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'mshopt.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      INTEGER NEWMOD(4)
      INTEGER IDNPS(*)
      INTEGER IDESS(*)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)

      INTEGER NDEFVW, IXVW
      LOGICAL MATSTR
      CHARACTER*(MXSTLN) WORD
      CHARACTER*(MXSTLN) MDEF, MNUM
      CHARACTER CDUM
      LOGICAL ISON
      INTEGER LTYP(-1:1)

      IF (MSHDEF(JVIEW) .EQ. 'EMPTY') THEN
         IF (.NOT. ((VERB .EQ. 'EMPTY') .OR.
     &      (VERB .EQ. 'DEFORM') .or. (verb .eq. 'UNDEFORM'))) THEN
            CALL PRTERR ('CMDERR', 'Specified view is empty')
            GOTO 160
         END IF
      END IF

      IF (VERB .EQ. 'EMPTY') THEN
         CALL FFADDC (VERB, INLINE)
         CALL SETMSH (IVIEW, 'EMPTY', CDUM, IDUM, LTYP,
     &      IDUM, IDUM, IDUM, IDUM, CDUM, CDUM, ISSNPS, ISSESS)

      ELSE IF ((VERB .EQ. 'DEFORM') .or. (verb .eq. 'UNDEFORM')) THEN
         if (verb .eq. 'UNDEFORM') then
            call prterr ('CMDREQ', 'Please use the DEFORM OFF command')
            verb = 'DEFORM'
            ison = .false.
         else
            CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *160)
         end if
         CALL FFADDC (VERB, INLINE)
         CALL FFADDO (ISON, INLINE)

         if ((mshdef(jview) .eq. 'DEFORM') .and. ison) then
            call prterr ('CMDWARN',
     &         'Use WIREFRAM to change the display mode')
         end if

         IF (ISON) THEN
            MDEF = 'DEFORM'
         ELSE
            MDEF = 'UNDEFORM'
         END IF

         IF (MSHDEF(JVIEW) .EQ. 'EMPTY') THEN
            CALL INIINT (3, 1, LTYP)
            CALL SETMSH (IVIEW, MDEF, 'NONE', MSHSEL, LTYP,
     &         0, IDUM, 0, IDUM, ' ', ' ', ISSNPS, ISSESS)
         END IF

         MSHDEF(JVIEW) = MDEF(:8)
         IF (IVIEW .EQ. 0) THEN
            DO 100 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               IF (JVIEW .NE. I) THEN
                  MSHDEF(I) = MSHDEF(JVIEW)
               END IF
  100       CONTINUE
         END IF

      ELSE IF (VERB .EQ. 'NUMBER') THEN
         CALL FFADDC (VERB, INLINE)
         CALL FFCHAR (IFLD, INTYP, CFIELD, 'ALL', WORD)
         IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE)
            MNUM = 'NONE'
         ELSE IF (MATSTR (WORD, 'NODES', 1)) THEN
            CALL FFADDC ('NODES', INLINE)
            MNUM = 'NODE'
         ELSE IF (MATSTR (WORD, 'ELEMENTS', 1)) THEN
            CALL FFADDC ('ELEMENTS', INLINE)
            MNUM = 'ELEMENT'
         ELSE IF (MATSTR (WORD, 'ALL', 1)) THEN
            CALL FFADDC ('ALL', INLINE)
            MNUM = 'ALL'
         ELSE
            MNUM = 'NONE'
            CALL PRTERR ('CMDERR',
     &         'Expected "NODES", "ELEMENTS", "ALL" or "OFF"')
         END IF

         MSHNUM(JVIEW) = MNUM(:8)
         IF (IVIEW .EQ. 0) THEN
            DO 110 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               IF (JVIEW .NE. I) THEN
                  MSHNUM(I) = MSHNUM(JVIEW)
               END IF
  110       CONTINUE
         END IF

      ELSE IF ((VERB .EQ. 'MLINES') .or. (verb .eq. 'OVERLAY')) THEN
         if (verb .eq. 'OVERLAY') then
            call prterr ('CMDREQ', 'Please use the MLINES command')
            verb = 'MLINES'
         end if
         CALL FFADDC (VERB, INLINE)

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
         IF (MATSTR (WORD, 'ON', 2)) THEN
            CALL FFADDC ('ON', INLINE)
            MLIN = MSHSEL
            ITYP = 1
         ELSE IF (MATSTR (WORD, 'DOTTED', 1)) THEN
            CALL FFADDC ('DOTTED', INLINE)
            MLIN = MSHSEL
            ITYP = 2
         ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE)
            MLIN = MSHDIV
            ITYP = 1
         ELSE
            CALL PRTERR ('CMDERR',
     &         'Expected "ON", "OFF" or "DOTTED"')
            GOTO 160
         END IF

         MSHLIN(JVIEW) = MLIN
         MLNTYP(1,JVIEW) = ISIGN (ITYP, MLNTYP(1,JVIEW))
         IF (IVIEW .EQ. 0) THEN
            DO 120 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               IF (JVIEW .NE. I) THEN
                  MSHLIN(I) = MSHLIN(JVIEW)
                  CALL CPYINT
     &               (3, MLNTYP(-1,JVIEW), MLNTYP(-1,I))
               END IF
  120       CONTINUE
         END IF

      ELSE IF (VERB .EQ. 'BOUNDARY') THEN
         CALL FFADDC (VERB, INLINE)
         IF (MLNTYP(1,JVIEW) .GT. 0) THEN
            CALL PRTERR ('CMDERR',
     &         'Command valid in painted mode only')
            GOTO 160
         END IF

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
         IF (MATSTR (WORD, 'ON', 2)) THEN
            CALL FFADDC ('ON', INLINE)
            IF (IABS (MLNTYP(1,JVIEW)) .EQ. 1) THEN
               MLIN = MSHSEL
            ELSE
               MLIN = MSHDIV
            END IF
            ITYP = 1
         ELSE IF (MATSTR (WORD, 'BLACK', 1)) THEN
            CALL FFADDC ('BLACK', INLINE)
            MLIN = MAX (MSHLIN(JVIEW), MSHDIV)
            ITYP = -1
         ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE)
            MLIN = MSHNON
            ITYP = 1
         ELSE
            CALL PRTERR ('CMDERR', 'Expected "ON", "OFF" or "BLACK"')
            GOTO 160
         END IF

         MSHLIN(JVIEW) = MLIN
         MLNTYP(-1,JVIEW) = ITYP
         MLNTYP( 0,JVIEW) = ITYP
         IF (IVIEW .EQ. 0) THEN
            DO 130 IVW = 1, NDEFVW (.FALSE.)
               I = IXVW (.FALSE., IVW)
               IF (JVIEW .NE. I) THEN
                  MSHLIN(I) = MSHLIN(JVIEW)
                  CALL CPYINT
     &               (3, MLNTYP(-1,JVIEW), MLNTYP(-1,I))
               END IF
  130       CONTINUE
         END IF

      ELSE IF ((VERB .EQ. 'NSETS') .OR. (VERB .EQ. 'SSETS')) THEN
         IF (VERB .EQ. 'NSETS') THEN
            CALL FFADDC (VERB, INLINE)
            CALL CKNONE (NUMNPS, .FALSE., 'node sets', *160)

            CALL RIXID (INLINE, IFLD, INTYP, CFIELD, IFIELD,
     &         'node set ID',
     &         NUMNPS, IDNPS, NNPSET(JVIEW), ISSNPS(1,JVIEW), *160)

            IF (IVIEW .EQ. 0) THEN
               DO 140 IVW = 1, NDEFVW (.FALSE.)
                  I = IXVW (.FALSE., IVW)
                  IF (JVIEW .NE. I) THEN
                     NNPSET(I) = NNPSET(JVIEW)
                     CALL CPYINT
     &                  (NNPSET(JVIEW), ISSNPS(1,JVIEW), ISSNPS(1,I))
                  END IF
  140          CONTINUE
            END IF

         ELSE IF (VERB .EQ. 'SSETS') THEN
            CALL FFADDC (VERB, INLINE)
            CALL CKNONE (NUMESS, .FALSE., 'side sets', *160)

            CALL RIXID (INLINE, IFLD, INTYP, CFIELD, IFIELD,
     &         'side set ID',
     &         NUMESS, IDESS, NESSET(JVIEW), ISSESS(1,JVIEW), *160)

            IF (IVIEW .EQ. 0) THEN
               DO 150 IVW = 1, NDEFVW (.FALSE.)
                  I = IXVW (.FALSE., IVW)
                  IF (JVIEW .NE. I) THEN
                     NESSET(I) = NESSET(JVIEW)
                     CALL CPYINT
     &                  (NESSET(JVIEW), ISSESS(1,JVIEW), ISSESS(1,I))
                  END IF
  150          CONTINUE
            END IF
         END IF

      ELSE
         INLINE = ' '
      END IF

      RETURN

  160 CONTINUE
      RETURN 1
      END
