C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDMVW (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   UNMESH, NEWMOD, ISSNPS, ISSESS, *)
C=======================================================================

C   --*** CMDMVW *** (MESH) Process multiple view commands
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   UNMESH - IN - the mesh limits
C   --   NEWMOD - OUT - the mode status of each view:
C   --     -1 = unchanged
C   --      0 = changed to default
C   --      n = changed to be like view n
C   --   ISSNPS - IN/OUT - the indices of the selected node sets
C   --   ISSESS - IN/OUT - the indices of the selected side sets
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Sets and uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET
C   --      of /MSHOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'mshopt.blk'
      include 'views.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      REAL UNMESH(KFAR)
      INTEGER NEWMOD(4)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      CHARACTER*(MXSTLN) WORD
      LOGICAL SYMCMD
      LOGICAL NEWSYM
      LOGICAL FFMATC, MATSTR

      IF ((VERB .EQ. 'XSYM') .OR. (VERB .EQ. 'XVIEW')) THEN
         CALL FFADDC (VERB, INLINE)
         SYMCMD = VERB(2:4) .EQ. 'SYM'
         IF (IS3DIM) THEN
            IF (SYMCMD) THEN
               CALL PRTERR ('CMDERR', 'Command allowed in 2D only')
               GOTO 100
            END IF
         END IF

         IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE)
            IF (MSHDEF(1) .EQ. 'NONE') THEN
               CALL PRTERR ('CMDERR',
     &            'Vertically divided views are not defined')
               GOTO 100
            END IF

C         --Move main views from left to right, if needed
            IF (XISSYM .AND. (.NOT. LFTSYM)) THEN
               CALL CPYMSH (2, 1, ISSNPS, ISSESS)
               IF (MSHDEF(3) .NE. 'NONE')
     &            CALL CPYMSH (4, 3, ISSNPS, ISSESS)
            END IF

C         --Delete left views
            CALL CPYMSH (1, 0, ISSNPS, ISSESS)
            IF (MSHDEF(3) .NE. 'NONE')
     &         CALL CPYMSH (3, 0, ISSNPS, ISSESS)

C         --Reset symmetric views if one view
            XISSYM = .FALSE.
            IF (MSHDEF(4) .EQ. 'NONE') MULTIM = .FALSE.

         ELSE
            IF (SYMCMD) THEN

C            --Get symmetry axis
               IF (LFTSYM) THEN
                  CALL FFCHAR (IFLD, INTYP, CFIELD, 'LEFT', WORD)
               ELSE
                  CALL FFCHAR (IFLD, INTYP, CFIELD, 'RIGHT', WORD)
               END IF
               IF (MATSTR (WORD, 'LEFT', 1)) THEN
                  CALL FFADDC ('LEFT', INLINE)
                  NEWSYM = .TRUE.
               ELSE IF (MATSTR (WORD, 'RIGHT', 1)) THEN
                  CALL FFADDC ('RIGHT', INLINE)
                  NEWSYM = .FALSE.
               ELSE
                  CALL PRTERR ('CMDERR',
     &               'Expected "LEFT", "RIGHT" or "OFF"')
                  GOTO 100
               END IF

               IF (FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)) THEN
                  CALL PICK2D ('point on symmetry axis', .TRUE.,
     &               .TRUE., IFLD, INTYP, RFIELD,
     &               XAXSYM, RDUM, *100)
               ELSE
                  IF (NEWSYM .EQV. LFTSYM) THEN
                     XS = XAXSYM
                  ELSE IF (NEWSYM) THEN
                     XS = UNMESH(KLFT)
                  ELSE IF (.NOT. NEWSYM) THEN
                     XS = UNMESH(KRGT)
                  END IF
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &               'symmetry axis', XS, XAXSYM, *100)
               END IF
               CALL FFADDR (XAXSYM, INLINE)

               XISSYM = .TRUE.
               LFTSYM = NEWSYM

            ELSE
               CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
               CALL FFADDC (WORD, INLINE)
               IF (WORD .NE. 'ON') THEN
                  CALL PRTERR ('CMDERR', 'Expected "ON" or "OFF"')
                  GOTO 100
               END IF

               XISSYM = .FALSE.
            END IF

C         --Copy right views to left views, if not defined
            IF (MSHDEF(1) .EQ. 'NONE') THEN
               CALL CPYMSH (1, 2, ISSNPS, ISSESS)
               IF (MSHDEF(4) .NE. 'NONE')
     &            CALL CPYMSH (3, 4, ISSNPS, ISSESS)
            END IF
         END IF

      ELSE IF ((VERB .EQ. 'YSYM') .OR. (VERB .EQ. 'YVIEW')) THEN
         CALL FFADDC (VERB, INLINE)
         SYMCMD = VERB(2:4) .EQ. 'SYM'
         IF (IS3DIM) THEN
            IF (SYMCMD) THEN
               CALL PRTERR ('CMDERR', 'Command allowed in 2D only')
               GOTO 100
            END IF
         END IF

         IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE)
            IF (MSHDEF(4) .EQ. 'NONE') THEN
               CALL PRTERR ('CMDERR',
     &            'Horizontally divided views are not defined')
               GOTO 100
            END IF

C         --Move main views from bottom to top, if needed
            IF (YISSYM .AND. (.NOT. BOTSYM)) THEN
               CALL CPYMSH (2, 4, ISSNPS, ISSESS)
               IF (MSHDEF(3) .NE. 'NONE')
     &            CALL CPYMSH (1, 3, ISSNPS, ISSESS)
            END IF

C         --Delete bottom views
            CALL CPYMSH (4, 0, ISSNPS, ISSESS)
            IF (MSHDEF(3) .NE. 'NONE')
     &         CALL CPYMSH (3, 0, ISSNPS, ISSESS)

C         --Reset symmetric views if one view
            YISSYM = .FALSE.
            IF (MSHDEF(1) .EQ. 'NONE') MULTIM = .FALSE.

         ELSE
            IF (SYMCMD) THEN

C            --Get symmetry axis
               IF (BOTSYM) THEN
                  CALL FFCHAR (IFLD, INTYP, CFIELD, 'BOTTOM', WORD)
               ELSE
                  CALL FFCHAR (IFLD, INTYP, CFIELD, 'TOP', WORD)
               END IF
               IF (MATSTR (WORD, 'BOTTOM', 1)) THEN
                  CALL FFADDC ('BOTTOM', INLINE)
                  NEWSYM = .TRUE.
               ELSE IF (MATSTR (WORD, 'TOP', 1)) THEN
                  CALL FFADDC ('TOP', INLINE)
                  NEWSYM = .FALSE.
               ELSE
                  CALL PRTERR ('CMDERR',
     &               'Expected "BOTTOM", "TOP" or "OFF"')
                  GOTO 100
               END IF

               IF (FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)) THEN
                  CALL PICK2D ('point on symmetry axis', .TRUE.,
     &               .TRUE., IFLD, INTYP, RFIELD,
     &               RDUM, YAXSYM, *100)
               ELSE
                  IF (NEWSYM .EQV. BOTSYM) THEN
                     YS = YAXSYM
                  ELSE IF (NEWSYM) THEN
                     YS = UNMESH(KBOT)
                  ELSE IF (.NOT. NEWSYM) THEN
                     YS = UNMESH(KTOP)
                  END IF
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &               'symmetry axis', YS, YAXSYM, *100)
               END IF
               CALL FFADDR (YAXSYM, INLINE)

               YISSYM = .TRUE.
               BOTSYM = NEWSYM

            ELSE
               CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
               CALL FFADDC (WORD, INLINE)
               IF (WORD .NE. 'ON') THEN
                  CALL PRTERR ('CMDERR', 'Expected "ON" or "OFF"')
                  GOTO 100
               END IF

               YISSYM = .FALSE.
            END IF

C         --Copy top views to bottom views, if not defined
            IF (MSHDEF(4) .EQ. 'NONE') THEN
               CALL CPYMSH (4, 2, ISSNPS, ISSESS)
               IF (MSHDEF(1) .NE. 'NONE')
     &            CALL CPYMSH (3, 1, ISSNPS, ISSESS)
            END IF
         END IF
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
