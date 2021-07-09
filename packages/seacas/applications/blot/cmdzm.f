C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDZM (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NEWZM, SETTIC, MAPND, A, *)
C=======================================================================

C   --*** CMDZM *** (MESH) Process scaling commands
C   --   Written by Amy Gilkey - revised 04/27/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NEWZM - IN/OUT - true iff a new zoom window or scaling is set
C   --   SETTIC - IN/OUT - true iff the tick interval is set by the user
C   --   A - the dynamic memory array
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses DFAC of /DEFORM/
C   --   Sets ZMMESH, RDMESH, TICMSH, MSCTYP, SQMESH of /MSHLIM/
C   --   Uses UNMESH, RDMESH of /MSHLIM/
C   --   Uses ROTMAT, ROTCEN, EYE of /ROTOPT/
C   --   Sets and uses NZMON,XZM,YZM,ZZM,RADZM,NDZMID OF /NODZOM/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'nodzom.blk'

      DIMENSION A(*)

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER       INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER       IFIELD(*)
      REAL          RFIELD(*)
      LOGICAL NEWZM
      LOGICAL SETTIC
      INTEGER       MAPND(*)

      CHARACTER*(MXSTLN) WORD
      REAL RNUM(KTOP)
      LOGICAL FFEXST, FFNUMB, FFMATC, MATSTR

      IF (VERB .EQ. 'ZOOM') THEN
         CALL FFADDC (VERB, INLINE)

         IF (FFEXST (IFLD, INTYP)) THEN
            IF (FFNUMB (IFLD, INTYP)) THEN
               WORD = 'limits'
            ELSE
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            END IF
         ELSE
            IF (.NOT. IS3DIM) THEN
               WORD = 'MESH'
            ELSE
               WORD = 'EACH'
            END IF
         END IF
         IF (MATSTR (WORD, 'RESET', 3)) THEN
            IF (.NOT. IS3DIM) THEN
               WORD = 'MESH'
            ELSE
               WORD = 'EACH'
            END IF
         END IF

         IF (MATSTR (WORD, 'EACH', 1)) THEN
            CALL FFADDC ('EACH', INLINE)
            MSCTYP = 'EACH'
            NZMON = .FALSE.
         ELSE IF (MATSTR (WORD, 'MESH', 1)) THEN
            IF (IS3DIM) THEN
               CALL PRTERR ('CMDERR', 'Command allowed in 2D only')
               GOTO 100
            END IF
            CALL FFADDC ('MESH', INLINE)
            MSCTYP = 'MESH'
            NZMON = .FALSE.

         ELSE IF (MATSTR (WORD, 'ROTATION', 3)) THEN
            IF (.NOT. IS3DIM) THEN
               CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
               GOTO 100
            END IF
            CALL FFADDC ('ROTATION', INLINE)
            MSCTYP = 'ROTATION'
            NZMON = .FALSE.

         ELSE IF (MATSTR (WORD, 'TRANSLAT', 1)) THEN
            CALL FFADDC ('TRANSLAT', INLINE)

            if (ffmatc (ifld, intyp, cfield, 'KEY', 1)) then
               call prterr ('CMDREQ', 'Please use CURSOR not KEY')
               ifld = ifld - 1
               cfield(ifld) = 'CURSOR'
            end if
            IF (FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)) THEN
               CALL PICK2D ('center of window', .TRUE.,
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            XCEN, YCEN, *100)
            ELSE
               IF ((MSCTYP .NE. 'ZOOM') .AND. (MSCTYP .NE. 'MESH')) THEN
                  CALL FFNEED (IFLD, INTYP, 'R', 2,
     &               'center window coordinates', *100)
               END IF
               XSVCEN = RDMESH(KLFT)
     &            + 0.5 * (RDMESH(KRGT) - RDMESH(KLFT))
               YSVCEN = RDMESH(KBOT)
     &            + 0.5 * (RDMESH(KTOP) - RDMESH(KBOT))
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'horizontal center', XSVCEN, XCEN, *100)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'vertical center', YSVCEN, YCEN, *100)
            END IF
            CALL FFADDR (XCEN, INLINE)
            CALL FFADDR (YCEN, INLINE)

            MSCTYP = 'ZOOM'
            X = XSVCEN - XCEN
            RDMESH(KLFT) = RDMESH(KLFT) - X
            RDMESH(KRGT) = RDMESH(KRGT) - X
            Y = YSVCEN - YCEN
            RDMESH(KBOT) = RDMESH(KBOT) - Y
            RDMESH(KTOP) = RDMESH(KTOP) - Y
            NZMON = .FALSE.

         ELSE IF (MATSTR (WORD, 'IN', 1)) THEN
            CALL FFADDC ('IN', INLINE)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'magnification', 1.0, R, *100)
            CALL FFADDR (R, INLINE)

            XCEN = RDMESH(KLFT)
     &         + 0.5 * (RDMESH(KRGT) - RDMESH(KLFT))
            YCEN = RDMESH(KBOT)
     &         + 0.5 * (RDMESH(KTOP) - RDMESH(KBOT))
            XTOT = (RDMESH(KRGT) - RDMESH(KLFT)) / R
            YTOT = (RDMESH(KTOP) - RDMESH(KBOT)) / R

            MSCTYP = 'ZOOM'
            RDMESH(KLFT) = XCEN - 0.5 * XTOT
            RDMESH(KRGT) = XCEN + 0.5 * XTOT
            RDMESH(KBOT) = YCEN - 0.5 * YTOT
            RDMESH(KTOP) = YCEN + 0.5 * YTOT
            NZMON = .FALSE.

         ELSE IF (MATSTR (WORD, 'limits', 1)
     &      .OR. MATSTR (WORD, 'CURSOR', 1)
     &      .or. matstr (word, 'KEY', 1)) THEN
            if (matstr (word, 'KEY', 1)) then
               call prterr ('CMDREQ', 'Please use CURSOR not KEY')
               word = 'CURSOR'
            end if
            IF (MATSTR (WORD, 'CURSOR', 1)) THEN
               CALL PICK2D ('bottom left corner', .TRUE.,
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            RNUM(KLFT), RNUM(KBOT), *100)
               CALL PICK2D ('top right corner', .TRUE.,
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            RNUM(KRGT), RNUM(KTOP), *100)
            ELSE

               IF (IS3DIM) THEN
                  CALL FFNEED (IFLD, INTYP, 'R', 4,
     &               'window coordinates of opposite corners', *100)
               END IF

C            --Calculate the default mesh limits
               IF (.NOT. IS3DIM) THEN
                  IF (DFAC .EQ. 0.0) THEN
                     CALL EXPLIM (NDIM, UNMESH, RNUM)
                  ELSE
                     CALL EXPLIM (NDIM, ALMESH, RNUM)
                  END IF
               END IF

               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'horizontal minimum', RNUM(KLFT), RNUM(KLFT), *100)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'horizontal maximum', RNUM(KRGT), RNUM(KRGT), *100)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'vertical minimum', RNUM(KBOT), RNUM(KBOT), *100)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'vertical maximum', RNUM(KTOP), RNUM(KTOP), *100)
            END IF
            IF (RNUM(KLFT) .GT. RNUM(KRGT)) THEN
               X = RNUM(KLFT)
               RNUM(KLFT) = RNUM(KRGT)
               RNUM(KRGT) = X
            END IF
            IF (RNUM(KBOT) .GT. RNUM(KTOP)) THEN
               X = RNUM(KBOT)
               RNUM(KBOT) = RNUM(KTOP)
               RNUM(KTOP) = X
            END IF
            IF ((RNUM(KLFT) .GE. RNUM(KRGT))
     &         .OR. (RNUM(KBOT) .GE. RNUM(KTOP))) THEN
               CALL PRTERR ('CMDERR', 'No limits are defined by ZOOM')
               GOTO 100
            END IF
            CALL FFADDR (RNUM(KLFT), INLINE)
            CALL FFADDR (RNUM(KRGT), INLINE)
            CALL FFADDR (RNUM(KBOT), INLINE)
            CALL FFADDR (RNUM(KTOP), INLINE)

            MSCTYP = 'ZOOM'
            CALL CPYREA (KTOP, RNUM, RDMESH)
            NZMON = .FALSE.

C ZOOM RADIUS COMMAND " ZOOM RADIUS XZM YZM (ZZM) RADZM "

         ELSE IF (MATSTR (WORD, 'RADIUS', 1)) THEN
            CALL FFADDC (WORD, INLINE)
C -- MAKE SURE THERE ARE VALUES IN THE INPUT
            IF(IS3DIM) THEN
               CALL FFNEED(IFLD, INTYP, 'R', 3,
     &                     'x,y,z center point',*100)
            ELSE
               CALL FFNEED(IFLD, INTYP, 'R', 2,
     &                     'x,y center point',*100)
            END IF
C -- GET X AND Y COORDINATES OF ZOOM CENTER
            CALL FFREAL (IFLD, INTYP, RFIELD, 'x', 0.0, XZM, *100)
            CALL FFADDR (XZM, INLINE)
            CALL FFREAL (IFLD, INTYP, RFIELD, 'y', 0.0, YZM, *100)
            CALL FFADDR (YZM, INLINE)
C -- IF 3D, GET THE Z COORDINATE
            IF(IS3DIM) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD, 'z', 0.0, ZZM, *100)
               CALL FFADDR (ZZM, INLINE)
            END IF
C -- GET RADIUS OF ZOOM FIELD
            CALL FFNEED(IFLD, INTYP, 'R', 1,
     &                  'zoom window radius',*100)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'radius', 0.0, RADZM, *100)
            CALL FFADDR (RADZM, INLINE)
            NODEZM = 0
            NZMON = .TRUE.
            MSCTYP = 'ZOOM'

C ZOOM NODE COMMAND " ZOOM NODE NODEID RADIUS" or "ZOOM NODE CURSOR"

         ELSE IF (MATSTR (WORD, 'NODE', 1)) THEN
            CALL FFADDC (WORD, INLINE)
C -- SEE IF CURSOR OPTION IS SELECTED
            IF (FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)) THEN
               CALL QNPICK('DISPLAYED', LDUM1, LDUM2, A, KXN,
     &                         KYN, KZN, KHIDEN, KNPSUR)
               IF(IS3DIM) THEN

C -- 3D PICK OF TWO NODES
                  CALL PICKN3('Select center node', NUMNPF, A(KXN),
     &               A(KYN), A(KZN), A(KHIDEN), .TRUE., NODEZM, *100)
                  CALL PICKN3('Select radius node', NUMNPF, A(KXN),
     &               A(KYN), A(KZN), A(KHIDEN), .FALSE., NRAD, *100)
               ELSE
C -- 2D PICK OF TWO NODES
                  CALL PICKN2('Select center node', NUMNPF, A(KXN),
     &                 A(KYN), .TRUE., NODEZM, *100)
                  CALL PICKN2('Select radius node', NUMNPF, A(KXN),
     &                 A(KYN), .FALSE., NRAD, *100)
               END IF

C --  GET DISTANCE FROM CENTER TO RADIUS POINT
               CALL GETDST(NODEZM, NRAD, A(KXN), A(KYN), A(KZN),
     &                        RADZM)

            ELSE
C -- INPUT NODE ID AND RADIUS
               CALL FFNEED(IFLD, INTYP, 'I', 1,
     &                     'node id',*100)
               CALL FFINTG (IFLD, INTYP, IFIELD, 'node id', 0,
     &                      INP, *100)
               CALL FFADDI (INP, INLINE)

C ... Convert global node id to local node offset
               NODEZM = locint(inp, numnp, mapnd)

               CALL FFNEED(IFLD, INTYP, 'R', 1,
     &                     'radius',*100)
               CALL FFREAL (IFLD, INTYP, RFIELD, 'radius', 0.0,
     &                      RADZM, *100)
               CALL FFADDR (RADZM, INLINE)

               XZM = 0.0
               YZM = 0.0
               ZZM = 0.0

            END IF
C -- SET OTHER VARIABLES
            NZMON = .TRUE.
            MSCTYP = 'ZOOM'

         ELSE
            IF (.NOT. IS3DIM) THEN
               CALL PRTERR ('CMDERR', 'Zoom options include'
     &            // ' "CURSOR", "TRANSLATE", "EACH", and "MESH"')
            ELSE
               CALL PRTERR ('CMDERR', 'Zoom options include'
     &            // ' "CURSOR", "TRANSLATE", "EACH", and "ROTATION"')
            END IF
            GOTO 100
         END IF

         NEWZM = .TRUE.
         IF (.NOT. SETTIC) TICMSH = 0.0

      ELSE IF (VERB .EQ. 'SCALE') THEN
         CALL FFADDC (VERB, INLINE)
         call prterr ('CMDREQ', 'Please use the ZOOM command')

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'EACH', WORD)
         IF (MATSTR (WORD, 'ROTATION', 3)) THEN
            IF (.NOT. IS3DIM) THEN
               CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
               GOTO 100
            END IF
            CALL FFADDC ('ROTATION', INLINE)
            MSCTYP = 'ROTATION'
         ELSE IF (MATSTR (WORD, 'EACH', 1)) THEN
            CALL FFADDC ('EACH', INLINE)
            MSCTYP = 'EACH'
         ELSE
            CALL PRTERR ('CMDERR', 'Expected "SET" or "EACH"')
            GOTO 100
         END IF

         NEWZM = .TRUE.
         IF (.NOT. SETTIC) TICMSH = 0.0

      ELSE IF (VERB .EQ. 'TICK') THEN
         CALL FFADDC (VERB, INLINE)

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'tick interval', 0.0, TICMSH, *100)
         CALL FFADDR (TICMSH, INLINE)

         SETTIC = .TRUE.

      ELSE IF (VERB .EQ. 'SQUARE') THEN
         CALL FFADDC (VERB, INLINE)

         CALL FFONOF (IFLD, INTYP, CFIELD, SQMESH, *100)
         CALL FFADDO (SQMESH, INLINE)
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
