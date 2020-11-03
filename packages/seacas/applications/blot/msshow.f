C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSSHOW (SHOTYP, NAMECO, NAMES, IDELB, IELBST,
     &   IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL)
C=======================================================================

C   --*** MSSHOW *** (MESH) Display mesh parameter information
C   --   Written by Amy Gilkey - revised 05/26/88
C   --
C   --MSSHOW displays the mesh plot parameters.
C   --
C   --The SHOW options with the items they display are:
C   --   PLOT     - display mode and plot variables
C   --   HARDCOPY -
C   --   VIEW     -
C   --   EMPTY    -
C   --   DEFORM   -
C   --   MLINES   - mesh lines displayed
C   --   BOUNDARY -
C   --   NSETS    - display the selected node sets
C   --   SSETS    - display the selected side sets
C   --   VISIBLE  - ON and selected element blocks
C   --   BLOCKS   -
C   --   DEATH    - element birth/death variable
C   --   XSYM     - vertical and horizontal window symmetry axis
C   --   YSYM     -    or windows defined for non-symmetric views
C   --   XVIEW    -
C   --   YVIEW    -
C   --   MULTTIME - display different time on each view flag
C   --   MAGNIFY  - displacement magnification factor and default
C   --   HIDDEN   - (3D) hidden line option
C   --   ZOOM     - total window coordinates and zoom window coordinates
C   --   TICK     - mesh axis tick interval
C   --   SQUARE   - square/non-square window scaling flag
C   --   ROTATE   - rotation matrix and eye position
C   --   EYE      -
C   --   CENTER   - center of rotation
C   --   CUT      - cutting plane
C   --   VECSCL   - vector/symbol scale factor
C   --   DEADNODE - dead node display option
C   --   BLKCOL string1 string2 ... stringi
C   --      where
C   --          stringi = block_id1 block_id2 ... block_idj block_coli
C   --                             Sets the color of the element blocks
C   --                             specified.  block_idi is the integer
C   --                             identifier of a block.  block_coli
C   --                             identifies the color to be given
C   --                             to the element blocks specified in
C   --                             stingi.  The color must be identified
C   --                             by its name, although it may be
C   --                             abbreviated to uniqueness within the
C   --                             list of colors.  The available colors
C   --                             are black, white, red, green, yellow,
C   --                             blue, cyan, and magenta.
C   --   LINETHIC - mesh line thickness specification
C   --   SPHERE or FSPHERE - specify if elements are to be displayed
C   --                       as spheres or filled spheres.
C   --
C   --Parameters:
C   --   SHOTYP - IN - the expanded SHOW option string
C   --   NAMECO - IN - the coordinate names
C   --   NAMES - IN - the variable names
C   --   IDELB - IN - the element block ID array
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   IDNPS - IN - the node set ID for each set
C   --   ISSNPS - IN - the indices of the selected node sets
C   --   IDESS - IN - the side set ID for each set
C   --   ISSESS - IN - the indices of the selected side sets
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --
C   --Common Variables:
C   --   Uses NDIM, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses VECSCL of /ETCOPT/
C   --   Uses DEFOK, DEFFAC, DFAC of /DEFORM/
C   --   Uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, NNPSET, NESSET,
C   --      IHIDOP, NALVAR, DEADNP of /MSHOPT/
C   --   Uses MULTIM, XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM of /VIEWS/
C   --   Uses ALMESH, ZMMESH, TICMSH, MSCTYP, SQMESH of /MSHLIM/
C   --   Uses ROTMAT, ROTCEN, EYE of /ROTOPT/
C   --   Uses ISCUT, CUTPLA of /CUTOPT/
C   --   Uses COLLST of /CLST/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'etcopt.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'cutopt.blk'
      include 'linthc.blk'
      include 'sphele.blk'
      include 'nodzom.blk'
      include 'light.blk'

      CHARACTER*(*) SHOTYP
      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IDELB(NELBLK)
      INTEGER IELBST(NELBLK)
      INTEGER IDNPS(*)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDESS(*)
      INTEGER BLKCOL(0:NELBLK)

      INTEGER MAINVW, NDEFVW, IXVW
      LOGICAL DEF
      CHARACTER*132 STRING
      CHARACTER*20 STR20
      CHARACTER CH
      CHARACTER*(MXSTLN) S1
      CHARACTER*(MXNAME) NAM(2)
      CHARACTER*132 ELBLIN
      CHARACTER*20 RSTR(9)
      REAL RNUM(9)
      CHARACTER*(MXSTLN) THKSPC

      include 'cmap-lst.blk'

C *** Display mode control and multiple views mode control ***

      IF ((SHOTYP .EQ. 'PLOT') .OR. (SHOTYP .EQ. 'HARDCOPY')
     &   .OR. (SHOTYP .EQ. 'VIEW') .OR. (SHOTYP .EQ. 'EMPTY')
     &   .OR. (SHOTYP .EQ. 'DEFORM')
     &   .OR. (SHOTYP .EQ. 'NUMBER')
     &   .OR. (SHOTYP .EQ. 'NSETS') .OR. (SHOTYP .EQ. 'SSETS')) THEN

C      --Find the "main" view on symmetric views
         MAIN = MAINVW ()

         DEF = (DFAC .NE. 0.0)
         NVIEWS = NDEFVW (.TRUE.)
         DO 100 IVW = 1, NVIEWS
            IVIEW = IXVW (.TRUE., IVW)

            IF (MSHDEF(IVIEW) .EQ. 'EMPTY') THEN
               STRING = 'Empty view'
            ELSE
               IF (DEF .AND. (MSHDEF(IVIEW) .EQ. 'DEFORM')) THEN
                  STRING = 'Deformed mesh'
               ELSE IF (MSHDEF(IVIEW) .EQ. 'UNDEFORM') THEN
                  STRING = 'Undeformed mesh'
               ELSE
                  STRING = 'Mesh'
               END IF
            END IF
            LSTR = LENSTR (STRING)

            IF (MSHDEF(IVIEW) .EQ. 'EMPTY') THEN
               CONTINUE
            ELSE IF (MSHNUM(IVIEW) .EQ. 'NONE') THEN
               STRING(LSTR+1:) = ' with no numbering'
            ELSE IF (MSHNUM(IVIEW) .EQ. 'NODE') THEN
               STRING(LSTR+1:) = ' with numbered nodes'
            ELSE IF (MSHNUM(IVIEW) .EQ. 'ELEMENT') THEN
               STRING(LSTR+1:) = ' with numbered elements'
            ELSE IF (MSHNUM(IVIEW) .EQ. 'ALL') THEN
               STRING(LSTR+1:) = ' with numbered nodes and elements'
            END IF
            LSTR = LENSTR (STRING)

            IF (MSHDEF(IVIEW) .EQ. 'EMPTY') THEN
               CONTINUE
            ELSE IF ((NNPSET(IVIEW) .GT. 0)
     &         .AND. (NESSET(IVIEW) .GT. 0)) THEN
               STRING(LSTR+1:) =
     &            ' with node sets and side sets'
            ELSE IF (NNPSET(IVIEW) .GT. 0) THEN
               STRING(LSTR+1:) =
     &            ' with node sets'
            ELSE IF (NESSET(IVIEW) .GT. 0) THEN
               STRING(LSTR+1:) =
     &            ' with side sets'
            END IF
            LSTR = LENSTR (STRING)

            IF (NVIEWS .GT. 1) THEN
               WRITE (CH, '(I1)') IVIEW
               IF (IVIEW .EQ. MAIN) THEN
                  WRITE (*, 10150) 'VIEW ', CH, ':* ', STRING(:LSTR)
               ELSE
                  WRITE (*, 10150) 'VIEW ', CH, ':  ', STRING(:LSTR)
               END IF
            ELSE
               WRITE (*, 10150) STRING(:LSTR)
            END IF
  100    CONTINUE

      ELSE IF ((SHOTYP .EQ. 'MLINES')
     &   .OR. (SHOTYP .EQ. 'BOUNDARY')) THEN

C      --Find the "main" view on symmetric views
         MAIN = MAINVW ()

         NVIEWS = NDEFVW (.TRUE.)
         DO 110 IVW = 1, NVIEWS
            IVIEW = IXVW (.TRUE., IVW)

            IF (MSHDEF(IVIEW) .EQ. 'EMPTY') THEN
               STRING = 'Empty view'
            ELSE IF (MSHLIN(IVIEW) .EQ. MSHNON) THEN
               STRING = 'Display mesh without boundaries'
            ELSE IF (MSHLIN(IVIEW) .EQ. MSHBOR) THEN
               STRING = 'Display mesh boundaries only'
            ELSE IF (MSHLIN(IVIEW) .EQ. MSHDIV) THEN
               STRING = 'Display element block divisions only'
            ELSE IF (MSHLIN(IVIEW) .EQ. MSHSEL) THEN
               IF (IABS (MLNTYP(1,IVIEW)) .LE. 1) THEN
                  STRING = 'Display mesh of elements'
     &               //' in selected element blocks'
               ELSE
                  STRING = 'Display mesh of elements'
     &               //' in selected element blocks'
     &               // ' with dotted lines'
               END IF
            ELSE IF (MSHLIN(IVIEW) .EQ. MSHALL) THEN
               IF (IABS (MLNTYP(1,IVIEW)) .LE. 1) THEN
                  STRING = 'Display mesh of all elements'
               ELSE
                  STRING = 'Display mesh of all elements'
     &               // ' with dotted lines'
               END IF
            END IF
            LSTR = LENSTR (STRING)

            IF (NVIEWS .GT. 1) THEN
               WRITE (CH, '(I1)') IVIEW
               IF (IVIEW .EQ. MAIN) THEN
                  WRITE (*, 10150) 'VIEW ', CH, ':* ', STRING(:LSTR)
               ELSE
                  WRITE (*, 10150) 'VIEW ', CH, ':  ', STRING(:LSTR)
               END IF
            ELSE
               WRITE (*, 10150) STRING(:LSTR)
            END IF

            IF ((MLNTYP(-1,IVIEW) .LT. 0)
     &         .AND. (MSHLIN(IVIEW) .GE. MSHBOR)) THEN
               WRITE (*, 10150) '   (mesh boundaries drawn in black)'
            END IF
  110    CONTINUE

      ELSE IF ((SHOTYP .EQ. 'NSETS') .OR. (SHOTYP .EQ. 'SSETS')) THEN
C      --Find the "main" view on symmetric views
         MAIN = MAINVW ()

         NVIEWS = NDEFVW (.FALSE.)
         DO 120 IVW = 1, NVIEWS
            IVIEW = IXVW (.FALSE., IVW)

            WRITE (STR20, 10000, IOSTAT=IDUM)
     &         NNPSET(IVIEW), NUMNPS
10000        FORMAT (I5, ' of ', I5)
            CALL SQZSTR (STR20, LSTR)
            STRING = 'Selected node sets:  ' // STR20(:LSTR)
            LSTR = LENSTR (STRING)
            IF (NVIEWS .GT. 1) THEN
               WRITE (CH, '(I1)') IVIEW
               IF (IVIEW .EQ. MAIN) THEN
                  WRITE (*, 10150) 'VIEW ', CH, ':* ', STRING(:LSTR)
               ELSE
                  WRITE (*, 10150) 'VIEW ', CH, ':  ', STRING(:LSTR)
               END IF
            ELSE
               WRITE (*, 10150) STRING(:LSTR)
            END IF
            IF (NNPSET(IVIEW) .GT. 0) WRITE (*, 10010)
     &         (IDNPS(ISSNPS(IX,IVIEW)), IX=1,NNPSET(IVIEW))
10010        FORMAT (4X, 12I6)

            WRITE (STR20, 10000, IOSTAT=IDUM)
     &         NESSET(IVIEW), NUMESS
            CALL SQZSTR (STR20, LSTR)
            STRING = 'Selected side sets:  ' // STR20(:LSTR)
            LSTR = LENSTR (STRING)
            IF (NVIEWS .GT. 1) THEN
               WRITE (CH, '(I1)') IVIEW
               IF (IVIEW .EQ. MAIN) THEN
                  WRITE (*, 10150) 'VIEW ', CH, ':* ', STRING(:LSTR)
               ELSE
                  WRITE (*, 10150) 'VIEW ', CH, ':  ', STRING(:LSTR)
               END IF
            ELSE
               WRITE (*, 10150) STRING(:LSTR)
            END IF
            IF (NESSET(IVIEW) .GT. 0) WRITE (*, 10010)
     &         (IDESS(ISSESS(IX,IVIEW)), IX=1,NESSET(IVIEW))
  120    CONTINUE

C *** Active element control ***

      ELSE IF ((SHOTYP .EQ. 'VISIBLE')
     &   .OR. (SHOTYP .EQ. 'BLOCKS') .OR. (SHOTYP .EQ. 'MATERIAL')) THEN
         CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
         IF (IS3DIM) THEN
            WRITE (STRING, 10020, IOSTAT=IDUM) NUMON, NELBLK
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10150) 'Visible element blocks:  ', STRING(:LSTR)
            N = 0
            ELBLIN = ' '
            DO 130 IELB = 1, NELBLK
               IF (IELBST(IELB) .GE. 0) THEN
                  N = N + 1
                  WRITE (ELBLIN((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     &               IDELB(IELB)
                  IF (N .GE. 12) THEN
                     WRITE (*, 10150) ELBLIN(1:LENSTR(ELBLIN))
                     N = 0
                     ELBLIN = ' '
                  END IF
               END IF
  130       CONTINUE
            IF (N .GT. 0) THEN
               WRITE (*, 10150) ELBLIN(1:LENSTR(ELBLIN))
            END IF
         END IF

         WRITE (STRING, 10020, IOSTAT=IDUM) NUMSEL, NELBLK
10020     FORMAT (I3, ' of ', I3)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10150) 'Active element blocks:  ', STRING(:LSTR)
         N = 0
         ELBLIN = ' '
         DO 140 IELB = 1, NELBLK
            IF (IELBST(IELB) .GT. 0) THEN
               N = N + 1
               WRITE (ELBLIN((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     &            IDELB(IELB)
               IF (N .GE. 12) THEN
                  WRITE (*, 10150) ELBLIN(1:LENSTR(ELBLIN))
                  N = 0
                  ELBLIN = ' '
               END IF
            END IF
  140    CONTINUE
         IF (N .GT. 0) THEN
            WRITE (*, 10150) ELBLIN(1:LENSTR(ELBLIN))
         END IF

      ELSE IF (SHOTYP .EQ. 'DEATH') THEN
         IF (NALVAR .GT. 0) THEN
           CALL NUMSTR1(4, ALIVAL, RSTR(1), LSTR)
            WRITE (*, 10150) 'Birth/Death variable: ',
     &       NAMES(NALVAR)(:LENSTR(NAMES(NALVAR))),
     *       ', Alive value = ', RSTR(1)(:LSTR)
         ELSE
            WRITE (*, 10150) 'Birth/Death not defined'
         END IF

C *** Multiple views control ***

      ELSE IF ((SHOTYP .EQ. 'XSYM') .OR. (SHOTYP .EQ. 'YSYM')
     &   .OR. (SHOTYP .EQ. 'XVIEW') .OR. (SHOTYP .EQ. 'YVIEW')) THEN
         IF (NDEFVW (.TRUE.) .GT. 1) THEN
            IF (XISSYM) THEN
               S1 = 'right'
               IF (LFTSYM) S1 = 'left'
               CALL NUMSTR1 (4, XAXSYM, RSTR(1), LSTR)
               WRITE (*, 10030) 'vertical', S1(:LENSTR(S1)),
     &            RSTR(1)(:LSTR)
            ELSE IF (MSHDEF(1) .NE. 'NONE') THEN
               WRITE (*, 10040) 'vertically divided'
            END IF
            IF (YISSYM) THEN
               S1 = 'top'
               IF (BOTSYM) S1 = 'bottom'
               CALL NUMSTR1 (4, YAXSYM, RSTR(1), LSTR)
               WRITE (*, 10030) 'horizontal', S1(:LENSTR(S1)),
     &            RSTR(1)(:LSTR)
            ELSE IF (MSHDEF(4) .NE. 'NONE') THEN
               WRITE (*, 10040) 'horizontally divided'
            END IF
10030        FORMAT (' Symmetric Views: ', A, ' axis on ', A, ' at ', A)
10040        FORMAT (' Non-symmetric Views: ', A)
         ELSE
            WRITE (*, 10150) 'Single view defined'
         END IF

      ELSE IF (SHOTYP .EQ. 'MULTTIME') THEN
         IF (MULTIM) THEN
            WRITE (*, 10150) 'Different time on each view'
         ELSE
            WRITE (*, 10150) 'Single times for all views'
         END IF

C *** Mesh control ***

      ELSE IF (SHOTYP .EQ. 'MAGNIFY') THEN
         IF (.NOT. DEFOK) THEN
            WRITE (*, 10150) 'No displacement variables found'
         ELSE IF (DEFFAC .LT. 0.0) THEN
            CALL NUMSTR1(4, DFAC, RSTR(1), L1)
            WRITE (*, 10150) 'Displacement magnification factor = ',
     &         RSTR(1)(:L1), ', default not calculated'
         ELSE
            CALL NUMSTR1(4, DFAC, RSTR(1), L1)
            CALL NUMSTR1(4, DEFFAC, RSTR(2), L2)
            WRITE (*, 10150) 'Displacement magnification factor = ',
     &         RSTR(1)(:L1), ', Calculated = ', RSTR(2)(:L2)
         END IF

       ELSE IF (SHOTYP .EQ. 'HIDDEN') THEN
         IF (.NOT. IS3DIM) GOTO 180

         IF (IHIDOP .EQ. 0) THEN
           WRITE (*, 10150) 'Display hidden surfaces and lines (0)'
         ELSE IF (IHIDOP .EQ. 1) THEN
           WRITE (*, 10150) 'Remove surfaces facing away (1)'
         ELSE IF (IHIDOP .EQ. 2) THEN
           WRITE (*, 10150) 'Remove surfaces with any hidden nodes (2)'
         ELSE IF (IHIDOP .EQ. 3) THEN
           WRITE (*, 10150) 'Remove partial hidden lines (3)'
         ELSE IF (IHIDOP .EQ. 4) THEN
           WRITE (*, 10150) 'Draw surfaces in reverse order (4)'
         ELSE IF (IHIDOP .EQ. 5) THEN
           WRITE (*, 10150) 'Draw surfaces in reverse order, shade (5)'
         ELSE IF (IHIDOP .GE. 6) THEN
           WRITE (*, 10150) 'Write polygons to rayshade file (6)'
         END IF

      ELSE IF ((SHOTYP .EQ. 'ZOOM') .or. (shotyp .eq. 'SCALE')) THEN
         IF (DFAC .EQ. 0.0) THEN
            CALL CPYREA (2*NDIM, UNMESH, RNUM)
         ELSE
            CALL CPYREA (2*NDIM, ALMESH, RNUM)
         END IF
         CALL NUMSTR (2*NDIM, 4, RNUM, RSTR, LSTR)
         L = 0
         DO 150 I = 1, NDIM
            L = MAX (L, LENSTR(NAMECO(I)))
  150    CONTINUE
         WRITE (STRING, 10050, IOSTAT=IDUM)
     &      'Mesh Limits', (NAMECO(I)(:L),
     &      RSTR(2*I-1)(:LSTR), RSTR(2*I)(:LSTR), I=1,NDIM)
10050     FORMAT (A, ':', 3 (3 (' ', A), :, ','))
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10150) STRING(:LSTR)

         IF (MSCTYP .NE. 'EACH') THEN
            IF (IS3DIM) THEN
               NAM(1) = 'H'
               NAM(2) = 'V'
            ELSE
               NAM(1) = NAMECO(1)
               NAM(2) = NAMECO(2)
            END IF
            CALL CPYREA (KTOP, ZMMESH, RNUM)
            CALL NUMSTR (KTOP, 4, RNUM, RSTR, LSTR)
            WRITE (STRING, 10050, IOSTAT=IDUM)
     &         'Zoom Limits', (NAM(I)(:L),
     &         RSTR(2*I-1)(:LSTR), RSTR(2*I)(:LSTR), I=1,2)
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10150) STRING(:LSTR)
         END IF

         IF (MSCTYP .EQ. 'ROTATION') THEN
            WRITE (*, 10150) 'Limits are adjusted for each rotation'
         ELSE IF (MSCTYP .EQ. 'MESH') THEN
            WRITE (*, 10150) 'Limits are the expanded mesh limits'
         ELSE IF (MSCTYP .EQ. 'EACH') THEN
            WRITE (*, 10150) 'Limits are adjusted for each plot'
         ELSE IF (NZMON) THEN
C -- WRITE MESH CENTERING INFORMATION
            IF(NODEZM .EQ. 0) THEN
               RNUM(1) = XZM
               RNUM(2) = YZM

               IF(IS3DIM) THEN
                  RNUM(3) = ZZM
               END IF
               CALL NUMSTR(NDIM, 4, RNUM, RSTR, LSTR)
               IF(IS3DIM) THEN
                   STRING = 'Mesh Center At: ' // 'X ' // RSTR(1)
     &                 // 'Y ' // RSTR(2) // 'Z ' // RSTR(3)
               ELSE
                   STRING = 'Mesh Center At: ' // 'X ' // RSTR(1)
     &                 // 'Y ' // RSTR(2)
               END IF
            ELSE
               WRITE(RSTR(1),'(I10)')NODEZM
               STRING = 'Mesh Center At Point: ' // rstr(1)
            END IF
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10150) STRING(:LSTR)
C -- WRITE ZOOM RADIUS INFORMATION
            CALL NUMSTR1(4, RADZM, RSTR, LSTR)
            STRING = 'Mesh Zoom Radius: ' // RSTR(1)
            CALL SQZSTR (STRING, LSTR)
            WRITE (*, 10150) STRING(:LSTR)
         END IF

      ELSE IF (SHOTYP .EQ. 'TICK') THEN
         IF (TICMSH .EQ. 0.0) THEN
            WRITE (*, 10150)
     &         'Mesh axis tick interval automatically scaled'
         ELSE
            CALL NUMSTR1(4, TICMSH, RSTR(1), LSTR)
            WRITE (*, 10150) 'Mesh axis tick interval = ',
     &         RSTR(1)(:LSTR)
         END IF

      ELSE IF (SHOTYP .EQ. 'SQUARE') THEN
         IF (SQMESH) THEN
            WRITE (*, 10150) 'Window will be square'
         ELSE
            WRITE (*, 10150) 'Window dimensions based on mesh'
         END IF

      ELSE IF ((SHOTYP .EQ. 'ROTATE') .OR. (SHOTYP .EQ. 'EYE')) THEN
         IF (.NOT. IS3DIM) GOTO 180

         CALL NUMSTR (3, 4, EYE, RSTR, LSTR)
         WRITE (STRING, 10140) 'Eye at ', (RSTR(I)(:LSTR), I=1,3)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10150) STRING(:LSTR)

      ELSE IF (SHOTYP .EQ. 'LIGHTS' .OR. SHOTYP .EQ. 'AMBIENT') THEN
         IF (.NOT. IS3DIM) GOTO 180
         if (nlit .eq. 0) then
           write (*, 10150)
     *       'No Light vectors defined'
         else
           WRITE (*, 10150)
     *       'Light vectors and brightness (screen coordinates) '
           do 155 ilit = 1, NLIT
             WRITE (*, 10170) ilit, (LITE(I,ILIT),I=1,4)
 155       continue
         end if
         call numstr1(4, AMBIENT, RSTR, LSTR)
         write (STRING, 10140) 'Ambient Light Intensity is ',
     *     RSTR(1)(:LSTR)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10150) STRING(:LSTR)

      ELSE IF (SHOTYP .EQ. 'CENTER') THEN
         IF (.NOT. IS3DIM) GOTO 180

         CALL NUMSTR (3, 4, ROTCEN, RSTR, LSTR)
         WRITE (STRING, 10140) 'Center of rotation = ',
     &      (RSTR(I)(:LSTR), I=1,3)
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10150) STRING(:LSTR)

      ELSE IF (SHOTYP .EQ. 'CUT') THEN
         IF (.NOT. IS3DIM) GOTO 180

         IF (ISCUT) THEN
            WRITE (*, 10150) 'Cutting Plane at point:'
            CALL NUMSTR (3, 4, CUTPT, RSTR, LSTR)
            WRITE (*, 10140) '    ', (RSTR(I)(:LSTR), I=1,3)
            WRITE (*, 10150) 'Normal to Cutting Plane :'
            CALL NUMSTR (3, 4, CUTNRM, RSTR, LSTR)
            WRITE (*, 10140) '    ', (RSTR(I)(:LSTR), I=1,3)
         ELSE
            WRITE (*, 10150) 'No Cutting Plane is defined'
         END IF

      ELSE IF (SHOTYP .EQ. 'LINETHIC') THEN
         WRITE (*, *)
         WRITE (*, *)
         WRITE (*, 10060)
10060     FORMAT (1X, 'Line thickness (specified as a real value ',
     &      'from 0. to 1000.)' /)
         IF (MSHBND .LT. ((THKNSS(3)+THKNSS(2))/2)) THEN
            THKSPC = 'THIN'
         ELSE IF (MSHBND .LT. ((THKNSS(2)+THKNSS(1))/2)) THEN
            THKSPC = 'MEDIUM'
         ELSE
            THKSPC = 'THICK'
         ENDIF
         LEN = LENSTR (THKSPC)
         WRITE (*, 10070) MSHBND, THKSPC(1:LEN)
10070     FORMAT (5X, 'Mesh Boundary -- ', F5.0, ' (',A,')')
         IF (BLKBND .LT. ((THKNSS(3)+THKNSS(2))/2)) THEN
            THKSPC = 'THIN'
         ELSE IF (BLKBND .LT. ((THKNSS(2)+THKNSS(1))/2)) THEN
            THKSPC = 'MEDIUM'
         ELSE
            THKSPC = 'THICK'
         ENDIF
         LEN = LENSTR (THKSPC)
         WRITE (*, 10080) BLKBND, THKSPC(1:LEN)
10080     FORMAT (5X, 'Element Block Boundaries -- ', F5.0, ' (',A,')')
         IF (ELEBND .LT. ((THKNSS(3)+THKNSS(2))/2)) THEN
            THKSPC = 'THIN'
         ELSE IF (ELEBND .LT. ((THKNSS(2)+THKNSS(1))/2)) THEN
            THKSPC = 'MEDIUM'
         ELSE
            THKSPC = 'THICK'
         ENDIF
         LEN = LENSTR (THKSPC)
         WRITE (*, 10090) ELEBND, THKSPC(1:LEN)
10090     FORMAT (5X, 'Element Boundaries -- ', F5.0, ' (',A,')')

      ELSE IF (SHOTYP .EQ. 'SPHERE'
     &   .OR.  SHOTYP .EQ. 'FSPHERE') THEN

         WRITE (*, *)
         WRITE (*, *)
         IF (SPHPLT .EQ. 0) THEN
            WRITE (*, 10100)
10100        FORMAT (1X, 'Elements will be plotted normally')
         ELSE IF (SPHPLT .GE. 1) THEN
            WRITE (*, 10110) SPHPLT
10110        FORMAT (1X,
     *        'Elements will be plotted as outlined spheres with ',i2,
     *        ' segments.')
         ELSE
            WRITE (*, 10120) -SPHPLT
10120        FORMAT (1x,
     *        'Elements will be plotted as filled spheres with ',i2,
     *        ' segments.')
         ENDIF
         if (SPHPLT .ne. 0 .and. defrad .ne. 1.0) then
           WRITE (*, 10121) DEFRAD
10121      FORMAT (1x,'   The default sphere radius is ',F10.3)
         ENDIF
C *** Display options ***

      ELSE IF (SHOTYP .EQ. 'VECSCL') THEN
         CALL NUMSTR1(4, VECSCL, RSTR, LSTR)
         WRITE (*, 10150) 'Vector/symbol scale factor = ',
     &      RSTR(1)(:LSTR)

      ELSE IF (SHOTYP .EQ. 'DEADNODE') THEN
         IF (DEADNP) THEN
            WRITE (*, 10150) 'Display dead nodes'
         ELSE
            WRITE (*, 10150) 'Do not display dead nodes'
         END IF

      ELSE IF (SHOTYP .EQ. 'BLKCOL') THEN

         WRITE (*, *)
         IF (BLKCOL(0) .EQ. -1) THEN
            WRITE (*, 10150) 'User specified block colors will not be',
     &         ' used on plots.'
         ELSE
            WRITE (*, 10150) 'User specified block colors will be',
     &         ' used on plots.'
         ENDIF
         WRITE (*, *)
         WRITE (*, 10150) 'User specified colors'
         DO 170 I = 1, NELBLK
            IF (BLKCOL(I) .NE. -2) THEN
               WRITE (*, 10130) I, COLLST(BLKCOL(I)+2)
10130           FORMAT (1X, 'Block ', I5, ' is assigned ', A)
            ENDIF
  170    CONTINUE

      END IF

  180 CONTINUE
      RETURN
10140  FORMAT (A, 3 (A, ' '))
10150  FORMAT (1X, 10A)
10170  FORMAT ('(',I2,') X: ',1pe10.3,' Y: ',1pe10.3,' Z: ',1pe10.3,
     *   ' Br: ',1pe10.3)
      END
