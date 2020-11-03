C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSCOMD (A, CURPRO, MSHTYP, INLINE, INVERB, IFLD,
     &           INTYP, CFIELD, IFIELD, RFIELD, NAMECO, NAMES,
     &           IDELB, NEWELB, IELBST, NEWMOD,  IDNPS, ISSNPS,
     &           IDESS, ISSESS, BLKCOL, SHDCOL, ISHDCL)
C=======================================================================

C   --*** MSCOMD *** (MESH) Process mesh commands
C   --   Modified by John Glick - 1/13/89
C   --   Written by Amy Gilkey - revised 05/31/88
C   --   Dennis Flanagan, 11/18/82
C   --
C   --MSCOMD interprets mesh commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Display mode control
C   --   EMPTY                     Set empty window
C   --   DEFORM   {ON/OFF}         Set deformed versus undeformed mesh
C   --   NUMBER   {NODES/ELEMENTS/ALL} Set to number nodes or elements
C   --   MLINES   {ON/OFF}         Set on/off mesh overlay option
C   --   BOUNDARY {ON/OFF}         Set on/off mesh and element block boundary
C   --                             option
C   --   NSETS    {id1,...}        Select node sets
C   --   SSETS    {id1,...}        Select side sets
C   --
C   --Active element control
C   --   VISIBLE  {element blocks} Set ON (displayed) element blocks (3D)
C   --   BLOCKS   {element blocks} Set selected element blocks
C   --   DEATH    {ON,variable}    Enable element birth/death
C   --   DEATH    {OFF}            Disable element birth/death
C   --
C   --Element block color control
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
C   --   MATCOL string1 string2 ... stringi
C   --                             Performs the same function as the
C   --                             BLKCOL command.
C   --
C   --Multiple views control
C   --   VIEW     {IVIEW,command}  Set view-dependent parameter mode for view
C   --   XSYM     {LEFT/RIGHT,xaxsym} Define vertical symmetry axis
C   --   XSYM     {OFF}            Delete vertically divided views
C   --   XVIEW                     Define vertical non-symmetric views
C   --   XVIEW    {OFF}            Same as XSYM OFF
C   --   YSYM     {BOTTOM/TOP,yaxsym} Define horizontal symmetry axis
C   --   YSYM     {OFF}            Delete horizontally divided views
C   --   YVIEW                     Define horizontal non-symmetric views
C   --   YVIEW    {OFF}            Same as YSYM OFF
C   --   MULTTIME {ON/OFF}         Define multiple time views
C   --
C   --Mesh control
C   --   MAGNIFY  {dfac}           Set displacement magnification factor
C   --   HIDDEN   {0/1/2/3}        Set hidden line removal option (3D)
C   --   ZOOM     {left,right,bottom,top} Set output window min/max
C   --   ZOOM     {EACH/SET}       Set window scaling
C   --   SQUARE   {ON/OFF}         Set square/non-square window scaling
C   --   TICK     {ticmsh}         Set mesh axis tick interval
C   --   ROTATE   {X/Y/Z/RESET,value} Rotate the mesh (cumulative) (3D)
C   --   EYE      {x,y,z}          Rotate the mesh to given position (3D)
C   --   CENTER   {x,y,z}          Change the center of rotation (3D)
C   --   CENTER   {RESET}          Reset the center of rotation (3D)
C   --   CUT      {x1,y1,x2,y2,x3,y3} Specify a cutting plane (3D)
C   --   CUT      {OFF}            Delete the cutting plane (3D)
C   --   LINETHICKNESS             Specify thickness for for mesh lines.
C   --   SPHERE or FSPHERE {ON/OFF}
C   --                             Display elements as spheres, taking
C   --                             their radius from the first
C   --                             element attribute.
C   --
C   --Display control
C   --   DEADNODE {ON/OFF}         Enable/Disable dead node display
C   --   color                     Set number of colors
C   --   spectrum                  Set number of spectrum colors
C   --
C   --Display
C   --   reset                     Reset default conditions
C   --   postmesh                  Initialize after a mesh plot
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   PLOT                      Exit to plot the plot set
C   --   HARDCOPY                  Exit to plot the plot set on hardcopy device
C   --   mesh                      Exit to plot the mesh
C   --
C   --Information
C   --   show     {option}         Display program information (based on type)
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A      - IN  - the dynamic memory base array
C   --   CURPRO - IN  - the current program name
C   --   MSHTYP - IN  - true iff the current program is a mesh program
C   --   INLINE - I/O - the parsed input lines for the log file
C   --   INVERB - I/O - the command verb
C   --   IFLD,  - I/O - the free-field reader index and fields
C   --   INTYP, - I/O - the free-field reader index and fields
C   --   CFIELD,- I/O - the free-field reader index and fields
C   --   IFIELD,- I/O - the free-field reader index and fields
C   --   RFIELD - I/O - the free-field reader index and fields
C   --   NAMECO - IN  - the coordinate names
C   --   NAMES  - IN  - the variable names
C   --   IDELB  - IN  - the element block ID array
C   --   NEWELB - I/O - the new element blocks flag:
C   --                  0 = no new element blocks (set elsewhere)
C   --                  1 = new selected element blocks
C   --                  2 = new displayed element blocks
C   --                      (implies new selected blocks)
C   --   IELBST - I/O - the element block status:
C   --                  -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NEWMOD - OUT - the mode status of each view:
C   --                  -1 = unchanged
C   --                   0 = changed to default
C   --                   n = changed to be like view n
C   --   IDNPS  - IN  - the node set ID for each set
C   --   ISSNPS - I/O - the indices of the selected node sets
C   --   IDESS  - IN  - the side set ID for each set
C   --   ISSESS - I/O - the indices of the selected side sets
C   --   BLKCOL - I/O - the user selected colors of the element blocks.
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
C   --   Uses NELBLK, NVARNP, NVAREL of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Sets DEFPRO, DFAC, DDFAC of /DEFORM/
C   --   Uses DEFOK, DEFFAC of /DEFORM/
C   --   Sets and uses MSHDEF, MSHNUM, MSHLIN, MLNTYP, IHIDOP, NALVAR, DEADNP
C   --      of /MSHOPT/
C   --   Sets and uses MULTIM, XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM
C   --      of /VIEWS/
C   --   Uses ALMESH of /MSHLIM/
C   --   Sets ZMMESH, RDMESH, TICMSH, MSCTYP, SQMESH of /MSHLIM/
C   --   Sets and uses NEWROT, ROTMAT, ROTCEN, EYE of /ROTOPT/
C   --   Sets NEWCUT, ISCUT, CUTPT, CUTNRM of /CUTOPT/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'views.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'cutopt.blk'
      include 'selne.blk'
      include 'sphele.blk'
      include 'light.blk'
      include 'debug.blk'

      DIMENSION A(*)
      CHARACTER*(*) CURPRO
      LOGICAL MSHTYP
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAMES(*)
      INTEGER IDELB(NELBLK)
      INTEGER NEWELB
      INTEGER IELBST(NELBLK)
      INTEGER NEWMOD(4)
      INTEGER IDNPS(*)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER IDESS(*)
      INTEGER BLKCOL(0:NELBLK)
      REAL    SHDCOL(7,NELBLK)
      INTEGER ISHDCL(3,NELBLK)

      LOGICAL FFMATC, FFEXST, MATSTR
      INTEGER NDEFVW, IXVW
      CHARACTER*(MXNAME) WORD
      CHARACTER CDUM
      LOGICAL LDUM1, LDUM2
      LOGICAL HELP
      LOGICAL VWCMD
      LOGICAL ISON
      INTEGER LTYP(-1:1)
      LOGICAL MESHOK

      CHARACTER*(MXNAME) VERB, VERB2
      CHARACTER*(MXNAME) SAVERB
C      --VERB - the command verb used in the SHOW
C      --VERB2 - the secondary SHOW command verb

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) OLDPRO
      SAVE OLDPRO
C      --OLDPRO - the name of the program that set the defaults for
C      --   the current parameters
      CHARACTER*8 MSCSV
      SAVE MSCSV
C      --MSCSV - the old MSCTYP if set to 'SELECTED'

      LOGICAL NEWZM, SETTIC
      SAVE NEWZM, SETTIC
C      --NEWZM - true iff a new zoom window has been specified
C      --SETTIC - true iff the axis tick interval has been specified by the user
      REAL DFRMAT(3,3), DFRCEN(3)
      SAVE DFRMAT, DFRCEN
C      --DFRMAT - the default rotation matrix, saved on first entry
C      --DFRCEN - the default rotation center, saved on first entry

      CHARACTER*(MXSTLN) CMDTBL(41)
      INTEGER IVWTBL, NVWTBL
      SAVE CMDTBL, IVWTBL, NVWTBL
C      --CMDTBL - the command table
C      --IVWTBL - the VIEW command table starting index
C      --NVWTBL - the number of entries in the VIEW command table

      DATA FIRST / .TRUE. /

      DATA IVWTBL, NVWTBL / 1, 0 /

C   --Command table follows.  The "VIEW" entry starts the VIEW
C   --command table, and all VIEW commands (and only these commands)
C   --must follow this entry.  The table is ended by a blank entry.
C   --Remember to change the dimensioned size when changing the table.
C   --Lowercase entries are for obsolete commands; they are converted to
C   --uppercase before use.
      DATA (CMDTBL(I),I=1,10) /
     1  'VISIBLE                         ',
     *  'BLOCKS                          ',
     *  'MATERIAL                        ',
     *  'DEATH                           ',
     2  'XSYM                            ',
     *  'YSYM                            ',
     *  'XVIEW                           ',
     *  'YVIEW                           ',
     *  'MULTTIME                        ',
     3  'MAGNIFY                         '/
      DATA (CMDTBL(I),I=11,20) /
     *  'HIDDEN                          ',
     4  'ZOOM                            ',
     *  'scale                           ',
     *  'SQUARE                          ',
     *  'TICK                            ',
     5  'ROTATE                          ',
     *  'EYE                             ',
     *  'CENTER                          ',
     *  'CUT                             ',
     *  'LIGHTS                          '/
      DATA (CMDTBL(I),I=21,30) /
     6  'DEADNODE                        ',
     *  'SHDCOLOR                        ',
     *  'AMBIENT                         ',
     7  'PLOT                            ',
     *  'HARDCOPY                        ',
     *  'BLKCOLOR                        ',
     *  'MATCOLOR                        ',
     *  'LINETHICKNESS                   ',
     8  'SPHERE                          ',
     *  'FSPHERE                         '/
      DATA (CMDTBL(I),I=31,41) /
     *  'VIEW                            ',
     *  'EMPTY                           ',
     *  'DEFORM                          ',
     9  'undeform                        ',
     *  'NUMBER                          ',
     *  'MLINES                          ',
     *  'overlay                         ',
     *  'BOUNDARY                        ',
     1  'NSETS                           ',
     *  'SSETS                           ',
     *  '                                ' /

C   --Find the command verb, which may be a variable name

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .NE. ' ') THEN
         VWCMD = (LOCSTR (VERB, NVWTBL, CMDTBL(IVWTBL)) .GT. 0)
      ELSE
         VERB = WORD
         VWCMD = .FALSE.
      END IF

      VERB2 = ' '

C   --Reset NEWMOD

      CALL INIINT (4, -999, NEWMOD)

C *** Initialization ***

      IF ((VERB .EQ. 'postmesh')
     &   .OR. (VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
     &   .OR. (VERB .EQ. 'initres') .OR. (VERB .EQ. 'reset')) THEN
         INVERB = VERB

C      --Initialize parameters first time through, then reset

         IF (FIRST) THEN

C         --Change the command table to upper case

            L = LOCSTR (' ', 999, CMDTBL) - 1
            DO 100 I = 1, L
               CALL EXUPCS (CMDTBL(I))
  100       CONTINUE

C         --Find the starting entry in the VIEW table

            IVWTBL = LOCSTR ('VIEW', 999, CMDTBL)
            NVWTBL = LOCSTR (' ', 999, CMDTBL(IVWTBL)) - 1

            DDFAC = DEFFAC
            IF (DDFAC .LT. 0.0) DDFAC = 0.0

C         --Set up the default rotation matrix and rotation center

            IF (IS3DIM) THEN
               CALL INIREA (3*3, 0.0, DFRMAT)
               DO 110 I = 1, 3
                  DFRMAT(I,I) = 1.0
  110          CONTINUE

               DFRCEN(1) = 0.5 * (UNMESH(KRGT)-UNMESH(KLFT))
               DFRCEN(2) = 0.5 * (UNMESH(KTOP)-UNMESH(KBOT))
               DFRCEN(3) = 0.5 * (UNMESH(KFAR)-UNMESH(KNEA))
            END IF

C         --Set values that do not change when reset

            LFTSYM = .TRUE.
            XAXSYM = UNMESH(KLFT)
            BOTSYM = .TRUE.
            YAXSYM = UNMESH(KBOT)

C         --Initialize MSCTYP (so SCALAX can be called)

            IF (.NOT. IS3DIM) THEN
               MSCTYP = 'MESH'
            ELSE
               MSCTYP = 'EACH'
            END IF

            VERB = 'initres'

C     --Initialize color map
            CALL STCLST('DETOUR')
         END IF

C      --Initialize for program change

         IF (VERB .EQ. 'initprog') THEN

C         --Mesh is undeformed (may be changed)

            DEFPRO = .FALSE.
            DFAC = 0.0
            CALL SCALAX

            IF ((.NOT. MSHTYP) .OR. (OLDPRO .NE. CURPRO)) THEN

C            --Set mode for all views

               CALL INIINT (3, 1, LTYP)
               CALL SETMSH (2, 'UNDEFORM', 'NONE', MSHSEL, LTYP,
     &            0, IDUM, 0, IDUM, 'WIREFRAM', ' ', ISSNPS, ISSESS)

C            --Set axis scale type and set axis scale

               IF (MSCTYP .EQ. 'SELECTED') THEN
                  SETTIC = .FALSE.
                  TICMSH = 0.0
                  NEWZM = .TRUE.
                  IF (.NOT. IS3DIM) THEN
                     MSCTYP = 'MESH'
                  ELSE
                     MSCTYP = 'EACH'
                  END IF

                  CALL SCALAX
               END IF
            END IF

            OLDPRO = CURPRO
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN

C         --Set active element control to all element blocks, no birth/death

            IF (FIRST) NEWELB = 0
            NEWELB = MAX (NEWELB, 1)
            DO 120 IELB = 1, NELBLK
               IF ((.NOT. FIRST) .AND. (IELBST(IELB) .LT. 0)) NEWELB = 2
               IELBST(IELB) = 1
  120       CONTINUE

            NALVAR = 0

C         --Mesh is undeformed (may be changed)

            DEFPRO = .FALSE.
            DFAC = 0.0
            CALL SCALAX

C         --Set rotation matrix and center and cutting plane

            IF (IS3DIM) THEN
               NEWROT = .TRUE.
               CALL CPYREA (3*3, DFRMAT, ROTMAT)
               CALL CPYREA (3, DFRCEN, ROTCEN)
               Z = UNMESH(KFAR) - UNMESH(KNEA)
               CALL UNROT (1, 1, ROTMAT, ROTCEN,
     &            0.0, 0.0, Z, EYE(1), EYE(2), EYE(3))

               ISCUT = .FALSE.
               CALL INIREA (3, 0.0, CUTPT)
               CALL INIREA (3, 0.0, CUTNRM)
            END IF

C         --Set wireframe mode on 1 view

            MULTIM = .FALSE.
            XISSYM = .FALSE.
            YISSYM = .FALSE.

            IF (.NOT. IS3DIM) THEN
               IHIDOP = 0
            ELSE
               if (nshl .gt. 100000) then
                  ihidop = 2
                  CALL PRTERR ('CMDREQ',
     &            'Setting "hidden 2" for faster shell plotting')
                  CALL PRTERR ('CMDREQ',
     &            'Reset via "hidden 3" when zoomed')
               else
                  IHIDOP = 3
               end if
            END IF
            CALL SETMSH (0, 'NONE', CDUM, IDUM, LTYP,
     &         IDUM, IDUM, IDUM, IDUM, 'NONE', CDUM, ISSNPS, ISSESS)
            CALL INIINT (3, 1, LTYP)
            CALL SETMSH (2, 'UNDEFORM', 'NONE', MSHSEL, LTYP,
     &         0, IDUM, 0, IDUM, 'WIREFRAM', ' ', ISSNPS, ISSESS)
            CALL SCOLOR (.TRUE., CDUM, IDUM, IDUM, RDUM, IDUM,
     *        CDUM, SHDCOL, ISHDCL, IDELB)

C         --Set display options

            DEADNP = .FALSE.

C         --Set axis scale type and set axis scale

            TICMSH = 0.0
            SETTIC = .FALSE.
            NEWZM = .TRUE.
            SQMESH = .TRUE.
            IF (.NOT. IS3DIM) THEN
               MSCTYP = 'MESH'
            ELSE
               MSCTYP = 'EACH'
            END IF

            CALL SCALAX

            OLDPRO = CURPRO

            FIRST = .FALSE.

C              reset line thicknesses
            CALL LINTHK (INLINE, IFLD, INTYP, IFIELD, RFIELD,
     &         CFIELD, .TRUE.)

C              reset whether to plot elements as spheres

            SPHPLT = 0
         END IF

C      --Initialize after mesh plot

         IF (VERB .EQ. 'postmesh') THEN

C         --Reset mesh scaling factor if selected

            IF (MSCTYP .EQ. 'SELECTED') MSCTYP = MSCSV
         END IF

C      --Initialize for new plot set

         IF (((VERB .EQ. 'postplot') .AND. MSHTYP)
     &      .OR. (VERB .EQ. 'postmesh')) THEN
            CALL QNPICK ('QUERY', ISON, LDUM2,
     &         A, IDUM, IDUM, IDUM, IDUM, IDUM)
            IF (ISON) THEN
               NEWZM = .FALSE.
               SETTIC = .FALSE.
            END IF
         END IF

         OLDPRO = CURPRO
         VERB = ' '

C *** Display Mode control ***

      ELSE IF (VWCMD) THEN
         SAVERB = INVERB
         INVERB = ' '
         IF ((VERB .EQ. 'DEFORM') .or. (verb .eq. 'UNDEFORM')) THEN
            IF (.NOT. DEFPRO) THEN
               CALL PRTERR ('CMDERR',
     &            'Command not allowed in this subprogram')
               GOTO 140
            END IF
         END IF

C      --Set up the view number
         CALL CMDVWC (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      IVIEW, JVIEW, *140)

         WORD = VERB
         CALL ABRSTR (VERB, WORD, CMDTBL(IVWTBL))
         IF (VERB .EQ. ' ') THEN
            INVERB = SAVERB
            GOTO 140
         END IF

         CALL CMDMSH (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      IVIEW, JVIEW, NEWMOD,
     &      IDNPS, ISSNPS, IDESS, ISSESS, *140)

         INVERB = 'newmod'

C *** Active element control ***

      ELSE IF ((VERB .EQ. 'VISIBLE')
     &   .OR. (VERB .EQ. 'BLOCKS') .OR. (VERB .EQ. 'MATERIAL')
     &   .OR. (VERB .EQ. 'DEATH')) THEN
         INVERB = ' '

         IF ((VERB .EQ. 'VISIBLE')
     &      .OR. (VERB .EQ. 'BLOCKS') .OR. (VERB .EQ. 'MATERIAL')) THEN
            if ((verb .eq. 'BLOCKS') .or. (verb .eq. 'MATERIAL')) then
               if (ffmatc (ifld, intyp, cfield, 'DISPLAY', 1)) then
                  call prterr ('CMDREQ',
     &               'Please use the VISIBLE command')
                  verb = 'VISIBLE'
               end if
            end if

            CALL CMDELB (VERB, INLINE(1), IFLD, INTYP, CFIELD,
     &                   IFIELD, IDELB, IELBST, NEWELB, *140)

         ELSE IF (VERB .EQ. 'DEATH') THEN
            IF (.NOT. DEFPRO) THEN
               CALL PRTERR ('CMDERR',
     &            'Command not allowed in this subprogram')
               GOTO 140
            END IF
            CALL DBVIX_BL ('E', 1, IEV)
            CALL CMDDEA (VERB, INLINE(1), IFLD, INTYP, CFIELD,
     &         RFIELD, NAMES(IEV), NALVAR, ALIVAL, *140)

            IF (NALVAR .GT. 0) VERB2 = 'DEADNODE'
         END IF

         INVERB = 'newele'

C *** Multiple views control ***

      ELSE IF ((VERB .EQ. 'XSYM') .OR. (VERB .EQ. 'XVIEW')
     &   .OR. (VERB .EQ. 'YSYM') .OR. (VERB .EQ. 'YVIEW')) THEN
         INVERB = ' '

         CALL CMDMVW (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      UNMESH, NEWMOD, ISSNPS, ISSESS, *140)
         VERB2 = 'VIEW'

         CALL SCALAX

         INVERB = 'newmod'

      ELSE IF (VERB .EQ. 'MULTTIME') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         IF (.NOT. DEFPRO) THEN
            CALL PRTERR ('CMDERR',
     &         'Command not allowed in this subprogram')
            GOTO 140
         END IF

         CALL FFONOF (IFLD, INTYP, CFIELD, MULTIM, *140)
         CALL FFADDO (MULTIM, INLINE(1))

C *** Mesh control ***

      ELSE IF (VERB .EQ. 'MAGNIFY') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         IF (.NOT. DEFPRO) THEN
            CALL PRTERR ('CMDERR',
     &         'Command not allowed in this subprogram')
            GOTO 140
         END IF
         IF (.NOT. DEFOK) GOTO 140
         IF (DEFFAC .LT. 0.0) THEN
            CALL PRTERR ('CMDERR',
     &         'Magnification is not initialized - Call DETOUR')
            GOTO 140
         END IF

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'magnification factor', DEFFAC, DFAC, *140)
         CALL FFADDR (DFAC, INLINE(1))
         IF (DFAC .NE. 0.0) DDFAC = DFAC

         CALL SCALAX

      ELSE IF (VERB .EQ. 'HIDDEN') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         IF (.NOT. IS3DIM) THEN
            CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
            GOTO 140
         END IF

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'hidden line removal level', IHIDOP, IHID, *140)
         IHIDOP = MAX (0, IHID)
         CALL FFADDI (IHIDOP, INLINE(1))

      ELSE IF ((VERB .EQ. 'ZOOM') .OR. (VERB .EQ. 'TICK')
     &   .or. (verb .eq. 'SCALE') .OR. (VERB .EQ. 'SQUARE')) THEN
         INVERB = ' '

         CALL MDFIND ('MAPND', KMAPND, LDUM)

         CALL CMDZM (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NEWZM, SETTIC, A(KMAPND), A, *140)

         CALL SCALAX

      ELSE IF ((VERB .EQ. 'ROTATE') .OR. (VERB .EQ. 'EYE')
     &   .OR. (VERB .EQ. 'CENTER')) THEN
         INVERB = ' '

         CALL CMDROT (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      DFRMAT, DFRCEN,
     &      NEWZM, A, *140)

         CALL SCALAX

      ELSE IF (VERB .EQ. 'LIGHTS') THEN
        INVERB = ' '
        CALL FFADDC (VERB, INLINE(1))
        if (FFMATC (IFLD, INTYP, CFIELD, 'ADD', 3)) THEN
          NLIT = NLIT + 1
          if (NLIT .GT. MAXLIT) THEN
            call prterr('WARNING',
     *        'Too many lights defined, overwriting last defined light')
            NLIT = MAXLIT
          end if
          CALL FFADDC ('ADD', INLINE(1))
        else if (FFMATC (IFLD, INTYP, CFIELD, 'DELETE', 3)) THEN
          CALL FFADDC ('DELETE', INLINE(1))
          IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &        'light number', 0, LITN, *125)
            IF (litn .gt. 0 .and. litn .le. nlit) then
              do 122 i = litn+1, nlit
                do 121 j = 1, 4
                  lite(j,i-1) = lite(j,i)
 121            continue
 122          continue
              nlit = nlit - 1
            ELSE
              call prterr ('WARNING', 'Undefined Light')
            END IF
          ELSE
            NLIT = 0
          END IF
          go to 125
        else
          nlit = 1
        end if

        CALL FFREAL (IFLD, INTYP, RFIELD,
     *    'X Component', 1.0, LITE(1,NLIT), *140)
        CALL FFADDR (LITE(1,NLIT), INLINE(1))
        CALL FFREAL (IFLD, INTYP, RFIELD,
     *    'Y Component', 1.0, LITE(2,NLIT), *140)
        CALL FFADDR (LITE(2,NLIT), INLINE(1))
        CALL FFREAL (IFLD, INTYP, RFIELD,
     *    'Z Component', 1.0, LITE(3,NLIT), *140)
        CALL FFADDR (LITE(3,NLIT), INLINE(1))
        CALL FFREAL (IFLD, INTYP, RFIELD,
     *    'Light Intensity', 1.0, LITE(4,NLIT), *140)
        CALL FFADDR (LITE(4,NLIT), INLINE(1))
 125    continue
C ... Normalize light vectors.
        do 127 i = 1, nlit
          vmag = sqrt(lite(1,i)**2 + lite(2,i)**2 + lite(3,i)**2)
          lite(5,i) = lite(1,i) / vmag
          lite(6,i) = lite(2,i) / vmag
          lite(7,i) = lite(3,i) / vmag
          lite(8,i) = lite(4,i)
 127    continue

      ELSE IF (VERB .EQ. 'AMBIENT') THEN
        INVERB = ' '
        CALL FFADDC (VERB, INLINE(1))
        CALL FFREAL (IFLD, INTYP, RFIELD,
     *    'Ambient Intensity', 0.2, AMBIENT, *140)
        CALL FFADDR (AMBIENT, INLINE(1))

      ELSE IF ((VERB .EQ. 'CUT') .OR. (VERB .EQ. 'CUT3')) THEN
         INVERB = ' '

         CALL CMDCUT (VERB, INLINE(1), IFLD, INTYP, CFIELD,
     &                RFIELD, A, *140)

      ELSE IF ((VERB .EQ. 'WHAT') .OR. (VERB .EQ. 'WHAT3')
     &   .OR. (VERB .EQ. 'WHERE')) THEN
         INVERB = ' '

         CALL QNPICK ('ORIGINAL', LDUM1, LDUM2,
     &      A, KXN, KYN, KZN, IDUM, IDUM)
         CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &      A, KXNN, KYNN, KZNN, KHIDEN, KNPSUR)
         CALL MDFIND ('XE', KXE, KLEN)
         CALL MDFIND ('YE', KYE, KLEN)
         IF (IS3DIM) THEN
           CALL MDFIND ('ZE', KZE, KLEN)
         ELSE
           KZE = 1
         END IF
         CALL MDFIND ('MAPEL', KMAPEL, LDUM)
         CALL MDFIND ('MAPND', KMAPND, LDUM)

         CALL CMDWHE (VERB, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      A(KXN), A(KYN), A(KZN),
     &      A(KXNN), A(KYNN), A(KZNN), A(KHIDEN),
     *      A(KXE), A(KYE), A(KZE), A(KMAPEL), A(KMAPND), *140)

      ELSE IF (VERB .EQ. 'LINETHICKNESS') THEN

         CALL MSCHK (.FALSE., MESHOK)
         IF (MESHOK) THEN
            CALL FFADDC (VERB, INLINE(1))
            CALL LINTHK (INLINE, IFLD, INTYP, IFIELD, RFIELD,
     &         CFIELD, .FALSE.)
            IF (INLINE(1) .EQ. ' ') VERB = ' '
         ELSE
            WRITE (*, 10010)
            VERB = ' '
         ENDIF
         INVERB = ' '

      ELSE IF (VERB .EQ. 'SPHERE'  .OR.  VERB .EQ. 'FSPHERE') THEN
         CALL FFADDC (VERB, INLINE(1))
         ISPHER = 0
         DEFRAD = 1.0
         IF (FFEXST (IFLD, INTYP)) THEN
C           Check that next field has characters in it.
            IF (INTYP(IFLD) .GE. 0) THEN
C              Check for 'ON' or 'OFF'
               IF (CFIELD(IFLD) .EQ. 'ON') THEN
                  IF (VERB .EQ. 'SPHERE') THEN
                     SPHPLT = 1
                  ELSE
                     SPHPLT = -1
                  ENDIF
                  CALL FFADDC ('ON', INLINE(1))
                  ISPHER = 1
                  IFLD = IFLD + 1
               ELSE IF (CFIELD(IFLD) .EQ. 'OFF') THEN
                  SPHPLT = 0
                  CALL FFADDC ('OFF', INLINE(1))
                  ISPHER = 1
                  IFLD = IFLD + 1
               ELSE
                  WRITE (*, 10000)
10000              FORMAT (1X, 'EXPECTED "ON" OR "OFF" AS',
     &               ' PARAMETER TO SPHERE/FSPHERE COMMAND'/)
               ENDIF
            ELSE
               WRITE (*, 10000)
            ENDIF
         ENDIF
         IF (ISPHER .EQ. 0) THEN
            IF (SPHPLT .EQ. 0) THEN
               IF (VERB .EQ. 'SPHERE') THEN
                  SPHPLT = 1
               ELSE
                  SPHPLT = -1
               ENDIF
            ELSE IF (SPHPLT .GE. 1) THEN
               IF (VERB .EQ. 'SPHERE') THEN
                  SPHPLT = 0
               ELSE
                  SPHPLT = -1
               ENDIF
            ELSE
               IF (VERB .EQ. 'SPHERE') THEN
                  SPHPLT = 1
               ELSE
                  SPHPLT = 0
               ENDIF
            ENDIF
         ENDIF
C ... Optional "RADIUS default_radius" option
         IF (FFEXST (IFLD, INTYP)) THEN
C          Check that next field has characters in it.
           IF (INTYP(IFLD) .GE. 0) THEN
             IF (MATSTR (CFIELD(IFLD), 'RADIUS', 1)) THEN
               IFLD = IFLD + 1
               CALL FFADDC ('RADIUS', INLINE)
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'default sphere radius multiplier', 1.0, DEFRAD, *140)
             END IF
           END IF
         END IF
C ... Optional integer specifying number of sphere segments
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &     'segments used to draw spheres', 12, ISMULT, *140)
         ISMULT = MIN(ISMULT, NPTSPX)
         SPHPLT = SPHPLT * ISMULT

         INVERB = ' '

C *** Display control ***

      ELSE IF (VERB .EQ. 'DEADNODE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         IF (.NOT. DEFPRO) THEN
            CALL PRTERR ('CMDERR',
     &         'Command not allowed in this subprogram')
            GOTO 140
         END IF

         CALL FFONOF (IFLD, INTYP, CFIELD, DEADNP, *140)
         CALL FFADDO (DEADNP, INLINE(1))

         VERB2 = 'DEATH'

      ELSE IF ((VERB .EQ. 'PLOT') .OR. (VERB .EQ. 'HARDCOPY') .OR.
     &   (VERB .EQ. 'mesh')) THEN
         IF (VERB .NE. 'mesh') CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (NDEFVW (.FALSE.) .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'No views are defined')
            VERB = 'PLOT'
            GOTO 140
         END IF

C      --Check if nodes/element selected for mesh plot

         IF (VERB .EQ. 'mesh') THEN

C         --Save mesh scaling type
            MSCSV = MSCTYP

            IF (SELOK) THEN

C            --Set MSCTYP to 'SELECTED' if nodes/element selected for mesh plot
               MSCTYP = 'SELECTED'

C            --Set selected numbering if nodes/element selected for mesh plot
               DO 130 IVW = 1, NDEFVW (.FALSE.)
                  IVIEW = IXVW (.FALSE., IVW)
                  IF (MSHNUM(IVIEW) .EQ. 'NONE')
     &               MSHNUM(IVIEW) = 'SELECTED'
  130          CONTINUE
            END IF
         END IF

C      --PLOT and HARDCOPY are to be passed as lower-case commands
         CALL LOWSTR (INVERB, VERB)
         VERB = ' '

      ELSE IF (VERB .EQ. 'BLKCOLOR'  .OR.  VERB .EQ. 'MATCOLOR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL MSCHK (.FALSE., MESHOK)
         VERB = 'BLKCOLOR'
         IF (MESHOK) THEN
            CALL BCOLOR (.FALSE., INLINE, IFLD, INTYP,
     &         IFIELD, CFIELD, BLKCOL)
         ELSE
            WRITE (*, 10010)
10010        FORMAT (1X, 'Command not valid unless valid mesh is',
     &         ' defined')
            VERB = ' '
         ENDIF
         INVERB = ' '

      ELSE IF (VERB .EQ. 'SHDCOLOR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL MSCHK (.FALSE., MESHOK)
         VERB = 'SHDCOLOR'
         IF (MESHOK) THEN
            CALL SCOLOR (.FALSE., INLINE, IFLD, INTYP,
     &         RFIELD, IFIELD, CFIELD, SHDCOL, ISHDCL, IDELB)
         ELSE
            WRITE (*, 10010)
            VERB = ' '
         ENDIF
         INVERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL ABRSTR (VERB, WORD, CMDTBL)
         IF (VERB .NE. ' ') THEN
            IF (.NOT. MSHTYP) THEN
               IF ((VERB .NE. 'PLOT') .AND. (VERB .NE. 'HARDCOPY')) THEN
                  CALL MSSHOW (VERB, NAMECO, NAMES, IDELB, IELBST,
     &               IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL)
                  INVERB = ' '
               END IF
            ELSE
               CALL MSSHOW (VERB, NAMECO, NAMES, IDELB, IELBST,
     &            IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL)
               IF ((VERB .NE. 'PLOT') .AND. (VERB .NE. 'HARDCOPY')
     &            .AND. (VERB .NE. 'VIEW')) THEN
                  INVERB = ' '
               ELSE
                  CFIELD(IFLD-1) = VERB
               END IF
            END IF
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON) THEN
            CALL SHOCMD ('Mesh Commands', CMDTBL)
         END IF
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 150

  140 CONTINUE
      INLINE(1) = ' '

  150 CONTINUE

      IF (VERB .NE. ' ') THEN
         CALL MSSHOW (VERB, NAMECO, NAMES, IDELB, IELBST,
     &      IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL)
      ELSE
         VERB2 = ' '
      END IF

      IF (VERB2 .NE. ' ') THEN
         CALL MSSHOW (VERB2, NAMECO, NAMES, IDELB, IELBST,
     &      IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL)
      END IF

      RETURN
      END
