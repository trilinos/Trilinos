C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DTCOMD (A, INLINE, INVERB, IFLD, INTYP, CFIELD,
     &                   IFIELD, RFIELD, NEWMOD, NAMES, IELBST,
     &                   ISSNPS, ISSESS, LIDSP, MAPEL, MAPND, NAMLEN)
C=======================================================================

C   --*** DTCOMD *** (DETOUR) Process DETOUR commands
C   --   Modified by John Glick - 11/1/88
C   --   Written by Amy Gilkey - revised 05/16/88
C   --   Dennis Flanagan, 11/18/82
C   --
C   --DTCOMD interprets DETOUR commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Display mode control
C   --   newmod
C   --   WIREFRAM                  Set wireframe mesh mode
C   --   SOLID                     Set solid mesh mode
C   --   NSETS    {id1,...}        Set sets mode and select node sets
C   --   SSETS    {id1,...}        Set sets mode and select side sets
C   --   CONTOUR  {variable}       Set line contour mode
C   --   PAINT    {variable}       Set paint contour mode
C   --   VECTOR   {2 or 3 variables} Set vector mode
C   --   SYMBOL   {type,variable}  Set element symbol mode
C   --   GAUSS    {4 variables}    Set element symbol gauss mode
C   --   variable                  Select contour variable, and set contour
C   --                             range; in vector mode, set the vector
C   --                             variables; in symbol mode, set the symbol
C   --                             variable
C   --
C   --Active element control
C   --   newele
C   --
C   --Contour control
C   --   NCNTRS   {ncntrs}         Set number of contours, round DELC
C   --   CRANGE   {cmin,cmax}      Set minimum and maximum contour value
C   --   CMIN     {cmin,cmax}      Set minimum contour value
C   --   CMAX     {cmax}           Set maximum contour value
C   --   CSHIFT   {cval}           Shift contour limits to fall on CVAL
C   --   DELCNTR  {delc,cmin}      Set contour interval, adjust limits
C   --   CINTV    {c1,c2,...}      Set specific contour intervals
C   --
C   --Display control
C   --   COPEN    {ON/OFF,ON/OFF}  Set contour minimum / maximum open or closed
C   --   CLABEL   {OFF/ON/labinc}  Label every LABINC interior mesh line
C   --                             with contour letter
C   --   CSYMBOLS {maxmin,maxmax}  Specify maximum number of min/max symbols
C   --                             to plot on contour plots.
C   --   VSCALE   {vecscl}         Set vector length / symbol scale factor
C   --   color                     Set number of colors
C   --   spectrum                  Set number of spectrum colors
C   --   DISPVAR [ADD] {var1,...}  Selects history variables, global
C   --                             variables, and/or TIME whose values
C   --                             will be displayed on the plot legend.
C   --
C   --Display
C   --   reset                     Reset default conditions
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   plot                      Exit to plot the plot set
C   --   hardcopy                  Exit to plot the plot set on hardcopy device
C   --
C   --Information
C   --   show     {option}         Display program information (based on type)
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A         - IN  - the dynamic memory base array
C   --   INLINE    - I/O - the parsed input lines for the log file
C   --   INVERB    - I/O - the command verb
C   --   IFLD,
C   --   INTYP,
C   --   CFIELD,
C   --   IFIELD,
C   --   RFIELD    - I/O - the free-field reader index and fields
C   --   NEWMOD    - IN  - the mode status of each view:
C   --                     -1 = unchanged
C   --                      0 = changed to default
C   --                      n = changed to be like view n
C   --   NAMES      - IN  - the variable names
C   --   IELBST     - IN  - the element block status:
C   --                      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   ISSNPS     - I/O - the indices of the selected node sets
C   --   ISSESS     - I/O - the indices of the selected side sets
C   --   LIDSP(0:*) - I/O - the indices of the selected variables
C   --                      whose values will be displayed on the plot legend.
C   --                      LIDSP(0) = the number of variables in the list.
C   --                      LIDSP(i) identifies the ith variable in the list.
C   --                      LIDSP(i) > 0, LIDSP(i) is the id of a history var.
C   --                      LIDSP(i) < 0, -LIDSP(i) is the id of a global var.
C   --                      LIDSP(i) = 0, TIME is displayed on the plot legend.
C   --
C   --Common Variables:
C   --   Uses NELBLK, NVARNP, NVAREL of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/
C   --   Uses NPTIMS of /TIMES/
C   --   Sets NTIMIN, HISTOK of /TIMES/
C   --   Sets VECSCL of /ETCOPT/
C   --   Sets DEFPRO, DFAC of /DEFORM/
C   --   Uses DDFAC of /DEFORM/
C   --   Sets and uses MSHDEF, MSHNUM, MSHLIN, MLNTYP of /MSHOPT/
C   --   Sets and uses MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR of /DETOPT/
C   --   Sets and uses CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC, CINTV,
C   --      NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX of /CNTR/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      EXTERNAL BLKDAT

      include 'debug.blk'
C      common /debugc/ cdebug
C      common /debugn/ idebug
C      character*8 cdebug

      include 'params.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'times.blk'
      include 'etcopt.blk'
      include 'deform.blk'
      include 'mshopt.blk'
      include 'detopt.blk'
      include 'cntr.blk'

C  Flag for exact contour values for each plot
      include 'icexct.blk'
      include 'axsplt.blk'

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      INTEGER NEWMOD(4)
      CHARACTER*(*) NAMES(*)
      INTEGER IELBST(NELBLK)
      INTEGER ISSNPS(NUMNPS,4)
      INTEGER ISSESS(NUMESS,4)
      INTEGER LIDSP(0:*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL FFNUMB, MATSTR
      INTEGER NUMMOD
      CHARACTER*(MXNAME) VERB, VERB2
C      --VERB - the command verb used in the SHOW
C      --VERB2 - the secondary SHOW command verb

      CHARACTER*(MXNAME) WORD
      LOGICAL HELP
      LOGICAL VWCMD
      LOGICAL ISON
      LOGICAL INL, INE
      CHARACTER TYP
      INTEGER LTYP(-1:1)

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      INTEGER NVOLD(4), NCOLD
      LOGICAL FIXCON
      SAVE NVOLD, NCOLD, FIXCON
C      --NVOLD - the last variables
C      --NCOLD - the last contour variable (for range)
C      --FIXCON - true iff contour parameters are fixed
      REAL FMINC, FMAXC
      SAVE FMINC, FMAXC
C      --FMINC, FMAXC - the minimum and maximum variable values,
C      --   saved to recalculate contour range

      CHARACTER*(MXSTLN) CMDTBL(29)
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
      DATA CMDTBL /
     1   'NCNTRS  ', 'CRANGE  ', 'CMIN    ', 'CMAX    ',
     2   'CSHIFT  ', 'DELCNTR ', 'CINTV   ', 'CLABEL  ',
     3   'lines   ', 'COPEN   ', 'DISPVAR ', 'opencntr',
     4   'CSYMBOLS', 'minmax  ', 'VSCALE  ', 'vecscl  ', 'CGLOBAL ',
     5   'VIEW    ', 'WIREFRAM', 'SOLID   ',
     6   'CONTOUR ', 'PAINT   ', 'EPAINT  ',
     7   'VECTOR  ', 'SYMBOL  ', 'GAUSS   ',
     8   'deform  ', 'undeform',
     9   '        ' /

C   --Find the command verb, which may be a variable name

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)

      IF (VERB .NE. ' ') THEN
         VWCMD = (LOCSTR (VERB, NVWTBL, CMDTBL(IVWTBL)) .GT. 0)
      ELSE
         CALL DBVIX_BL ('N', 1, INV)
         CALL DBVIX_BL ('E', 1, IEV)
         if (inv .gt. 0) then
           inl = LOCSTR (WORD, NVARNP, NAMES(INV)) .GT. 0
         else
           inl = .FALSE.
         end if
         if (iev .gt. 0) then
           ine = LOCSTR (WORD, NVAREL, NAMES(IEV)) .GT. 0
         else
           ine = .FALSE.
         end if
         if (inl .or. ine) then
            VERB = ' '
            VWCMD = .TRUE.
            IFLD = IFLD - 1
         ELSE
            VERB = WORD
            VWCMD = .FALSE.
         END IF
      END IF

      VERB2 = ' '

C -- CHECK FOR AXIS ONLY PLOTTING
      IF(VERB .EQ. 'plot') THEN
          IF( MATSTR(CFIELD(IFLD),'AXIS', 2)) THEN
              AXONLY = .TRUE.
          END IF
      END IF

C *** Initialization ***

      IF ((VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
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

            VERB = 'initres'
         END IF

C      --Initialize time step selection default

         IF ((VERB .EQ. 'initprog') .OR. (VERB .EQ. 'initres')) THEN

C         --Set default time step selection to 10 times with delta interval

            HISTOK = .FALSE.
            NTIMIN = 10
         END IF

C      --Initialize for program change

         IF (VERB .EQ. 'initprog') THEN

C         --Set default magnification factor

            DEFPRO = .TRUE.
            IF (DFAC .EQ. 0.0) DFAC = DDFAC
            CALL SCALAX

C         --Set deformed wireframe mode on all views

            CALL INIINT (3, 1, LTYP)
            CALL SETMSH (0, 'DEFORM', 'NONE', MSHSEL, LTYP,
     &         0, IDUM, 0, IDUM, 'WIREFRAM', ' ', ISSNPS, ISSESS)
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN

            IF (VERB .EQ. 'reset') THEN
C         --Set list of display variables
               CALL DISPV (.TRUE., INLINE, IFLD, INTYP,
     &            CFIELD, NAMES, LIDSP, NAMLEN)
            ENDIF

C         --Set default magnification factor

            DEFPRO = .TRUE.
            IF (DFAC .EQ. 0.0) DFAC = DDFAC
            CALL SCALAX

C         --Set wireframe on all defined views

            CALL INIINT (3, 1, LTYP)
            CALL SETMSH (0, 'DEFORM', 'NONE', MSHSEL, LTYP,
     &         0, IDUM, 0, IDUM, 'WIREFRAM', ' ', ISSNPS, ISSESS)

            CALL INIINT (4, 0, IDTVAR)
            CALL INIINT (4, 0, NVOLD)

            NCOLD = 0
            FIXCON = .FALSE.
            CINTOK = .FALSE.
            CALL INIREA (256, 0.0, CINTV)
            LINCON = .TRUE.
            NCNTR = 6
            CMIN = 0.0
            CMAX = 0.0
            DELC = 0.0

C         --Set display options

            LABINC = 0
            NOCMIN = .FALSE.
            NOCMAX = .FALSE.
            MAXMIN = 10
            MAXMAX = 10

            VECSCL = 1.0

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

         VERB = ' '

C *** Display mode control ***

      ELSE IF (VWCMD) THEN
         INVERB = ' '

C      --Set up the view number
         CALL CMDVWC (VERB, INLINE,
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      IVIEW, JVIEW, *150)

         IF (VERB .NE. ' ') THEN
            WORD = VERB
            CALL ABRSTR (VERB, WORD, CMDTBL(IVWTBL))
            IF (VERB .EQ. ' ') THEN
               IF ((LOCSTR (WORD, NVARNP, NAMES(INV)) .GT. 0)
     &            .OR. (LOCSTR (WORD, NVAREL, NAMES(IEV)) .GT. 0)) THEN
                  VERB = ' '
                  IFLD = IFLD - 1
               ELSE
                  VERB = 'VIEW'
                  CALL PRTERR ('CMDERR',
     &               'Expected view-dependent command')
                  GOTO 150
               END IF
            END IF
         END IF

         CALL CMDMOD (VERB, VERB2, INLINE(1),
     &      IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      IVIEW, JVIEW, NAMES, NVOLD, NCOLD, FIXCON, VECSCL,
     &      ISSNPS, ISSESS, *150)

      ELSE IF (VERB .EQ. 'newmod') THEN
         INVERB = ' '

         ISON = .FALSE.
         DO 110 I = 1, 4
            IF (NEWMOD(I) .LT. 0) THEN
               CONTINUE
            ELSE IF (NEWMOD(I) .EQ. 0) THEN
               ISON = .TRUE.
               IF ((MSHDEF(I) .EQ. 'NONE')
     &            .OR. (MSHDEF(I) .EQ. 'EMPTY')) THEN
                  CALL SETMOD (I, 'NONE', ' ')
               ELSE
                  CALL SETMOD (I, 'WIREFRAM', ' ')
               END IF
            ELSE
               ISON = .TRUE.
               J = NEWMOD(I)
               IF ((MSHDEF(J) .EQ. 'NONE')
     &            .OR. (MSHDEF(J) .EQ. 'EMPTY')) THEN
                  CALL SETMOD (I, 'NONE', ' ')
               ELSE
                  CALL SETMOD (I, MODDET(J), MODTYP(J))
               END IF
            END IF
  110    CONTINUE

         VERB = ' '
         IF (ISON) VERB = 'PLOT'
         VWCMD = .TRUE.

      ELSE IF (VERB .EQ. 'newele') THEN
         INVERB = ' '

         VERB = ' '
         VWCMD = .TRUE.

C      --Force all element variables to be re-scaled
         DO 120 I = 1, 4
            CALL DBVTYP_BL (IDTVAR(I), TYP, IDUM)
            IF (TYP .EQ. 'E') NVOLD(I) = 0
  120    CONTINUE

C      --Force the contour range to be re-calculated (if element variable)
         CALL DBVTYP_BL (IDTVAR(1), TYP, IDUM)
         IF ((.NOT. FIXCON) .AND. (TYP .EQ. 'E')) NCOLD = 0

C *** Contour control ***

      ELSE IF ((VERB .EQ. 'NCNTRS') .OR. (VERB .EQ. 'CRANGE')
     &   .OR. (VERB .EQ. 'CMIN') .OR. (VERB .EQ. 'CMAX')
     &   .OR. (VERB .EQ. 'CSHIFT') .OR. (VERB .EQ. 'DELCNTR')
     &   .OR. (VERB .EQ. 'CINTV')) THEN
         INVERB = ' '

         CALL CMDCON (VERB, INLINE(1), IFLD, INTYP,
     &                IFIELD, RFIELD, *150)

C      --Fix the contour type
         FIXCON = .TRUE.
         IEXCON=0

      ELSE IF (VERB .EQ. 'CGLOBAL')THEN
C  CGLOBAL toggles the flag for using the globally defined contours
C  or setting to contour limits to the exact (plus/minus epsilon)
C  values for the plot being drawn
         IF(IEXCON.EQ.0)THEN
          IEXCON=1
          CALL PRTERR('CMDSPEC','Using exact contour values'
     &    //' for each plot.')
         ELSEIF(IEXCON.EQ.1)THEN
          CALL PRTERR('CMDSPEC','Using global contour values'
     &    //' for each plot.')
          IEXCON=0
         ENDIF
         INVERB = ' '

C *** Display control ***

      ELSE IF (VERB .EQ. 'DISPVAR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL DISPV (.FALSE., INLINE, IFLD, INTYP,
     &      CFIELD, NAMES, LIDSP, NAMLEN)
         INVERB = ' '

      ELSE IF ((VERB .EQ. 'CLABEL') .or. (verb .eq. 'LINES')) THEN
         if (verb .eq. 'LINES') then
            call prterr ('CMDREQ', 'Please use the CLABEL command')
            verb = 'CLABEL'
         end if
         CALL FFADDC (VERB, INLINE(1))

         INVERB = ' '

         IF (FFNUMB (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'labeling density', 0, LABINC, *150)
            CALL FFADDI (LABINC, INLINE(1))
         ELSE
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
            IF (MATSTR (WORD, 'ON', 2)) THEN
               CALL FFADDC ('ON', INLINE(1))
               LABINC = 0
            ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
               CALL FFADDC ('OFF', INLINE(1))
               LABINC = -1
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "ON", "OFF" or contour labeling density')
               GOTO 150
            END IF
         END IF

      ELSE IF ((VERB .EQ. 'COPEN') .or. (verb .eq. 'OPENCNTR')) THEN
         if (verb .eq. 'OPENCNTR') then
            call prterr ('CMDREQ', 'Please use the command COPEN')
            verb = 'COPEN'
         end if

         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFONOF (IFLD, INTYP, CFIELD, NOCMIN, *150)
         CALL FFADDO (NOCMIN, INLINE(1))
         CALL FFONOF (IFLD, INTYP, CFIELD, NOCMAX, *150)
         CALL FFADDO (NOCMAX, INLINE(1))

      ELSE IF ((VERB .EQ. 'CSYMBOLS') .or. (verb .eq. 'MINMAX')) THEN
         if (verb .eq. 'MINMAX') then
            call prterr ('CMDREQ', 'Please use the command CSYMBOLS')
            verb = 'CSYMBOLS'
         end if

         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'maximum number of minimums', 10, MM, *150)
         CALL FFADDI (MM, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'maximum number of maximums', MM, MAXMAX, *150)
         CALL FFADDI (MAXMAX, INLINE(1))
         MAXMIN = MM

      ELSE IF ((VERB .EQ. 'VSCALE') .or. (VERB .EQ. 'VECSCL')) THEN
         if (verb .eq. 'VECSCL') then
            call prterr ('CMDREQ', 'Please use the command VSCALE')
            verb = 'VSCALE'
         end if

         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'vector/symbol scale factor', 1.0, VECSCL, *150)
         CALL FFADDR (VECSCL, INLINE(1))

      ELSE IF ((VERB .EQ. 'color') .OR. (VERB .EQ. 'spectrum')) THEN
         VERB = ' '

C      --Reset the number of contours when the colors change

         IF (.NOT. CINTOK) THEN
            NCOL = 0
            IF (VERB .EQ. 'color') THEN
               DO 130 IDEV = 1, 2
                  CALL GRGPARD ('DEVICE', IDEV, ISON, WORD)
                  IF (ISON) THEN
                     CALL GRGPAR ('SPECTRUM', IDEV, IFIELD, WORD)
                     IF (IFIELD(3) .GT. 0) NCOL = -999
                  END IF
  130          CONTINUE
            END IF
            IF (NCOL .EQ. 0) THEN
               DO 140 IDEV = 1, 2
                  CALL GRGPARD ('DEVICE', IDEV, ISON, WORD)
                  IF (ISON) THEN
                     CALL GRGPAR ('SPECTRUM', IDEV, IFIELD, WORD)
                     NCOL = MAX (NCOL, IFIELD(1))
                  END IF
  140          CONTINUE
            END IF
            IF ((NCOL .GT. 0) .AND. (NCNTR .NE. NCOL)) THEN
               NCNTR = NCOL
               IF (NCNTR .GE. 256) THEN
                  CALL PRTERR ('CMDWARN',
     &               'Number of contours reduced to 255')
                  NCNTR = 255
               END IF
               IF (.NOT. FIXCON) THEN
                  VWCMD = .TRUE.
                  NCOLD = 0
               ELSE
                  CALL ADJCON (.FALSE.)
                  VERB = 'NCNTRS'
               END IF
            END IF
         END IF

      ELSE IF ((VERB .EQ. 'plot') .OR. (VERB .EQ. 'hardcopy')) THEN

         IF (NPTIMS .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'No times are selected')
            VERB = ' '
            GOTO 150
         END IF
         CALL CHKVAR (MODDET, MODTYP, -999, IDTVAR, NVOLD,
     &      NNDVAR, NEDVAR, ISON)
         IF (.NOT. ISON) THEN
            VERB = 'PLOT'
            GOTO 150
         END IF

         VERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN
         VERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF ((WORD .NE. 'PLOT') .AND. (WORD .NE. 'HARDCOPY')
     &      .AND. (WORD .NE. 'VIEW')) THEN
            CALL ABRSTR (VERB, WORD, CMDTBL)
         ELSE
            VERB = WORD
         END IF
         IF (VERB .NE. ' ') THEN
            CALL DTSHOW (VERB, NAMES, LIDSP)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN
         VERB = ' '

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON) THEN
            CALL SHOCMD ('DETOUR Commands', CMDTBL)
         END IF

      ELSE
         VERB = ' '
      END IF

      GOTO 160

  150 CONTINUE
      INLINE(1) = ' '

  160 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL DTSHOW (VERB, NAMES, LIDSP)
      ELSE
         VERB2 = ' '
      END IF

      IF (VWCMD) THEN

C      --Scale requested variable

         NNEDVAR = MAX (NNDVAR, NEDVAR)

         DO 170 IVAR = 1, NNEDVAR
            IF (IDTVAR(IVAR) .NE. NVOLD(IVAR)) THEN
               IF (IDTVAR(IVAR) .EQ. 0) THEN
                  FMIN = 0.0
                  FMAX = 0.0
               ELSE
                 CALL SCALER (A, A, 1, NAMES(IDTVAR(IVAR)),
     *             IDTVAR(IVAR),
     &             .TRUE., IELBST, NALVAR, FMIN, FMAX, MAPEL, MAPND)
               END IF

               IF (IVAR .EQ. 1) THEN
                  FMINC = FMIN
                  FMAXC = FMAX
               END IF

               NVOLD(IVAR) = IDTVAR(IVAR)
            END IF
  170    CONTINUE

C      --Compute contour range, if changed

         IF (((NUMMOD (MODDET, ' ', 'CONTOUR', ' ') .GE. 1)
     &      .OR. (NUMMOD (MODDET, ' ', 'ELEMCONT', ' ') .GE. 1))
     &      .AND. (IDTVAR(1) .GT. 0)) THEN
            ISON = (NUMMOD (MODDET, MODTYP, 'CONTOUR', 'LINE') .GE. 1)
     &         .OR. (NUMMOD (MODDET, MODTYP, 'ELEMCONT', 'LINE') .GE. 1)
            IF ((IDTVAR(1) .NE. NCOLD) .OR.
     &         ((.NOT. FIXCON) .AND. (ISON .NEQV. LINCON))) THEN
C            --Recompute and display contour parameters
               CINTOK = .FALSE.
               CALL CONRNG (ISON, FMINC, FMAXC, NCNTR, DELC, CMIN, CMAX)
               NCOLD = IDTVAR(1)
               LINCON = ISON
               IF (VERB2 .EQ. ' ') VERB2 = 'CMIN'
            ELSE IF (ISON .NEQV. LINCON) THEN
C            --Recompute CMAX for current parameters if changed from
C            --line to painted contours
               LINCON = ISON
               CALL ADJCON (.TRUE.)
            END IF
         END IF
      END IF

      IF (VERB2 .NE. ' ') THEN
         CALL DTSHOW (VERB2, NAMES, LIDSP)
      END IF

      RETURN
      END
