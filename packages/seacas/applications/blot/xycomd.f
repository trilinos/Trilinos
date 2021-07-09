C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE XYCOMD (A, CURPRO, INLINE,
     &   INVERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   MESHOK)
C=======================================================================

C   --*** XYCOMD *** (XYPLOT) Process XY-plot commands
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --XYCOMD interprets XY-plot commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Axis size and scaling
C   --   RATIOXY  {ratio}          Set X to Y axis ratio
C   --   XSCALE   {xmin,xmax}      Set X axis scale
C   --   YSCALE   {ymin,ymax}      Set Y axis scale
C   --   XTICK    {xtick}          Set X axis tick interval
C   --   YTICK    {ytick}          Set Y axis tick interval
C   --   XLABEL                    Set X axis label (on following line)
C   --   YLABEL                    Set Y axis label (on following line)
C   --   xaxis                     Reset X axis parameters
C   --   yaxis                     Reset Y axis parameters
C   --   xyaxis                    Reset X and Y axis parameters
C   --
C   --Neutral file options
C   --   ACURVE   {crvnam}         Set neutral file curve name
C   --   NCURVE   {numcrv,inccrv}  Set neutral file curve number and increment
C   --
C   --Display options
C   --   GRID     {ON/OFF}         Draw grid on plot if ON
C   --   LINES    {ON/OFF/VARY}    Plot line, none or vary line
C   --   SYMBOLS  {ON/OFF/#}       Plot varying symbols, none or specific symbol
C   --   CRVNUM   {FIRST/LAST/MIDDLE/OFF} Set curve numbering position
C   --   OVERLAY  {ON/OFF/VARIABLE/TIME} Curves for all variables or all time
C   --                             will be overlaid on one plot
C   --   SAMESCAL {ON/OFF}         Curves will have the same scale if ON
C   --   NORMAL   {ON/OFF}         Curves will be normalized if ON
C   --
C   --Display
C   --   reset                     Reset plot parameters
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   plot                      Exit to plot the plot set
C   --   hardcopy                  Exit to plot the plot set on hardcopy device
C   --   neutral                   Exit to write the plot set to neutral file
C   --
C   --Mesh display commands
C   --   MESH                      Display mesh with numbering, zoom
C   --
C   --Information
C   --   show     {option}         Display plot parameters and information
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INVERB - IN/OUT - the command verb
C   --   CURPRO - IN - the current program name
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   MESHOK - IN - true iff mesh can be displayed
C   --
C   --Common Variables:
C   --   Sets DOGRID, LINTYP, ISYTYP, LABSID, OVERLY, OVERTM,
C   --      IAXSCA of /XYOPT/
C   --   Sets ASPECT, IXSCAL, IYSCAL, XMIN, XMAX, YMIN, YMAX, XTICK, YTICK
C   --      of /XYLIM/
C   --   Sets XLAB, YLAB of /XYLAB/
C   --   Sets NUMCRV, INCCRV and CRVNAM of /NEUTR./

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (NUMSYM = 6, NUMLIN = 6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'xyopt.blk'
      include 'xylim.blk'
      include 'xylab.blk'
      include 'neutr.blk'

      DIMENSION A(*)
      CHARACTER*(*) CURPRO
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      LOGICAL MESHOK

      LOGICAL FFNUMB, MATSTR

      CHARACTER*(MXSTLN) VERB, WORD
      LOGICAL HELP
      LOGICAL ISON

      INTEGER IDUM
      REAL RDUM
      CHARACTER*(MXSTLN) CDUM
      CHARACTER*8 PROMPT
      INTEGER LPROM

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      LOGICAL XSCSET, YSCSET, XTCSET, YTCSET, XLBSET, YLBSET
      SAVE XSCSET, YSCSET, XTCSET, YTCSET, XLBSET, YLBSET
C      --These flags show whether an axis parameter has been set since
C      --   the start of the plot set

      CHARACTER*(MXSTLN) CMDTBL(18)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST /.TRUE./
      DATA IDUM / 1 /

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1  'RATIOXY                         ',
     *  'XSCALE                          ',
     *  'YSCALE                          ',
     *  'XTICK                           ',
     *  'YTICK                           ',
     2  'XLABEL                          ',
     *  'YLABEL                          ',
     3  'ACURVE                          ',
     *  'NCURVE                          ',
     4  'GRID                            ',
     *  'LINES                           ',
     *  'SYMBOLS                         ',
     *  'CRVNUM                          ',
     5  'OVERLAY                         ',
     *  'SAMESCALE                       ',
     *  'NORMAL                          ',
     6  'MESH                            ',
     7  '                                ' /

C   --Get the command verb

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

C *** Initialization ***

      IF ((VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
     &   .OR. (VERB .EQ. 'initres') .OR. (VERB .EQ. 'reset')) THEN
         INVERB = VERB

C      --Initialize parameters first time through, then reset

         IF (FIRST) THEN
            IAXSCA = 'PLOT'
            IXSCAL = IAXSCA
            IYSCAL = IAXSCA

            CRVNAM = ' '
            NUMCRV = 1
            INCCRV = 1

            VERB = 'initres'
         END IF

C      --Initialize for program change

         IF (VERB .EQ. 'initprog') THEN

C         --Leave display options the same

C         --Reset axis scaling and labeling

            ASPECT = 1.0
            IXSCAL = IAXSCA
            IYSCAL = IAXSCA
            XMIN = 0.0
            XMAX = 0.0
            YMIN = 0.0
            YMAX = 0.0
            XTICK = 0.0
            YTICK = 0.0
            XLAB = ' '
            YLAB = ' '
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN
            DOGRID = .FALSE.
            LINTYP = 1
            ISYTYP = 0
            LABSID = 'LAST'
            OVERLY = .FALSE.
            OVERTM = .FALSE.
            IAXSCA = 'PLOT'

            ASPECT = 1.0
            IXSCAL = IAXSCA
            IYSCAL = IAXSCA
            XMIN = 0.0
            XMAX = 0.0
            YMIN = 0.0
            YMAX = 0.0
            XTICK = 0.0
            YTICK = 0.0
            XLAB = ' '
            YLAB = ' '

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

         XSCSET = .FALSE.
         YSCSET = .FALSE.
         XTCSET = .FALSE.
         YTCSET = .FALSE.
         XLBSET = .FALSE.
         YLBSET = .FALSE.

         VERB = ' '

C *** Axis size and scaling ***

      ELSE IF (VERB .EQ. 'RATIOXY') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'X to Y axis length ratio', 1.0, ASPECT, *100)
         CALL FFADDR (ASPECT, INLINE(1))

      ELSE IF (VERB .EQ. 'XSCALE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
            IXSCAL = IAXSCA
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'minimum axis value', XMIN, RMIN, *100)
            CALL FFADDR (ASPECT, INLINE(1))
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'maximum axis value', XMAX, RMAX, *100)
            CALL FFADDR (ASPECT, INLINE(1))
            IF (RMIN .GE. RMAX) THEN
               CALL PRTERR ('CMDERR',
     &            'Axis minimum must be less than axis maximum')
               GOTO 100
            END IF
            XMIN = RMIN
            XMAX = RMAX
            IXSCAL = 'SET'
         END IF

C      --Mark the axis scaling as having been set
         XSCSET = .TRUE.

      ELSE IF (VERB .EQ. 'YSCALE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (.NOT. FFNUMB (IFLD, INTYP)) THEN
            IYSCAL = IAXSCA
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'minimum axis value', YMIN, RMIN, *100)
            CALL FFADDR (ASPECT, INLINE(1))
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'maximum axis value', YMAX, RMAX, *100)
            CALL FFADDR (ASPECT, INLINE(1))
            IF (RMIN .GE. RMAX) THEN
               CALL PRTERR ('CMDERR',
     &            'Axis minimum must be less than axis maximum')
               GOTO 100
            END IF
            YMIN = RMIN
            YMAX = RMAX
            IYSCAL = 'SET'
         END IF

C      --Mark the axis scaling as having been set
         YSCSET = .TRUE.

      ELSE IF (VERB .EQ. 'XTICK') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'tick interval', 0.0, XTICK, *100)
         CALL FFADDR (ASPECT, INLINE(1))

C      --Mark the axis tic mark as having been set
         XTCSET = .TRUE.

      ELSE IF (VERB .EQ. 'YTICK') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &      'tick interval', 0.0, YTICK, *100)
         CALL FFADDR (ASPECT, INLINE(1))

C      --Mark the axis tic mark as having been set
         YTCSET = .TRUE.

      ELSE IF (VERB .EQ. 'XLABEL') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         PROMPT = 'LABEL> '
         LPROM = LENSTR (PROMPT) + 1
         CALL GETINS ('line', IDUM, IDUM, IDUM, CDUM,
     &      IDUM, RDUM, XLAB, IOSTAT,
     &      PROMPT, LPROM, *100)
C        CALL GETINP (0, 0, 'LABEL> ', XLAB, IOSTAT)
         INLINE(2) = XLAB

C      --Mark the axis label as having been set
         XLBSET = .TRUE.

      ELSE IF (VERB .EQ. 'YLABEL') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         PROMPT = 'LABEL> '
         LPROM = LENSTR (PROMPT) + 1
         CALL GETINS ('line', IDUM, IDUM, IDUM, CDUM,
     &      IDUM, RDUM, YLAB, IOSTAT,
     &      PROMPT, LPROM, *100)
C        CALL GETINP (0, 0, 'LABEL> ', YLAB, IOSTAT)
         INLINE(2) = YLAB

C      --Mark the axis label as having been set
         YLBSET = .TRUE.

      ELSE IF ((VERB .EQ. 'xaxis') .OR. (VERB .EQ. 'yaxis')
     &   .OR. (VERB .EQ. 'xyaxis')) THEN
         INVERB = ' '

         IF ((VERB .EQ. 'xaxis') .OR. (VERB .EQ. 'xyaxis')) THEN
            IF ((.NOT. XSCSET) .AND. (IXSCAL .EQ. 'SET')) THEN
               CALL PRTERR ('CMDWARN',
     &            'Automatic scaling selected for X axis')
               IXSCAL = IAXSCA
            END IF
            IF ((.NOT. XTCSET) .AND. (XTICK .NE. 0.0)) THEN
               CALL PRTERR ('CMDWARN',
     &            'Automatic scaling selected for X axis tick interval')
               XTICK = 0.0
            END IF
            IF ((.NOT. XLBSET) .AND. (XLAB .NE. ' ')) THEN
               CALL PRTERR ('CMDWARN',
     &            'Default label selected for X axis')
               XLAB = ' '
            END IF
         END IF

         IF ((VERB .EQ. 'yaxis') .OR. (VERB .EQ. 'xyaxis')) THEN
            IF ((.NOT. YSCSET) .AND. (IYSCAL .EQ. 'SET')) THEN
               CALL PRTERR ('CMDWARN',
     &            'Automatic scaling selected for Y axis')
               IYSCAL = IAXSCA
            END IF
            IF ((.NOT. YTCSET) .AND. (YTICK .NE. 0.0)) THEN
               CALL PRTERR ('CMDWARN',
     &            'Automatic scaling selected for Y axis tick interval')
               YTICK = 0.0
            END IF
            IF ((.NOT. YLBSET) .AND. (YLAB .NE. ' ')) THEN
               CALL PRTERR ('CMDWARN',
     &            'Default label selected for Y axis')
               YLAB = ' '
            END IF
         END IF

C *** Neutral file options ***

      ELSE IF (VERB .EQ. 'ACURVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', CRVNAM)
         CALL FFADDC (CRVNAM, INLINE(1))

      ELSE IF (VERB .EQ. 'NCURVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'curve number', 1, NUMCRV, *100)
         CALL FFADDI (NUMCRV, INLINE(1))
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &      'curve number increment', 1, INCCRV, *100)
         CALL FFADDI (NUMCRV, INLINE(1))

C *** Display options ***

      ELSE IF (VERB .EQ. 'GRID') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFONOF (IFLD, INTYP, CFIELD, DOGRID, *100)

      ELSE IF (VERB .EQ. 'LINES') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
         IF (MATSTR (WORD, 'ON', 2)) THEN
            CALL FFADDC ('ON', INLINE(1))
            LINTYP = 1
         ELSE IF (MATSTR (WORD, 'VARY', 1)) THEN
            CALL FFADDC ('VARY', INLINE(1))
            LINTYP = -1
         ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE(1))
            LINTYP = 0
            IF (ISYTYP .EQ. 0) ISYTYP = -1
         ELSE
            CALL PRTERR ('CMDERR', 'Expected "ON", "VARY" or "OFF"')
            GOTO 100
         END IF

      ELSE IF (VERB .EQ. 'SYMBOLS') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (FFNUMB (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'symbol number', 1, ISY, *100)
            CALL FFADDI (NUMCRV, INLINE(1))
            ISYTYP = ISY
         ELSE
            CALL FFCHAR (IFLD, INTYP, CFIELD, 'ON', WORD)
            IF (MATSTR (WORD, 'ON', 2)) THEN
               CALL FFADDC ('ON', INLINE(1))
               ISYTYP = -1
            ELSE IF (MATSTR (WORD, 'VARY', 1)) THEN
               CALL FFADDC ('VARY', INLINE(1))
               ISYTYP = -1
            ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
               CALL FFADDC ('OFF', INLINE(1))
               ISYTYP = 0
               IF (LINTYP .EQ. 0) LINTYP = 1
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "VARY", "OFF" or symbol number')
               GOTO 100
            END IF
         END IF

      ELSE IF (VERB .EQ. 'CRVNUM') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'LAST', WORD)
         IF (MATSTR (WORD, 'FIRST', 1)) THEN
            CALL FFADDC ('FIRST', INLINE(1))
            LABSID = 'FIRST'
         ELSE IF (MATSTR (WORD, 'MIDDLE', 1)) THEN
            CALL FFADDC ('MIDDLE', INLINE(1))
            LABSID = 'MIDDLE'
         ELSE IF (MATSTR (WORD, 'LAST', 1)) THEN
            CALL FFADDC ('LAST', INLINE(1))
            LABSID = 'LAST'
         ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE(1))
            LABSID = 'NONE'
         ELSE
            CALL PRTERR ('CMDERR',
     &         'Expected "FIRST", "MIDDLE", "LAST" or "OFF"')
            GOTO 100
         END IF

      ELSE IF (VERB .EQ. 'OVERLAY') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, 'VARIABLE', WORD)
         IF (MATSTR(WORD, 'VARIABLE', 1)
     &      .OR. MATSTR (WORD, 'ON', 2)) THEN
            CALL FFADDC ('VARIABLE', INLINE(1))
            OVERLY = .TRUE.
            OVERTM = .FALSE.
         ELSE IF (MATSTR (WORD, 'TIME', 1)) THEN
            CALL FFADDC ('TIME', INLINE(1))
            OVERLY = .FALSE.
            OVERTM = .TRUE.
         ELSE IF (MATSTR (WORD, 'OFF', 3)) THEN
            CALL FFADDC ('OFF', INLINE(1))
            OVERLY = .FALSE.
            OVERTM = .FALSE.
         ELSE
            CALL PRTERR ('CMDERR',
     &         'Expected "VARIABLE", "TIME" or "OFF"')
            GOTO 100
         END IF

      ELSE IF (VERB .EQ. 'SAMESCAL') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *100)
         IF (ISON) THEN
            IAXSCA = 'ALL'
         ELSE
            IAXSCA = 'PLOT'
         END IF

         IF (IXSCAL .NE. 'SET') IXSCAL = IAXSCA
         IF (IYSCAL .NE. 'SET') IYSCAL = IAXSCA

      ELSE IF (VERB .EQ. 'NORMAL') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *100)
         IF (ISON) THEN
            IAXSCA = 'CURVE'
         ELSE
            IAXSCA = 'PLOT'
         END IF

         IF (IXSCAL .NE. 'SET') IXSCAL = IAXSCA
         IF (IYSCAL .NE. 'SET') IYSCAL = IAXSCA

C *** Curve display commands ***

      ELSE IF ((VERB .EQ. 'plot') .OR. (VERB .EQ. 'hardcopy')
     &   .OR. (VERB .EQ. 'neutral')) THEN
         CONTINUE

C *** Mesh display commands ***

      ELSE IF (VERB .EQ. 'MESH') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (.NOT. MESHOK) THEN
            CALL PRTERR ('CMDERR', 'No mesh is defined')
            GOTO 100
         END IF

         CALL LOWSTR (INVERB, VERB)
         VERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL ABRSTR (VERB, WORD, CMDTBL)
         IF (VERB .NE. ' ') THEN
            CALL XYSHOW (VERB)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON)
     &      CALL SHOCMD ('General XY-plot', CMDTBL)
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 110

  100 CONTINUE
      INLINE(1) = ' '

  110 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL XYSHOW (VERB)
      END IF

      RETURN
      END
