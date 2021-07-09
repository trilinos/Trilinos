C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLCOMD (A, CURPRO, XYTYPE, MESHOK, INLINE,
     &           INVERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &           TIMES, WHOTIM, IPTIMS)
C=======================================================================

C   --*** PLCOMD *** (BLOT) Process general graphics commands
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --PLCOMD interprets general graphics commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Time step control
C   --   TMIN     {tmin}           Set minimum selected time
C   --   TMAX     {tmax}           Set maximum selected time
C   --   DELTIME  {delt}           Set selected time interval
C   --   NINTV    {nintv}          Change selected time interval
C   --   ZINTV    {nintv}          Change selected time interval (zero interval)
C   --   ALLTIMES                  Select all times
C   --   TIMES    {t1,t2...}       Select specified times
C   --   STEPS    {n1,n2...}       Select specified steps
C   --   HISTORY  {ON/OFF}         Allow history steps to be selected
C   --
C   --Display control
C   --   QA       {ON/OFF}         Set plot QA information on label flag
C   --   AXIS     {ON/OFF}         Set axis numbering flag
C   --   LEGEND   {ON/OFF}         Set legend information on label flag
C   --   CAPTION  {line_number}    Set plot caption (on following lines)
C   --   SOFTCHAR {ON/OFF,device}  Select software or hardware characters
C   --   FONT     {STICK/SANSERIF/ROMAN,device} Set font
C   --   COLOR    {numcol}         Set number of colors
C   --   SPECTRUM {numcol}         Set number of spectrum colors
C   --   SNAP     {nsnap,device}   Set number of frames to snap for movies
C   --   AUTO     {ON/OFF,device}  Select automatic or user-directed plotting
C   --
C   --Display
C   --   reset                     Reset default conditions
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   postinit                  Initialize for program change (second pass)
C   --
C   --Information
C   --   show     {option}         Display program information (based on type)
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   CURPRO - IN - the current program name
C   --   XYTYPE - IN - true iff current program is an XY curve versus mesh plot
C   --   MESHOK - IN - true iff mesh can be displayed
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INVERB - IN/OUT - the command verb
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   IPTIMS - IN/OUT - the selected time steps
C   --
C   --Common Variables:
C   --   Uses NSTEPS of /DBNUMS/
C   --   Sets DOQA, DOLEG, DOAXIS, CAPTN of /LEGOPT/
C   --   Sets and uses NPTIMS, NPTIMW, TMIN, TMAX, DELT, NINTV, NTIMIN,
C   --      WHONLY, HISTOK of /TIMES/

      EXTERNAL BLKDAT

      include 'params.blk'
      include 'dbnums.blk'
      include 'times.blk'
      include 'legopt.blk'
      include 'icrnbw.blk'

      DIMENSION A(*)
      CHARACTER*(*) CURPRO
      LOGICAL XYTYPE, MESHOK
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)

      CHARACTER*(MXSTLN) VERB
C      --VERB - the command verb used in the SHOW

      CHARACTER*(MXSTLN) WORD
      LOGICAL HELP
      LOGICAL ISON

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(24)
      SAVE CMDTBL
C      --CMDTBL - the command table

      DATA FIRST /.TRUE./

C   --Command table follows.  The table is ended by a blank entry.
C   --Remember to change the dimensioned size when changing the table.
      DATA (CMDTBL(I),I=1,10) /
     1  'TMIN                            ',
     *  'TMAX                            ',
     *  'DELTIME                         ',
     *  'NINTV                           ',
     *  'ZINTV                           ',
     2  'ALLTIMES                        ',
     *  'TIMES                           ',
     *  'STEPS                           ',
     *  'HISTORY                         ',
     3  'QA                              '/

      DATA (CMDTBL(I),I=11,20) /
     *  'AXIS                            ',
     *  'LEGEND                          ',
     *  'CAPTION                         ',
     *  'OUTLINE                         ',
     4  'SOFTCHARACTERS                  ',
     *  'FONT                            ',
     *  'COLOR                           ',
     *  'SPECTRUM                        ',
     *  'SNAP                            ',
     5  'AUTO                            '/

      DATA (CMDTBL(I),I=21,24) /
     *  'BACKGROUND                      ',
     *  'FOREGROUND                      ',
     *  'RAINBOW                         ',
     6  '                                ' /

C   --Get the command verb

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

C *** Initialization ***

      IF ((VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
     &   .OR. (VERB .EQ. 'initres') .OR. (VERB .EQ. 'reset')
     &   .OR. (VERB .EQ. 'postinit')) THEN
         INVERB = VERB

C      --Initialize parameters first time through, then reset

         IF (FIRST) THEN

C         --Set default time step selection to 10 times with delta interval
            HISTOK = .FALSE.
            WHONLY = .NOT. HISTOK
            NTIMIN = 10

C         --Set to one time step if no time steps
            IF (NSTEPS .LE. 0) THEN
               TIMES(1) = 0.0
               WHOTIM(1) = .TRUE.
               NPTIMS = 1
               NPTIMW = 1
               IPTIMS(1) = 1
               DELT = 0.0
               NINTV = 1
            END IF

            VERB = 'initres'
         END IF

C      --Initialize for program change

         IF (VERB .EQ. 'initprog') THEN

C         --Change the time step selection when 'postinit'

C         --Set display options

C         --DOQA, DOAXIS, and DOLEG are not changed

c            CALL INISTR (3, ' ', CAPTN(1,1))
c            CALL INISTR (3, ' ', CAPTN(1,2))
         END IF

C      --Initialize for program change or program reset (part 2)

         IF (VERB .EQ. 'postinit') THEN

C         --Set time step selection to selected default

            WHONLY = .NOT. HISTOK
            IF (NSTEPS .GT. 0) THEN
               IF ((.NOT. WHONLY) .OR. (NSTEPW .GT. 0)) THEN
                  CALL INITIM (NTIMIN, WHONLY, NSTEPS, TIMES, WHOTIM,
     &               TMIN, TMAX, DELT, NINTV, NPTIMS, IPTIMS)
                  NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)
               ELSE
                  NPTIMS = 1
                  NPTIMW = 1
                  IPTIMS(1) = 1
                  DELT = 0.0
                  NINTV = 1
               END IF
            END IF
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN

C         --Change the time step selection when 'postinit'

C         --Set time step selection according to default

C         --Set display options

            DOQA(1) = .TRUE.
            DOQA(2) = .TRUE.
            DOAXIS(1) = (NDIM .LE. 2)
            DOAXIS(2) = .TRUE.
            DOLEG(1) = .TRUE.
            DOLEG(2) = .TRUE.
            DOBOX    = .TRUE.

            if (verb .EQ. 'reset') then
              CALL INISTR (3, ' ', CAPTN(1,1))
              CALL INISTR (3, ' ', CAPTN(1,2))
            end if

C              reset background color

            CALL SETBCK (3, INLINE, IFLD, INTYP, CFIELD, *100)

C              reset foreground color

            CALL SETFOR (2, INLINE, IFLD, INTYP, CFIELD, *100)

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

         IF (VERB .NE. 'postplot') THEN
            CONTINUE
         END IF

         VERB = ' '

C *** Time step control ***

      ELSE IF ((VERB .EQ. 'TMIN') .OR. (VERB .EQ. 'TMAX')
     &   .OR. (VERB .EQ. 'DELTIME')
     &   .OR. (VERB .EQ. 'NINTV') .OR. (VERB .EQ. 'ZINTV')
     &   .OR. (VERB .EQ. 'ALLTIMES')
     &   .OR. (VERB .EQ. 'TIMES') .OR. (VERB .EQ. 'STEPS')) THEN
         INVERB = ' '
         CALL CKNONE (NSTEPS, .FALSE., 'time steps', *100)
         IF (WHONLY) THEN
            CALL CKNONE (NSTEPW, .FALSE., 'whole time steps', *100)
         END IF

         CALL CMDTIM (INLINE,
     &      VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NSTEPS, TIMES, WHOTIM,
     &      WHONLY, TMIN, TMAX, DELT, NINTV, NPTIMS, IPTIMS)
         NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)

      ELSE IF (VERB .EQ. 'HISTORY') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         CALL CKNONE (NSTEPS, .FALSE., 'time steps', *100)
         IF (.NOT. HISTOK) THEN
            CALL PRTERR ('CMDERR',
     &         'History time steps cannot be referenced'
     &         // ' by this subprogram')
            GOTO 100
         END IF

         CALL FFONOF (IFLD, INTYP, CFIELD, ISON, *100)
         CALL FFADDO (ISON, INLINE(1))
         WHONLY = .NOT. ISON

         IF ((.NOT. WHONLY) .OR. (NSTEPW .GT. 0)) THEN
            NPTIMS = -999
            CALL CALTIM (WHONLY, TMIN, TMAX, DELT, NINTV,
     &         NSTEPS, TIMES, WHOTIM, NPTIMS, IPTIMS)
            NPTIMW = NWHSEL (NPTIMS, IPTIMS, WHOTIM)
         ELSE
            CALL PRTERR ('CMDWARN',
     &         'There are no whole time steps defined')
            NPTIMS = 1
            NPTIMW = 1
            IPTIMS(1) = 1
            DELT = 0.0
            NINTV = 1
         END IF

         CALL PLSHOW (VERB, XYTYPE, MESHOK, TIMES, WHOTIM, IPTIMS)
         VERB = 'TMIN'

      ELSE IF ((VERB .EQ. 'QA') .OR. (VERB .EQ. 'LEGEND')
     &   .OR. (VERB .EQ. 'AXIS') .OR. (VERB .EQ. 'CAPTION')
     *   .OR. (VERB .EQ. 'OUTLINE')) THEN
         INVERB = ' '

         CALL CMDLEG (INLINE, VERB, IFLD, INTYP, CFIELD, IFIELD,
     &        XYTYPE, MESHOK, DOQA, DOLEG, DOAXIS, DOBOX, CAPTN, *100)

      ELSE IF ((VERB .EQ. 'SOFTCHARACTERS') .OR. (VERB .EQ. 'FONT')
     &   .OR. (VERB .EQ. 'COLOR') .OR. (VERB .EQ. 'SPECTRUM')
     &   .OR. (VERB .EQ. 'SNAP') .OR. (VERB .EQ. 'AUTO')
     &   .OR. (VERB .EQ. 'RAINBOW')) THEN

C  Check for RAINBOW command, which will be reset to 'SPECTRUM'
         IF (VERB .EQ. 'RAINBOW')THEN
C  RAINBOW sets the flag for using the rainbow spectrum rather
C  than the default blue-brown-red spectrum
          IRAINB=1
          CALL PRTERR('CMDSPEC','Using rainbow color scale.')
          VERB='SPECTRUM'
         ELSEIF (VERB .EQ. 'SPECTRUM')THEN
          IRAINB=0
         ENDIF

         INVERB = ' '

         CALL CMDDEV (INLINE,
     &      VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD, *100)

C      --COLOR and SPECTRUM are to be passed as lower-case commands
         IF ((VERB .EQ. 'COLOR') .OR. (VERB .EQ. 'SPECTRUM'))
     &      CALL LOWSTR (INVERB, VERB)

      ELSE IF (VERB .EQ. 'BACKGROUND') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         CALL SETBCK (1, INLINE, IFLD, INTYP, CFIELD, *100)

      ELSE IF (VERB .EQ. 'FOREGROUND') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '
         CALL SETFOR (1, INLINE, IFLD, INTYP, CFIELD, *100)

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL ABRSTR (VERB, WORD, CMDTBL)
         IF (VERB .NE. ' ') THEN
            CALL PLSHOW (VERB, XYTYPE, MESHOK, TIMES, WHOTIM, IPTIMS)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON)
     &      CALL SHOCMD ('General Graphics Commands', CMDTBL)
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 110

  100 CONTINUE
      INLINE(1) = ' '

  110 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL PLSHOW (VERB, XYTYPE, MESHOK, TIMES, WHOTIM, IPTIMS)
      END IF

      RETURN
      END
