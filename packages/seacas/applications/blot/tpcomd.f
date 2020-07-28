C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPCOMD (A, INLINE,
     &   INVERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD, MAXFLD,
     &   NAMES, ISEVOK, IE2ELB, NENUM, NNENUM, MAPEL, MAPND)
C=======================================================================

C   --*** TPCOMD *** (TPLOT) Process TPLOT commands
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --TPCOMD interprets TPLOT commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Plot curve selection
C   --   ADD                       Save plot curves, otherwise a new plot set
C   --                             is started if a curve is requested
C   --   REMOVE   {n}              Delete plot curve n
C   --   TYPLOT                    Specify time plot
C   --            (variable, # on following line)
C   --   XYPLOT                    Specify X-Y plot
C   --            (variable, # on following two lines)
C   --   variable, #               Specify time plot
C   --
C   --Display
C   --   reset                     Reset plot parameters
C   --   postmesh                  Initialize after mesh plot
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   PLOT                      Exit to plot the plot set
C   --   HARDCOPY                  Exit to plot the plot set on hardcopy device
C   --   NEUTRAL                   Exit to write the plot set to neutral file
C   --
C   --Information
C   --   show     {option}         Display plot parameters and information
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INVERB - IN/OUT - the command verb
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD, MAXFLD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NAMES - IN - the variable names
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   IE2ELB - IN - the element block for each element
C   --   NENUM - SCRATCH - size = 2 * NNENUM
C   --   NNENUM - used to properly dimension NENUM
C   --
C   --Common Variables:
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NPTIMS, NPTIMW of /TIMES/
C   --   Sets NTIMIN, HISTOK of /TIMES/
C   --   Sets NTPCRV, NTPVAR, TIMPLT, ITVID of /TPVARS/
C   --   Sets OVERLY of /XYOPT/

      include 'params.blk'
      include 'dbnums.blk'
      include 'times.blk'
      include 'tpvars.blk'
      include 'xyopt.blk'

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMES(*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      INTEGER NENUM(NNENUM,2)
      INTEGER IE2ELB(*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)

      LOGICAL FFEXST
      LOGICAL BATCH

      CHARACTER*(MXNAME) VERB, WORD
      LOGICAL HELP
      LOGICAL ADDCRV
      LOGICAL ISON
      LOGICAL TYPLOT
      CHARACTER*5 STRA
      CHARACTER*(MXNAME) NAME
      CHARACTER CDUM
      CHARACTER TYP(2)
      INTEGER IID(2)
      CHARACTER TYPX, TYPY
      CHARACTER*12 PROMPT
      INTEGER LPROM

      LOGICAL FIRST
      SAVE FIRST, ADDCRV
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(13)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST /.TRUE./
      DATA LASTX /-1/
      DATA LASTY /-1/

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1  'ADD                             ',
     *  'REMOVE                          ',
     *  'delete                          ',
     *  'TYPLOT                          ',
     *  'XYPLOT                          ',
     2  'PLOT                            ',
     *  'HARDCOPY                        ',
     *  'NEUTRAL                         ',
     *  'XMGR                            ',
     *  'CSV                             ',
     *  'RAW                             ',
     *  'GRAFAID                         ',
     3  '                                ' /

C   --Get the command verb, which may be a variable name

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') THEN
         IF (LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES) .GT. 0)
     &      THEN
            VERB = 'TYPLOT'
            IFLD = IFLD - 1
         ELSE
            VERB = WORD
         END IF
      END IF

C *** Initialization ***

      IF ((VERB .EQ. 'postmesh')
     &   .OR. (VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
     &   .OR. (VERB .EQ. 'initres') .OR. (VERB .EQ. 'reset')) THEN
         INVERB = VERB

C      --Initialize parameters first time through, then reset

         IF (FIRST) THEN
            VERB = 'initres'
         END IF

C      --Initialize time step selection default

         IF ((VERB .EQ. 'initprog') .OR. (VERB .EQ. 'initres')) THEN

C         --Set default time step selection to all times
            HISTOK = (NVARHI .GT. 0)
            NTIMIN = 0
         END IF

C      --Initialize for program change

         IF (VERB .EQ. 'initprog') THEN
            CONTINUE
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN
            NTPVAR = 0
            NTPCRV = 0
            TIMPLT = .TRUE.

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

C      --Initialize after mesh plot

         IF (VERB .EQ. 'postmesh') THEN
            CONTINUE
         END IF

         ADDCRV = .FALSE.
         IF (NTPVAR .GT. 0) THEN
            IF (TIMPLT) THEN
               LASTX = 0
            ELSE
               LASTX = ITVID(NTPVAR-1)
            END IF
            LASTY = ITVID(NTPVAR)
         ELSE
            LASTX = -1
            LASTY = -1
         END IF

         VERB = ' '

C *** Plot curve selection ***

      ELSE IF (VERB .EQ. 'ADD') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         VERB = 'PLOT'
         ADDCRV = .TRUE.

      ELSE IF (VERB .EQ. 'REMOVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         VERB = 'PLOT'

  100    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'number of the curve to delete', 0, ICRV, *120)
            IF (ICRV .LT. 0) THEN
               NDEL = MIN (-ICRV, NTPCRV)
               CALL FFADDI (-NDEL, INLINE(1))
               IF (TIMPLT) THEN
                  NTPVAR = NTPVAR - NDEL
               ELSE
                  NTPVAR = NTPVAR - 2*NDEL
               END IF
               NTPCRV = NTPCRV - NDEL
            ELSE
               IF ((ICRV .LT. 1) .OR. (ICRV .GT. NTPCRV)) THEN
                  CALL INTSTR (1, 0, ICRV, STRA, L)
                  CALL PRTERR ('CMDWARN', 'Curve number ' // STRA(:L) //
     &               ' is invalid, ignored')
                  GOTO 110
               END IF
               CALL FFADDI (ICRV, INLINE(1))
               IF (TIMPLT) THEN
                  ITVID(ICRV) = 0
               ELSE
                  ITVID(2*ICRV-1) = 0
                  ITVID(2*ICRV) = 0
               END IF
            END IF
  110       CONTINUE
            GOTO 100
         END IF

  120    CONTINUE
         N = 0
         ICRV = 0
         DO 130 NP = 1, NTPCRV
            N = N + 1
            IF (.NOT. TIMPLT) N = N + 1
            IF (ITVID(N) .GT. 0) THEN
               ICRV = ICRV + 1
               IF (ICRV .NE. NP) THEN
                  IF (TIMPLT) THEN
                     ITVID(ICRV) = ITVID(N)
                     ITVNE(ICRV) = ITVNE(NP)
                  ELSE
                     ITVID(2*ICRV-1) = ITVID(N-1)
                     ITVNE(2*ICRV-1) = ITVNE(N-1)
                     ITVID(2*ICRV) = ITVID(N)
                     ITVNE(2*ICRV) = ITVNE(N)
                  END IF
               END IF
            END IF
  130    CONTINUE
         IF (TIMPLT) THEN
            NTPVAR = ICRV
         ELSE
            NTPVAR = 2*ICRV
         END IF
         NTPCRV = ICRV

      else if (verb .eq. 'DELETE') then
         call prterr ('CMDREQ', 'Please use the REMOVE command')

         call ffaddc (verb, inline(1))
         inverb = ' '

         verb = 'PLOT'
         call ffintg (ifld, intyp, ifield,
     &      'number of curves to delete', 1, ndel, *210)
         ndel = max (0, min (ndel, ntpcrv))
         call ffaddi (ndel, inline(1))
         if (timplt) then
            ntpvar = ntpvar - ndel
         else
            ntpvar = ntpvar - 2*ndel
         end if
         ntpcrv = ntpcrv - ndel

      ELSE IF ((VERB .EQ. 'TYPLOT') .OR. (VERB .EQ. 'XYPLOT')) THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

C      --Delete existing curves, if not adding and set to adding

         IF (.NOT. ADDCRV) THEN
            NTPVAR = 0
            NTPCRV = 0
            ADDCRV = .TRUE.
         END IF

C      --Determine if X-Y plot or time-Y plot

         TYPLOT = (VERB .EQ. 'TYPLOT')
         IF ((TIMPLT .NEQV. TYPLOT) .AND. (NTPCRV .GT. 0)) THEN
            CALL PRTERR ('CMDERR',
     &         'Time curves and X-Y curves must be defined separately')
            GOTO 210
         END IF

         IF (TYPLOT) THEN
            ISTART = 2
         ELSE
            ISTART = 1
         END IF

         NOK = ISTART - 1

         DO 160 IXY = ISTART, 2
            IF (ISTART .EQ. 1) THEN
               ILIN = IXY + 1
            ELSE
               ILIN = IXY
            END IF
            INLINE(ILIN) = ' '

C         --Input the variable line

            IF ((IXY .GT. ISTART)
     &         .OR. (.NOT. FFEXST (IFLD, INTYP))) THEN
               IF (IXY .EQ. 1) THEN
                  PROMPT = 'X VARIABLE> '
                  LPROM = LENSTR (PROMPT)
                  CALL GETINS ('parse', MAXFLD, NUMFLD, INTYP,
     &               CFIELD, IFIELD, RFIELD, ' ', IOSTAT, PROMPT,
     &               LPROM, *140)
               ELSE
                  PROMPT = 'Y VARIABLE> '
                  LPROM = LENSTR (PROMPT)
                  CALL GETINS ('parse', MAXFLD, NUMFLD, INTYP,
     &               CFIELD, IFIELD, RFIELD, ' ', IOSTAT, PROMPT,
     &               LPROM, *140)
               END IF
  140          CONTINUE
               IF (IOSTAT .LT. 0) NUMFLD = 0
               INTYP(MIN(NUMFLD,MAXFLD)+1) = -999
               IFLD = 1
            END IF

C         --Get the variable name

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)
            CALL FFADDC (NAME, INLINE(ILIN))
            IID(IXY) = LOCSTR (NAME, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
            IF (IID(IXY) .LE. 0) THEN
               CALL PRTERR ('CMDERR', '"' // NAME(:LENSTR(NAME))
     &            // '" is an invalid variable name')
               GOTO 150
            END IF

C         --Get the variable type

            CALL DBVTYP_BL (IID(IXY), TYP(IXY), ID)

            IF ((NSTEPW .LE. 0) .AND. (TYP(IXY) .NE. 'H')) THEN
               CALL PRTERR ('CMDERR',
     &            'Only history variables are defined on this database')
               GOTO 150
            END IF

C         --Get the variable range

            IF ((TYP(IXY) .EQ. 'N') .OR. (TYP(IXY) .EQ. 'E')) THEN
               NNUM = 0
               IF (FFEXST (IFLD, INTYP)) THEN
                  IF (TYP(IXY) .EQ. 'N') THEN
                     CALL RMIXINT (INLINE(ILIN),
     &                  IFLD, INTYP, CFIELD, IFIELD,
     &                  'node number', NUMNP, NNUM, NENUM(1,IXY),
     &                  MAPND, *150)
                  ELSE IF (TYP(IXY) .EQ. 'E') THEN
                     CALL RMIXINT (INLINE(ILIN),
     &                IFLD, INTYP, CFIELD, IFIELD,
     &                'element number', NUMEL, NNUM, NENUM(1,IXY),
     &                MAPEL, *150)
                  END IF
               END IF
               IF (NNUM .LE. 0) THEN
                  CALL PRTERR ('CMDERR', 'Expected node/element number')
                  GOTO 150
               END IF
            ELSE
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (WORD .NE. ' ') THEN
                  CALL PRTERR ('CMDERR',
     &               'A number is given for a history/global variable')
                  GOTO 150
               END IF
               NNUM = 1
            END IF

            IF (IXY .EQ. 1) THEN
               NX = NNUM
            ELSE
               NY = NNUM
            END IF

            NOK = NOK + 1
  150       CONTINUE
            IF ((.NOT. BATCH ()) .AND. (NOK .LT. IXY)) GOTO 210
  160    CONTINUE
         IF (NOK .NE. 2) GOTO 210

C      --Get the number of X and Y parts

         IF (ISTART .EQ. 2) NX = NY
         IF ((TYP(2) .EQ. 'N') .OR. (TYP(2) .EQ. 'E')) THEN
            IF ((NX .NE. NY) .AND. (NX .NE. 1) .AND. (NY .NE. 1)) THEN
               CALL PRTERR ('CMDERR',
     &            'Ranges of X and Y variables do not match')
               GOTO 210
            END IF
         END IF

C      --Add the curves to the list (deleted if error)

         NPTMP = NTPVAR
         DO 190 IENT = 1, MAX (NX, NY)

C         --Check that element variable exists for element block

            DO 170 IXY = ISTART, 2
               IF (TYP(IXY) .EQ. 'E') THEN
                  IEL = NENUM(IENT,IXY)
                  IELB = IE2ELB(IEL)
                  CALL DBVTYP_BL (IID(IXY), CDUM, ID)
                  IF (.NOT. ISEVOK (IELB, ID)) THEN
                     NAME = NAMES(IID(IXY))
                     CALL INTSTR (1, 0, IEL, STRA, LSTRA)
                     CALL PRTERR ('CMDWARN',
     &                  'Variable ' // NAME(:LENSTR(NAME)) //
     &                  ' does not exist for element ' // STRA(:LSTRA)
     &                  // ', curve ignored')
                     GOTO 190
                  END IF
               END IF
  170       CONTINUE

C         --Add the curve to the list

            IF ((NPTMP + (2 - ISTART + 1)) .GT. MXTVAR) THEN
               CALL PRTERR ('CMDERR', 'Too many plot variables')
               GOTO 210
            END IF

            DO 180 IXY = ISTART, 2
               NPTMP = NPTMP + 1

               ITVID(NPTMP) = IID(IXY)
               IF ((TYP(IXY) .EQ. 'N') .OR. (TYP(IXY) .EQ. 'E')) THEN
                  ITVNE(NPTMP) = NENUM(IENT,IXY)
               ELSE
                  ITVNE(NPTMP) = 0
               END IF
  180       CONTINUE
  190    CONTINUE

         IF (NTPVAR .EQ. NPTMP) GOTO 210

         NTPVAR = NPTMP
         TIMPLT = TYPLOT
         NTPCRV = NTPCRV + MAX (NX, NY)

C      --Reset the X and Y axis if the variables have changed

         IF (TYPLOT) IID(1) = 0
         NEWX = IID(1)
         NEWY = IID(2)

         IF (LASTX .NE. NEWX) THEN
            INVERB = 'xaxis'
            LASTX = NEWX
         END IF
         IF (LASTY .NE. NEWY) THEN
            IF (INVERB .EQ. ' ') THEN
               INVERB = 'yaxis'
            ELSE
               INVERB = 'xyaxis'
            END IF
            LASTY = NEWY
         END IF

C *** Display ***

      ELSE IF ((VERB .EQ. 'PLOT') .OR. (VERB .EQ. 'HARDCOPY') .OR.
     &         (VERB .EQ. 'NEUTRAL') .OR. (VERB .EQ. 'GRAFAID') .OR.
     *         (VERB .EQ. 'XMGR')    .OR. (VERB .EQ. 'CSV') .OR.
     *         (VERB .EQ. 'RAW' )) THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (NTPVAR .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No curves are defined')
            GOTO 210
         END IF
         IF (NPTIMS .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No times are selected')
            VERB = ' '
            GOTO 210
         END IF

C      --Check that at least one whole time step is selected for non-history
C      --variables

         IF (NPTIMW .LE. 0) THEN
            NBAD = 0
            IF (TIMPLT) TYPX = 'H'
            N = 0
            DO 200 NP = 1, NTPCRV
               IF (.NOT. TIMPLT) THEN
                  N = N + 1
                  CALL DBVTYP_BL (ITVID(N), TYPX, IDUM)
               END IF
               N = N + 1
               CALL DBVTYP_BL (ITVID(N), TYPY, IDUM)
               IF ((TYPX .NE. 'H') .OR. (TYPY .NE. 'H')) THEN
                  NBAD = NBAD + 1
               END IF
  200       CONTINUE
            IF (NBAD .GT. 0) THEN
               IF (NBAD .GE. NTPCRV) THEN
                  CALL PRTERR ('CMDERR',
     &               'All of the curves are undefined'
     &               // ' for the selected time steps')
                  GOTO 210
               ELSE
                  CALL PRTERR ('CMDERR',
     &               'Some of the curves are undefined'
     &               // ' for the selected time steps')
               END IF
            END IF
         END IF

C      --PLOT, HARDCOPY, and NEUTRAL are to be passed as lower-case commands
         CALL LOWSTR (INVERB, VERB)
         VERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL ABRSTR (VERB, WORD, CMDTBL)
         IF (VERB .NE. ' ') THEN
            CALL TPSHOW (VERB, NAMES, MAPEL, MAPND)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON)
     &      CALL SHOCMD ('TPLOT Commands', CMDTBL)
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 220

  210 CONTINUE
      INLINE(1) = ' '

  220 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL TPSHOW (VERB, NAMES, MAPEL, MAPND)
      END IF

      RETURN
      END
