C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LNCOMD (A, INLINE,
     &   INVERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NAMES, ISEVOK, IE2ELB, NENUM)
C=======================================================================
C TODO: Support for Mapped node/element ids

C   --*** LNCOMD *** (PATHLN) Process PATHLINE commands
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --LNCOMD interprets PATHLINE commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Plot curve selection
C   --   ADD                       Save plot curves, otherwise a new plot set
C   --                             is started if a curve is requested
C   --   REMOVE   {n}              Delete plot curve n
C   --   LOCATION {xvar,yvar,zvar,#} Specify pathline plot
C   --
C   --Display
C   --   reset                     Reset plot parameters
C   --   postplot                  Initialize after plot
C   --   initprog                  Initialize for program change
C   --   initres                   Initialize for program change and reset
C   --   plot                      Exit to plot the plot set
C   --   hardcopy                  Exit to plot the plot set on hardcopy device
C   --
C   --Information
C   --   show     {option}         Display plot parameters and information
C   --   help     {option}         Display system dependent HELP
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INVERB - IN/OUT - the command verb
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NAMES - IN - the global, nodal, and element variable names
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   IE2ELB - IN - the element block for each element
C   --   NENUM - SCRATCH - size = MAX (NUMNP, NUMEL)
C   --
C   --Common Variables:
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NPTIMS, NPTIMW of /TIMES/
C   --   Sets NTIMIN, HISTOK of /TIMES/
C   --   Sets NLNCRV, ILVNE, ILVID of /LNVARS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'times.blk'
      include 'lnvars.blk'

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMES(*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      INTEGER IE2ELB(*)
      INTEGER NENUM(*)

      LOGICAL FFEXST, FFMATC

      CHARACTER*(MXNAME) VERB, WORD
      LOGICAL HELP
      LOGICAL ADDCRV
      LOGICAL ISON
      LOGICAL ISPATH, NODVAR
      CHARACTER*5 STRA
      CHARACTER*(MXNAME) NAME
      CHARACTER CDUM
      CHARACTER TYP, T
      INTEGER IID(3)

      LOGICAL FIRST
      SAVE FIRST, ADDCRV
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(4)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST / .TRUE. /

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1   'ADD     ', 'REMOVE  ', 'LOCATION  ',
     2   '        ' /

C   --Get the command verb, which may be a variable name

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') THEN
         IF (LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES) .GT. 0)
     &      THEN
            VERB = 'LOCATION'
            IFLD = IFLD - 1
         ELSE
            VERB = WORD
         END IF
      END IF

C *** Initialization ***

      IF ((VERB .EQ. 'postplot') .OR. (VERB .EQ. 'initprog')
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
            NLNCRV = 0

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

         ADDCRV = .FALSE.

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
               NDEL = MIN (-ICRV, NLNCRV)
               CALL FFADDI (-NDEL, INLINE(1))
               NLNCRV = NLNCRV - NDEL
            ELSE
               IF ((ICRV .LT. 1) .OR. (ICRV .GT. NLNCRV)) THEN
                  CALL INTSTR (1, 0, ICRV, STRA, L)
                  CALL PRTERR ('CMDWARN', 'Curve number ' // STRA(:L) //
     &               ' is invalid, ignored')
                  GOTO 110
               END IF
               CALL FFADDI (ICRV, INLINE(1))
               ILVID(1,ICRV) = 0
            END IF
  110       CONTINUE
            GOTO 100
         END IF

  120    CONTINUE
         ICRV = 0
         DO 140 NP = 1, NLNCRV
            IF (ILVID(1,NP) .GT. 0) THEN
               ICRV = ICRV + 1
               IF (ICRV .NE. NP) THEN
                  ILVNE(ICRV) = ILVNE(NP)
                  DO 130 I = 1, NDIM
                     ILVID(I,ICRV) = ILVID(I,NP)
  130             CONTINUE
               END IF
            END IF
  140    CONTINUE
         NLNCRV = ICRV

      ELSE IF (VERB .EQ. 'LOCATION') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

C      --Delete existing curves, if not adding and set to adding

         IF (.NOT. ADDCRV) THEN
            NLNCRV = 0
            ADDCRV = .TRUE.
         END IF

         DO 150 IXY = 1, NDIM

C         --Get the variable name

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', NAME)
            IF (NAME .EQ. ' ') THEN
               CALL PRTERR ('CMDERR', 'Expected variable name')
               GOTO 200
            END IF
            CALL FFADDC (NAME, INLINE(1))
            IID(IXY) = LOCSTR (NAME, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
            IF (IID(IXY) .LE. 0) THEN
               CALL PRTERR ('CMDERR', '"' // NAME(:LENSTR(NAME))
     &            // '" is an invalid variable name')
               GOTO 200
            END IF

C         --Get the variable type

            CALL DBVTYP_BL (IID(IXY), T, IDUM)

            IF (IXY .EQ. 1) THEN
               TYP = T
               IF (TYP .EQ. 'E') THEN
                  CALL PRTERR ('CMDERR',
     &               'Pathlines not allowed on element variables')
                  GOTO 200
               END IF
            END IF
            IF (TYP .NE. T) THEN
               CALL PRTERR ('CMDERR',
     &            'Variables are of different types')
               GOTO 200
            END IF

            IF ((NSTEPW .LE. 0) .AND. (TYP .NE. 'H')) THEN
               CALL PRTERR ('CMDERR',
     &            'Only history variables are defined on this database')
               GOTO 200
            END IF
  150    CONTINUE

C      --Get the variable range

         IF ((TYP .EQ. 'N') .OR. (TYP .EQ. 'E')) THEN
            NNUM = 0
            IF (FFEXST (IFLD, INTYP)) THEN
               IF (TYP .EQ. 'N') THEN
                  CALL RIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &               'node number', NUMNP, NNUM, NENUM, *200)
               ELSE IF (TYP .EQ. 'E') THEN
                  CALL RIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &               'element number', NUMEL, NNUM, NENUM, *200)
               END IF
            END IF
            IF (NNUM .LE. 0) THEN
               CALL PRTERR ('CMDERR', 'Expected node/element number')
               GOTO 200
            END IF
         ELSE
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (WORD .NE. ' ') THEN
               CALL PRTERR ('CMDERR',
     &            'A number is given for a history/global variable')
               GOTO 200
            END IF
            NNUM = 1
         END IF

C      --Add the curves to the list (deleted if error)

         NPTMP = NLNCRV
         DO 180 NP = 1, NNUM

C         --Check that element variable exists for element block

            IF (TYP .EQ. 'E') THEN
               IEL = NENUM(NP)
               IELB = IE2ELB(IEL)
               DO 160 IXY = 1, NDIM
                  CALL DBVTYP_BL (IID(IXY), CDUM, ID)
                  IF (.NOT. ISEVOK (IELB, ID)) THEN
                     NAME = NAMES(IID(IXY))
                     CALL INTSTR (1, 0, IEL, STRA, LSTRA)
                     CALL PRTERR ('CMDERR',
     &                  'Variable ' // NAME(:LENSTR(NAME)) //
     &                  ' does not exist for element ' // STRA(:LSTRA)
     &                  // ', curve ignored')
                     GOTO 180
                  END IF
  160          CONTINUE
            END IF

C         --Add the curve to the list (deleted if error)

            IF ((NPTMP + 1) .GT. MXLCRV) THEN
               CALL PRTERR ('CMDERR', 'Too many plot variables')
               GOTO 200
            END IF

            NPTMP = NPTMP + 1
            ILVNE(NPTMP) = NENUM(NP)
            DO 170 IXY = 1, NDIM
               ILVID(IXY,NPTMP) = IID(IXY)
  170       CONTINUE
  180    CONTINUE

         IF (NLNCRV .EQ. NPTMP) GOTO 200

         NLNCRV = NPTMP

      ELSE IF ((VERB .EQ. 'NODES') .OR. (VERB .EQ. 'ELEMENTS')) THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         ISON = (VERB .EQ. 'NODES')
C      --Change the verb to lower-case for SHOW
         WORD = VERB
         CALL LOWSTR (VERB, WORD)

         ISPATH = FFMATC (IFLD, INTYP, CFIELD, 'PATH', 1)
         IF (ISPATH) CALL FFADDC ('PATH', INLINE(1))

C      --Make a list of the selected nodes

         IF (FFEXST (IFLD, INTYP)) THEN
            IF (ISON) THEN
               CALL RIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &            'node number', NUMNP, NNENUM, NENUM, *200)
            ELSE
               CALL RIXINT (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &            'element number', NUMEL, NNENUM, NENUM, *200)
            END IF
         END IF

         IF (ISON .NEQV. NODVAR) NSPVAR = 0
         NODVAR = ISON

C      --Determine a path between the nodes/element, if requested

         IF (ISPATH) THEN
C         --Get memory for path efficiency information
            CALL MDFIND ('IPATH', KIPATH, L)
            IF (L .LE. 0) THEN
               CALL MDLONG ('IPATH', KIPATH, 2*NUMNP)
            END IF
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 220

            IF (NODVAR) THEN
               CALL FNDPTH (NODVAR, LENE, NLNKE, LINKE, XN, YN, ZN,
     &            A(KIPATH), NUMNP, NNENUM, NENUM)
            ELSE
               CALL FNDPTH (NODVAR, LENE, NLNKE, LINKE, XE, YE, ZE,
     &            A(KIPATH), NUMEL, NNENUM, NENUM)
            END IF
         END IF

C *** Display ***

      ELSE IF ((VERB .EQ. 'plot') .OR. (VERB .EQ. 'hardcopy')) THEN
         CALL FFADDC (VERB, INLINE(1))

         IF (NLNCRV .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No curves are defined')
            GOTO 200
         END IF
         IF (NPTIMS .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No times are selected')
            VERB = ' '
            GOTO 200
         END IF

C      --Check that at least one whole time step is selected for non-history
C      --variables

         IF (NPTIMW .LE. 0) THEN
            NBAD = 0
            DO 190 NP = 1, NLNCRV
               CALL DBVTYP_BL (ILVID(1,NP), TYP, IDUM)
               IF (TYP .NE. 'H') THEN
                  NBAD = NBAD + 1
               END IF
  190       CONTINUE
            IF (NBAD .GT. 0) THEN
               IF (NBAD .GE. NLNCRV) THEN
                  CALL PRTERR ('CMDERR',
     &               'All of the curves are undefined'
     &               // ' for the selected time steps')
                  GOTO 200
               ELSE
                  CALL PRTERR ('CMDERR',
     &               'Some of the curves are undefined'
     &               // ' for the selected time steps')
               END IF
            END IF
         END IF

         VERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         IF ((WORD .NE. 'PLOT') .AND. (WORD .NE. 'HARDCOPY')
     &      .AND. (WORD .NE. 'VIEW')) THEN
            CALL ABRSTR (VERB, WORD, CMDTBL)
         ELSE
            VERB = WORD
         END IF
         IF (VERB .NE. ' ') THEN
            CALL LNSHOW (VERB, NAMES)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON)
     &      CALL SHOCMD ('PATHLINE Commands', CMDTBL)
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 210

  200 CONTINUE
      INLINE(1) = ' '

  210 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL LNSHOW (VERB, NAMES)
      END IF

  220 CONTINUE
      RETURN
      END
