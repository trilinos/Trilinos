C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SPCOMD (A, INLINE,
     &   INVERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   NAMES, LENE, NLNKE, LINKE, XN, YN, ZN, XE, YE, ZE,
     &   ISEVOK, IE2ELB, NENUM, LIDSP, MAPEL, MAPND, NAMLEN)
C=======================================================================

C   --*** SPCOMD *** (SPLOT) Process SPLOT commands
C   --   Modified by John Glick - 11/1/88
C   --   Written by Amy Gilkey - revised 05/31/88
C   --
C   --SPCOMD interprets SPLOT commands.
C   --
C   --The commands are listed below with their parameters and their
C   --function.
C   --
C   --Plot curve selection
C   --   ADD                       Save plot curves, otherwise a new plot set
C   --                             is started if a curve is defined
C   --   REMOVE   {n}              Delete plot curve n
C   --   SYPLOT   {variable}       Specify plot variable
C   --   variable                  Specify time plot
C   --   NODES    {# range}        Select node numbers
C   --   ELEMENTS {# range}        Select element numbers
C   --
C   --Display Control
C   --   DISPVAR [ADD] {var1,...}  Selects history variables, global
C   --                             variables, and/or TIME whose values
C   --                             will be displayed on the plot legend.
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
C   --Mesh Display
C   --   ECHO                      Plot mesh with selected nodes/elements
C   --
C   --Information
C   --   show     {option}         Display plot parameters and information
C   --   help     {option}         Display system dependent HELP
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   IPATH - optimization information for FNDPTH
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   INLINE - IN/OUT - the parsed input lines for the log file
C   --   INVERB - IN/OUT - the command verb
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   NAMES - IN - the variable names
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the element connectivity
C   --   XN, YN, ZN - IN - the nodal mesh coordinates
C   --   XE, YE, ZE - IN - the element centroid mesh coordinates
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   IE2ELB - IN - the element block for each element
C   --   NENUM - IN/OUT - the selected node/element numbers
C   --   LIDSP(0:*)  - IN/OUT - the indices of the selected variables
C   --          whose values will be displayed on the plot legend.
C   --          LIDSP(0) = the number of variables in the list.
C   --          LIDSP(i) identifies the ith variable in the list.
C   --          If LIDSP(i) > 0, LIDSP(i) is the id of a history variable.
C   --          If LIDSP(i) < 0, -LIDSP(i) is the id of a global variable.
C   --          If LIDSP(i) = 0, TIME is to be displayed on the plot legend.
C   --
C   --Common Variables:
C   --   Uses NVARNP, NVAREL of /DBNUMS/
C   --   Uses NPTIMS of /TIMES/
C   --   Sets NTIMIN, HISTOK of /TIMES/
C   --   Sets NODVAR, NNENUM of /SELNE/
C   --   Sets NSPVAR of /SPVARS/
C   --   Sets OVERLY, OVERTM of /XYOPT/

      include 'params.blk'
      include 'dbnums.blk'
      include 'times.blk'
      include 'selne.blk'
      include 'spvars.blk'
      include 'xyopt.blk'

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      CHARACTER*(*) INVERB
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      CHARACTER*(*) NAMES(*)
      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(NELBLK)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      INTEGER IE2ELB(*)
      INTEGER NENUM(*)
      INTEGER LIDSP(0:*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)

      LOGICAL FFEXST, FFMATC

      CHARACTER*(MXNAME) VERB, WORD
      CHARACTER*5 STRA
      CHARACTER TYP
      LOGICAL HELP
      LOGICAL ADDCRV
      LOGICAL ISON, ISPATH

      LOGICAL FIRST
      SAVE FIRST, ADDCRV
C      --FIRST - true iff first time through routine

      CHARACTER*(MXSTLN) CMDTBL(15)
      SAVE CMDTBL
C      --CMDTBL - the valid commands table

      DATA FIRST /.TRUE./

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1  'ADD                             ',
     *  'REMOVE                          ',
     *  'delete                          ',
     2  'SYPLOT                          ',
     *  'NODES                           ',
     *  'ELEMENTS                        ',
     *  'DISPVAR                         ',
     3  'PLOT                            ',
     *  'HARDCOPY                        ',
     *  'NEUTRAL                         ',
     *  'GRAFAID                         ',
     *  'XMGR                            ',
     *  'CSV                             ',
     *  'ECHO                            ',
     4  '                                ' /

C   --Get the command verb, which may be a variable name

      WORD = INVERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') THEN
         CALL DBVIX_BL ('N', 1, INV)
         CALL DBVIX_BL ('E', 1, IEV)
         IF ((LOCSTR (WORD, NVARNP, NAMES(INV)) .GT. 0)
     &      .OR. (LOCSTR (WORD, NVAREL, NAMES(IEV)) .GT. 0)) THEN
            VERB = 'SYPLOT'
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

C         --Change the command table to upper case

            L = LOCSTR (' ', 999, CMDTBL) - 1
            DO 100 I = 1, L
               CALL EXUPCS (CMDTBL(I))
  100       CONTINUE

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
            CONTINUE
         END IF

C      --Reset parameters

         IF ((VERB .EQ. 'reset') .OR. (VERB .EQ. 'initres')) THEN

            IF (VERB .EQ. 'reset') THEN
C           --reset display variable list
               CALL DISPV (.TRUE., INLINE, IFLD, INTYP,
     &            CFIELD, NAMES, LIDSP, NAMLEN)
            ENDIF

            SELOK = .FALSE.
            NODVAR = .TRUE.
            NNENUM = 0
            NSPVAR = 0

            FIRST = .FALSE.
         END IF

C      --Initialize for new plot set

         IF (VERB .EQ. 'postplot') THEN
            CONTINUE
         END IF

C      --Initialize after mesh plot

         IF (VERB .EQ. 'postmesh') THEN
            SELOK = .FALSE.
         END IF

         ADDCRV = .FALSE.

         VERB = ' '

C *** Plot variable selection ***

      ELSE IF (VERB .EQ. 'ADD') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         VERB = 'PLOT'
         ADDCRV = .TRUE.

      ELSE IF (VERB .EQ. 'REMOVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         VERB = 'PLOT'

  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'number of the curve to delete', 0, ICRV, *130)
            IF (ICRV .LT. 0) THEN
               NDEL = MIN (-ICRV, NSPVAR)
               CALL FFADDI (-NDEL, INLINE(1))
               NSPVAR = NSPVAR - NDEL
            ELSE
               IF ((ICRV .LT. 1) .OR. (ICRV .GT. NSPVAR)) THEN
                  CALL INTSTR (1, 0, ICRV, STRA, L)
                  CALL PRTERR ('CMDWARN', 'Curve number ' // STRA(:L) //
     &               ' is invalid, ignored')
                  GOTO 120
               END IF
               CALL FFADDI (ICRV, INLINE(1))
               ISVID(ICRV) = 0
            END IF
  120       CONTINUE
            GOTO 110
         END IF

  130    CONTINUE
         ICRV = 0
         DO 140 NP = 1, NSPVAR
            IF (ISVID(NP) .GT. 0) THEN
               ICRV = ICRV + 1
               IF (ICRV .NE. NP) THEN
                  ISVID(ICRV) = ISVID(NP)
               END IF
            END IF
  140    CONTINUE
         NSPVAR = ICRV

      else if (verb .eq. 'DELETE') then
         call prterr ('CMDREQ', 'Please use the REMOVE command')

         call ffaddc (verb, inline(1))
         inverb = ' '

         verb = 'PLOT'
         call ffintg (ifld, intyp, ifield,
     &      'number of curves to delete', 1, ndel, *180)
         ndel = max (0, min (ndel, nspvar))
         call ffaddi (ndel, inline(1))
         nspvar = nspvar - ndel

      ELSE IF (VERB .EQ. 'SYPLOT') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (.NOT. ADDCRV) THEN
            NSPVAR = 0
            ADDCRV = .TRUE.
         END IF

         IF (NNENUM .EQ. 0) THEN
            CALL PRTERR ('CMDERR',
     &         'A NODES/ELEMENTS command must be issued')
            GOTO 180
         END IF

C      --Get variable name

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL FFADDC (WORD, INLINE(1))

         IVAR = LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
         IF (IVAR .LE. 0) THEN
            CALL PRTERR ('CMDERR', '"' // WORD(:LENSTR(WORD))
     &         // '" is an invalid variable name')
            GOTO 180
         END IF

C      --Get variable type

         CALL DBVTYP_BL (IVAR, TYP, IDUM)

         IF (NODVAR) THEN
            IF (TYP .NE. 'N') THEN
               CALL PRTERR ('CMDERR',
     &            'Only nodal variables may be specified')
               GOTO 180
            END IF
         ELSE
            IF (TYP .NE. 'E') THEN
               CALL PRTERR ('CMDERR',
     &            'Only element variables may be specified')
               GOTO 180
            END IF
         END IF

C      --Add variable to the plot variable list

         IF (NSPVAR .GE. MXSVAR) THEN
            CALL PRTERR ('CMDERR', 'Too many plot variables')
            GOTO 180
         END IF
         NSPVAR = NSPVAR + 1
         ISVID(NSPVAR) = IVAR

C      --Set Y axis to automatic scaling, etc if not set
         INVERB = 'yaxis'

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
               CALL RXINTA (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &            'node number', NUMNP, NNENUM, NENUM, MAPND, *180)
            ELSE
               CALL RXINTA (INLINE(1), IFLD, INTYP, CFIELD, IFIELD,
     &            'element number', NUMEL, NNENUM, NENUM, MAPEL, *180)
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
            IF (NERR .GT. 0) GOTO 200

            IF (NODVAR) THEN
               CALL FNDPTH (NODVAR, LENE, NLNKE, LINKE, XN, YN, ZN,
     &            A(KIPATH), NUMNP, NNENUM, NENUM)
            ELSE
               CALL FNDPTH (NODVAR, LENE, NLNKE, LINKE, XE, YE, ZE,
     &            A(KIPATH), NUMEL, NNENUM, NENUM)
            END IF
         END IF

C      --Set X axis to automatic scaling, etc if not set
         INVERB = 'xaxis'

C *** Display control ***

      ELSE IF (VERB .EQ. 'DISPVAR') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL DISPV (.FALSE., INLINE, IFLD, INTYP,
     &      CFIELD, NAMES, LIDSP, NAMLEN)
         INVERB = ' '

C *** Display ***

      ELSE IF ((VERB .EQ. 'PLOT') .OR. (VERB .EQ. 'HARDCOPY') .OR.
     &         (VERB .EQ. 'NEUTRAL') .OR. (VERB .EQ. 'GRAFAID') .OR.
     *         (VERB .EQ. 'XMGR')    .OR. (VERB .EQ. 'CSV')) THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (NNENUM .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No nodes or elements are selected')
            GOTO 180
         END IF
         IF (NSPVAR .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No curves are defined')
            GOTO 180
         END IF
         IF (NPTIMS .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No times are selected')
            VERB = ' '
            GOTO 180
         END IF

C      --Check that at least one selected element is defined for each selected
C      --element variable

         IF (.NOT. NODVAR) THEN
            NBAD = 0
            DO 170 IP = 1, NSPVAR
               CALL DBVTYP_BL (ISVID(IP), TYP, IVAR)
               IF (NUMEQL (.FALSE., NELBLK, ISEVOK(1,IVAR)) .GT. 0) THEN
                  DO 150 IX = 1, NNENUM
                     IF (ISEVOK(IE2ELB(NENUM(IX)),IVAR)) GOTO 160
  150             CONTINUE
                  NBAD = NBAD + 1
  160             CONTINUE
               END IF
  170       CONTINUE
            IF (NBAD .GT. 0) THEN
               IF (NBAD .GE. NSPVAR) THEN
                  CALL PRTERR ('CMDERR',
     &               'All of the curve variables are undefined'
     &               // ' for the requested elements')
                  GOTO 180
               ELSE
                  CALL PRTERR ('CMDERR',
     &               'Some of the curve variables are undefined'
     &               // ' for the requested elements')
               END IF
            END IF
         END IF

C      --PLOT, HARDCOPY, and NEUTRAL are to be passed as lower-case commands
         CALL LOWSTR (INVERB, VERB)
         VERB = ' '

C *** Mesh display commands ***

      ELSE IF (VERB .EQ. 'ECHO') THEN
         CALL FFADDC (VERB, INLINE(1))
         INVERB = ' '

         IF (NNENUM .EQ. 0) THEN
            CALL PRTERR ('CMDERR', 'No nodes or elements are selected')
            GOTO 180
         END IF
         SELOK = .TRUE.

         INVERB = 'mesh'
         VERB = ' '

C *** Information ***

      ELSE IF (VERB .EQ. 'show') THEN

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL ABRSTR (VERB, WORD, CMDTBL)
         IF (VERB .NE. ' ') THEN
            CALL SPSHOW (VERB, NAMES, NENUM, LIDSP)
            INVERB = ' '
         END IF
         VERB = ' '

      ELSE IF (VERB .EQ. 'help') THEN

         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON)
     &      CALL SHOCMD ('SPLOT Commands', CMDTBL)
         VERB = ' '

      ELSE
         VERB = ' '
      END IF

      GOTO 190

  180 CONTINUE
      INLINE(1) = ' '

  190 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL SPSHOW (VERB, NAMES, NENUM, LIDSP)
      END IF

  200 CONTINUE
      RETURN
      END
