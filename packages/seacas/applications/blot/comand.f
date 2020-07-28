C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE COMAND (A, CURPRO,
     &   QAREC, INFREC, NAMECO, NAMELB, NAMES, TIMES, WHOTIM, IPTIMS,
     &   MAPEL, MAPND, IDELB, NEWELB, IELBST, IE2ELB,
     &   LENE, NLNKE, LINKE, XN, YN, ZN, XE, YE, ZE,
     &   ISEVOK,
     &   ISSNPS, IDNPS, NNNPS, ISSESS, IDESS, NEESS, NNESS,
     &   NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &   LISHV, LISGV, LISNV, LISEV, LIDSP, BLKCOL,
     &   NENUM, NEUTRL, NEWPRO, SHDCOL, ISHDCL,
     $   EBNAME, NSNAME, SSNAME, NAMLEN, *)
C=======================================================================

C   --*** COMAND *** (BLOT) Input and process commands
C   --   Modified by John Glick - 11/21/88
C   --   Written by Amy Gilkey - revised 05/20/88
C   --   Dennis Flanagan, 11/18/82
C   --
C   --COMAND inputs and executes command lines.  It returns when a plot
C   --set is defined and the plots are requested.
C   --
C   --The command line uses the 1520 free-field reader.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   CURPRO - IN - the current program name
C   --   QAREC - IN - the QA records
C   --   INFREC - IN - the information records
C   --   NAMECO - IN - the coordinate names
C   --   NAMELB - IN - the element block names
C   --   NAMES - IN - the variable names
C   --   TIMES - IN - the database times
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   IPTIMS - IN/OUT - the selected time steps
C   --   IDELB - IN - the element block ID array
C   --   NEWELB - IN/OUT - the new element blocks flag
C   --      0 = no new element blocks (set elsewhere)
C   --      1 = new selected element blocks
C   --      2 = new displayed element blocks (implies new selected blocks)
C   --   IELBST - IN/OUT - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   IE2ELB - IN - the element block for each element
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the element connectivity
C   --   XN, YN, ZN - IN - the nodal mesh coordinates
C   --   XE, YE, ZE - IN - the element centroid mesh coordinates
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   ISSNPS - IN/OUT - the indices of the selected node sets (SETS)
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   ISSESS - IN/OUT - the indices of the selected side sets (SETS)
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   NCSTEP - IN/OUT - the current step number for display
C   --   LISNP - IN/OUT - the indices of the selected nodes
C   --   NLISEL - IN/OUT - the indices of the selected element blocks
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)
C   --   LISNPS - IN/OUT - the indices of the selected node sets
C   --   LISESS - IN/OUT - the indices of the selected side sets
C   --   LISHV - IN/OUT - the indices of the selected history variables
C   --   LISGV - IN/OUT - the indices of the selected global variables
C   --   LISNV - IN/OUT - the indices of the selected nodal variables
C   --   LISEV - IN/OUT - the indices of the selected element variables
C   --   LIDSP - IN/OUT - the indices of the variables selected for
C   --                    display on the plot legend.
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
C   --   NENUM - IN/OUT - the selected node/element numbers
C   --   NEUTRL - OUT - The type of neutral file to write.
C   --   NEWPRO - IN/OUT - the program indicator:
C   --      'N' = new program, initialize for program
C   --      'R' = new program, initialize for program with reset
C   --      'P' = old program, initialize for plot set
C   --
C   --Common Variables:
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL of /DBNUMS/
C   --   Uses NPTIMS of /TIMES/
C   --   Uses NALVAR of /MSHOPT/

      PARAMETER (MAXFLD = 80)
      include 'params.blk'
      include 'neutral.blk'
      include 'dbnums.blk'
      include 'times.blk'
      include 'mshopt.blk'
      include 'dbase.blk'
      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      DIMENSION A(*)
      CHARACTER*(*) CURPRO
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(MXLNLN) INFREC(*)
      CHARACTER*(MXSTLN) NAMECO(*)
      CHARACTER*(MXSTLN) NAMELB(*)
      CHARACTER*(NAMLEN) NAMES(*)
      CHARACTER*(NAMLEN) EBNAME(*), NSNAME(*), SSNAME(*)
      REAL TIMES(*)
      LOGICAL WHOTIM(*)
      INTEGER IPTIMS(*)
      INTEGER MAPEL(*)
      INTEGER MAPND(*)
      INTEGER IDELB(*)
      INTEGER NEWELB
      INTEGER IELBST(*)
      INTEGER IE2ELB(*)
      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(*)
      REAL XN(*), YN(*), ZN(*)
      REAL XE(*), YE(*), ZE(*)
      LOGICAL ISEVOK(*)
      INTEGER ISSNPS(*)
      INTEGER IDNPS(*), NNNPS(*)
      INTEGER ISSESS(*)
      INTEGER IDESS(*), NEESS(*), NNESS(*)
      INTEGER LISNP(0:*)
      INTEGER NLISEL(0:*)
      INTEGER LISEL(0:*)
      INTEGER LISNPS(0:*), LISESS(0:*)
      INTEGER LISHV(0:*), LISGV(0:*), LISNV(0:*), LISEV(0:*)
      INTEGER LIDSP(0:*)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER NENUM(*)
      INTEGER IDUMA(1)
      REAL SHDCOL(7, NELBLK)
      INTEGER ISHDCL(3, NELBLK)
      CHARACTER NEWPRO
      INTEGER IOSTAT
      SAVE IOSTAT

      LOGICAL MATSTR
      LOGICAL OKABRT

      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*(MXNAME) CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
C      --INTYP, CFIELD, IFIELD, RFIELD - the free-field reader types,
C      --   character fields, integer fields, and real fields

      CHARACTER*(MXNAME) VERB, INVERB
      CHARACTER*1024 INLINE(5)
      CHARACTER*10 PROMPT
      CHARACTER*(MXNAME) WORD
      LOGICAL HELP
      LOGICAL ISON
      LOGICAL OK
      INTEGER NEWMOD(4)

      LOGICAL FIRST
      SAVE FIRST
C      --FIRST - true iff first time through routine

      LOGICAL SAVLOG
      SAVE SAVLOG
C      --SAVLOG - true iff the log file is to be saved

      CHARACTER*(MXSTLN) SAVPRO, MEMPRO
      SAVE SAVPRO, MEMPRO
C      --SAVPRO - the current program to be restored if mesh plot
C      --MEMPRO - the current program used to reserve memory

      LOGICAL MSHTYP, XYTYPE
      SAVE MSHTYP, XYTYPE
C      --MSHTYP - true iff mesh plot program
C      --XYTYPE - true iff XY plot program

      LOGICAL MESHOK
      SAVE MESHOK
C      --MESHOK - true iff a mesh is defined

      CHARACTER*(MXSTLN) CMDTBL(18)
      SAVE CMDTBL
C      --CMDTBL - the command table

      DATA FIRST / .TRUE. /

C   --Command table follows.  The table is ended by a blank entry.
C   --Remember to change the dimensioned size when changing the table.
      DATA (CMDTBL(I),I=1,10) /
     1  'DETOUR                          ',
     *  'PATHLINE                        ',
     *  'SETS                            ',
     *  'TPLOT                           ',
     *  'SPLOT                           ',
     2  'RANGE                           ',
     *  'SELECT                          ',
     *  'LIST                            ',
     *  'PRINT                           ',
     3  'SHOW                            '/
      DATA (CMDTBL(I),I=11,18) /
     *  'HELP                            ',
     *  'SAVELOG                         ',
     *  'LOG                             ',
     4  'RESET                           ',
     *  'EXIT                            ',
     *  'END                             ',
     *  'QUIT                            ',
     5  '                                ' /

      IF (FIRST) THEN
         CURPRO = 'BLOT'
         NEWPRO = '*'
         MSHTYP = .FALSE.
         XYTYPE = .FALSE.

         CALL MSCHK (.FALSE., MESHOK)

C      --Open the log file
         NLOG = 99
         CALL OPNLOG (NLOG)
         SAVLOG = .FALSE.
c        savlog = (nlog .gt. 0)
         if (cdebug .ne. ' ') savlog = (nlog .gt. 0)

         FIRST = .FALSE.
      END IF

      NEUTRL = NONE
C   --Select primary device (ready for any plot)
      CALL GRSPAR ('DEVICE', 0, IDUM, WORD)

C   --Set up to initialize for new plot set or new program

      IF (NEWPRO .EQ. '*') THEN
         INVERB = ' '
      ELSE IF (NEWPRO .EQ. 'N') THEN
         INVERB = 'initprog'
      ELSE IF (NEWPRO .EQ. 'R') THEN
         INVERB = 'initres'
      ELSE
C      --Restore saved program if mesh plot
         IF (CURPRO .EQ. 'MESH') THEN
            CURPRO = SAVPRO
            INVERB = 'postmesh'
         ELSE
            INVERB = 'postplot'
         END IF
      END IF

C   --Reserve memory needed in this routine

      MEMPRO = CURPRO
      IF (MEMPRO .EQ. 'PATHLINE') THEN
         CALL MDRSRV ('LNSCR', KLNSCR, MAX (NUMNP, NUMEL))
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160
      ELSE IF (MEMPRO .EQ. 'TPLOT') THEN
         NTPSCR = MAX (NUMNP, NUMEL)
         CALL MDRSRV ('TPSCR', KTPSCR, 2 * NTPSCR)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160
      END IF

C   --Initialize for new plot set or new program

      NEWPRO = ' '
      IF (INVERB .NE. ' ') THEN
         IFLD = 1
         INTYP(1) = -999
         INLINE(1) = ' '
C      --Initialize the interrupt flag
         ISON = OKABRT (.TRUE.)
         GOTO 110
      END IF

C   --Get next command line and extract verb

  100 CONTINUE
      WRITE (*, *)
      PROMPT = CURPRO
      LPROM = LENSTR (PROMPT) + 2
      PROMPT(LPROM-1:LPROM) = '> '
      CALL GETINS ('parse', MAXFLD, NUMFLD, INTYP, CFIELD,
     &   IFIELD, RFIELD, ' ', IOSTAT, PROMPT, LPROM, *170)
      IF (IOSTAT .NE. 0) THEN
         CALL PRTERR ('CMDWARN', 'Error reading command')
         GOTO 170
      END IF
      IF (NUMFLD .EQ. 0) GOTO 100
      IF (NUMFLD .GT. MAXFLD)
     &   CALL PRTERR ('CMDWARN', 'Too many fields input')
      INTYP(MIN(NUMFLD,MAXFLD)+1) = -999

C   --Initialize the interrupt flag
      ISON = OKABRT (.TRUE.)

      CALL INISTR (5, CHAR(0), INLINE)
      INLINE(1) = ' '
      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', VERB)
      IF (VERB .EQ. ' ') GOTO 150

C   --Get the command verb

C   --Save the command verb - it may be restored later
      INVERB = VERB

      WORD = VERB
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

      IF ((VERB(3:) .NE. ' ') .AND. (WORD(3:) .EQ. ' ')) THEN
         CALL PRTERR ('CMDERR',
     &      'Enter at least 3 letters for BLOT-level commands')
         GOTO 150
      END IF

C   --These commands are valid at any level, including BLOT

C   --Commands to switch programs

      IF (VERB .EQ. 'SETS') THEN
         CALL PRTERR ('CMDERR',
     &      'This function has been moved into DETOUR')
         GOTO 150

      ELSE IF ((VERB .EQ. 'DETOUR') .OR. (VERB .EQ. 'PATHLINE')
     &   .OR. (VERB .EQ. 'TPLOT') .OR. (VERB .EQ. 'SPLOT')) THEN
C     Commented out code that forces the user to type the whole
C     program name on the command line to transfer to another program
c         IF ((CURPRO .NE. 'BLOT') .AND. (CURPRO .NE. VERB)
c     &      .AND. (INVERB .NE. VERB)) THEN
c            CALL PRTERR ('CMDERR',
c     &         'Program name cannot be abbreviated')
c            GOTO 150
c         END IF
         CALL FFADDC (VERB, INLINE(1))

         IF (CURPRO .NE. VERB) THEN
            IF (VERB .EQ. 'DETOUR') THEN
               CALL DTCHK (.TRUE., OK)
            ELSE IF (VERB .EQ. 'PATHLINE') THEN
               CALL LNCHK (.TRUE., OK)
            ELSE IF (VERB .EQ. 'TPLOT') THEN
               CALL TPCHK (.TRUE., OK)
            ELSE IF (VERB .EQ. 'SPLOT') THEN
               CALL SPCHK (.TRUE., OK)
            END IF
            IF (.NOT. OK) GOTO 150

            MSHTYP = (VERB .EQ. 'DETOUR') .OR. (VERB .EQ. 'PATHLINE')
            XYTYPE = (VERB .EQ. 'TPLOT') .OR. (VERB .EQ. 'SPLOT')

            CURPRO = VERB(1:8)
            VERB = ' '

            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               CALL FFADDC ('RESET', INLINE(1))
               NEWPRO = 'R'
            ELSE
               NEWPRO = 'N'
            END IF

         ELSE
            VERB = ' '
         END IF

C   --Commands to exit

      ELSE IF ((VERB .EQ. 'EXIT') .OR. (VERB .EQ. 'END')
     &   .OR. (VERB .EQ. 'QUIT')) THEN
         CALL FFADDC (VERB, INLINE(1))
         VERB = ' '

         CALL SCNEOF
         NEWPRO = 'E'

C   --Informational commands

      ELSE IF ((VERB .EQ. 'LOG') .or. (verb .eq. 'SAVELOG')) THEN
         if (verb .eq. 'SAVELOG') then
            call prterr ('CMDSPEC', 'Please use the LOG command')
            verb = 'LOG'
         end if
         IF (NLOG .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Log file cannot be opened')
            GOTO 150
         END IF
         VERB = ' '
         SAVLOG = .TRUE.
         IF (SAVLOG) THEN
            WRITE (*, 10000) 'Log file will be saved'
         ELSE
            WRITE (*, 10000) 'Log file will NOT be saved'
         END IF

      ELSE IF (VERB .EQ. 'SELECT') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)

         CALL DBSEL (A, A, INLINE,
     &      WORD, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NAMES, TIMES, WHOTIM, NPTIMS, IPTIMS,
     &      IDELB, LENE, IDNPS, IDESS,
     &      NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &      LISHV, LISGV, LISNV, LISEV, MAPEL, MAPND)
         VERB = ' '
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160

      ELSE IF ((VERB .EQ. 'LIST') .OR. (VERB .EQ. 'PRINT')) THEN
         CALL DBLIST (A, INLINE,
     &      VERB, IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NAMECO, NAMELB, NAMES, QAREC, INFREC,
     &      TIMES, WHOTIM, NPTIMS, IPTIMS, XN, YN, ZN,
     &      IDELB, LENE, NLNKE, LINKE, ISEVOK,
     &      IDNPS, NNNPS, IDESS, NEESS, NNESS,
     &      NCSTEP, LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &      LISHV, LISGV, LISNV, LISEV, EBNAME, NSNAME, SSNAME, NAMLEN,
     *      MAPEL, MAPND)
         VERB = ' '
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160

      ELSE IF (VERB .EQ. 'RANGE') THEN
         CALL FFADDC (VERB, INLINE(1))
         VERB = ' '

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL FFADDC (WORD, INLINE(1))
         IVAR = LOCSTR (WORD, NVARHI+NVARGL+NVARNP+NVAREL, NAMES)
         IF (IVAR .LE. 0) THEN
            CALL PRTERR ('CMDERR', 'Expected variable name')
            GOTO 150
         END IF

         IF (MSHTYP) THEN
            CALL SCALER (A, A, 2, NAMES(IVAR), IVAR,
     &         .TRUE., IELBST, NALVAR, DUMMIN, DUMMAX, MAPEL, MAPND)
         ELSE
            CALL SCALER (A, A, 2, NAMES(IVAR), IVAR,
     &         .FALSE., IDUMA, 0, DUMMIN, DUMMAX, MAPEL, MAPND)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160

      ELSE IF ((VERB .EQ. 'SHOW') .AND. (CURPRO .NE. 'BLOT')) THEN
         CALL LOWSTR (INVERB, VERB)

      ELSE IF (VERB .EQ. 'HELP') THEN
         ISON = HELP ('BLOT', 'COMMANDS', CFIELD(IFLD))
         IF (.NOT. ISON) THEN
            CALL SHOCMD ('BLOT Commands', CMDTBL)
            CALL LOWSTR (INVERB, VERB)
         ELSE
            VERB = ' '
         END IF

C   --Reset command

      ELSE IF (VERB .EQ. 'RESET') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL LOWSTR (INVERB, VERB)

      else if (verb .eq. 'DEBUG') then
         call ffaddc (verb, inline(1))
         verb = ' '
         call ffchar (ifld, intyp, cfield, ' ', cdebug)
         call ffaddc (cdebug, inline(1))
         call ffintg (ifld, intyp, ifield,
     &      'debug constant', 0, idebug, *150)
         call ffaddi (idebug, inline(1))

      ELSE IF (CURPRO .EQ. 'BLOT') THEN
         VERB = ' '
         CALL PRTERR ('CMDERR', 'Expected a subprogram name')
         GOTO 150
      END IF

      IF (VERB .EQ. ' ') GOTO 120

  110 CONTINUE

C   --Restore the command verb
      VERB = INVERB

C   --Perform general command

      IF (MSHTYP .OR. XYTYPE) THEN
         IIFLD = IFLD
         CALL PLCOMD (A, CURPRO, XYTYPE, MESHOK, INLINE,
     &      VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      TIMES, WHOTIM, IPTIMS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160
         IF (VERB .EQ. ' ') GOTO 120
      END IF

      IF (MSHTYP) THEN

C      --Perform MESH command

         IIFLD = IFLD
         CALL MSCOMD (A, CURPRO, MSHTYP, INLINE,
     &      VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      NAMECO, NAMES, IDELB, NEWELB, IELBST, NEWMOD,
     &      IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL, SHDCOL, ISHDCL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160
         IF (VERB .EQ. ' ') GOTO 120

C      --Perform DETOUR command

         IF (CURPRO .EQ. 'DETOUR') THEN
            IIFLD = IFLD
            CALL DTCOMD (A, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NEWMOD, NAMES, IELBST,
     &         ISSNPS, ISSESS, LIDSP, MAPEL, MAPND, NAMLEN)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF

C      --Perform PATHLINE command

         IF (CURPRO .EQ. 'PATHLINE') THEN
            IIFLD = IFLD
            CALL LNCOMD (A, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, ISEVOK, IE2ELB, A(KLNSCR))
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF

      ELSE IF (XYTYPE) THEN

C      --Perform TPLOT command

         IF (CURPRO .EQ. 'TPLOT') THEN
            IIFLD = IFLD
            CALL TPCOMD (A, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD, MAXFLD,
     &         NAMES, ISEVOK, IE2ELB, A(KTPSCR), NTPSCR, MAPEL, MAPND)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF

C      --Perform SPLOT command

         IF (CURPRO .EQ. 'SPLOT') THEN
            IIFLD = IFLD
            CALL SPCOMD (A, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMES, LENE, NLNKE, LINKE, XN, YN, ZN, XE, YE, ZE,
     &         ISEVOK, IE2ELB, NENUM, LIDSP, MAPEL, MAPND, NAMLEN)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF

C      --Perform XY command

         IIFLD = IFLD
         CALL XYCOMD (A, CURPRO, INLINE,
     &      VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &      MESHOK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160
         IF (VERB .EQ. ' ') GOTO 120

C      --Perform MESH command

         IF (MESHOK) THEN
            IIFLD = IFLD
            CALL MSCOMD (A, CURPRO, MSHTYP, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         NAMECO, NAMES, IDELB, NEWELB, IELBST, NEWMOD,
     &         IDNPS, ISSNPS, IDESS, ISSESS, BLKCOL, SHDCOL, ISHDCL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF
      END IF

C   --Make second pass through PLCOMD on reset to set options that may be
C   --controlled by xxCOMD after PLCOMD is called the first time

      IF (MSHTYP .OR. XYTYPE) THEN
         IF ((VERB .EQ. 'initprog') .OR. (VERB .EQ. 'initres')
     &      .OR. (VERB .EQ. 'reset')) THEN
            VERB = 'postinit'
            IIFLD = 0
            CALL PLCOMD (A, CURPRO, XYTYPE, MESHOK, INLINE,
     &         VERB, IIFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &         TIMES, WHOTIM, IPTIMS)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 160
            IF (VERB .EQ. ' ') GOTO 120
         END IF
      END IF

C   --Check command verb

      IF (VERB .EQ. 'show') THEN
         CALL PRTERR ('CMDERR', 'Invalid SHOW option')
         VERB = ' '

      ELSE IF ((VERB .EQ. 'plot') .OR. (VERB .EQ. 'hardcopy') .OR.
     *     (VERB .EQ. 'neutral')  .OR. (VERB .EQ. 'grafaid')  .OR.
     *     (VERB .EQ. 'csv')      .OR. (VERB .EQ. 'xmgr')     .OR.
     *     (VERB .EQ. 'raw')      .OR. (VERB .EQ. 'mesh')) THEN

C      --Set neutral file flag
        IF (VERB .EQ. 'neutral') neutrl = XMGR
        IF (VERB .EQ. 'xmgr')    neutrl = XMGR
        IF (VERB .EQ. 'grafaid') neutrl = GRAF
        IF (VERB .EQ. 'csv')     neutrl = CSV
        IF (VERB .EQ. 'raw')     neutrl = RAW

C      --Switch to hardcopy device
         IF (VERB .EQ. 'hardcopy') THEN
            CALL GRGPARD ('DEVICE', 2, ISON, WORD)
            IF (.NOT. ISON) VERB = 'plot'
         END IF
         IF (VERB .EQ. 'hardcopy')
     &      CALL GRSPAR ('DEVICE', 2, IDUM, WORD)

C      --If mesh display, save current program to be restored
         IF (VERB .EQ. 'mesh') THEN
            SAVPRO = CURPRO
            CURPRO = 'MESH'
         END IF

         NEWPRO = 'P'
         VERB = ' '

      ELSE IF ((VERB(1:1) .GE. 'a') .AND. (VERB(1:1) .LE. 'z')) THEN
         VERB = ' '
      END IF

      IF (VERB .EQ. ' ') GOTO 120
      CALL PRTERR ('CMDERR',
     &   '"' // VERB(:LENSTR(VERB)) // '" is an invalid command')
      VERB = ' '
      GOTO 150

  120 CONTINUE
      IF ((NLOG .GT. 0) .AND. (INLINE(1) .NE. ' ')) THEN
         DO 130 I = 1, 5
            IF (INLINE(I) .EQ. CHAR(0)) GOTO 140
            WRITE (NLOG, '(A)') INLINE(I)(:LENSTR(INLINE(I)))
  130    CONTINUE
  140    CONTINUE
      END IF

  150 CONTINUE
      IF (NEWPRO .EQ. ' ') GOTO 100

C   --Release memory needed in this routine

      IF (MEMPRO .EQ. 'PATHLINE') THEN
         CALL MDDEL ('LNSCR')
      ELSE IF (MEMPRO .EQ. 'TPLOT') THEN
         CALL MDDEL ('TPSCR')
      END IF

      IF (NEWPRO .EQ. 'E') GOTO 170

C   --Return to plot or change programs

  160 CONTINUE
      RETURN

C   --Exit program

  170 CONTINUE

C   --Close the log file (delete if not saved)

      IF (NLOG .GT. 0) THEN
         IF (SAVLOG) THEN
            CLOSE (NLOG, IOSTAT=IDUM)
         ELSE
            CLOSE (NLOG, STATUS='DELETE', IOSTAT=IDUM)
         END IF
      END IF

      RETURN 1
10000  FORMAT (1X, 5A)
      END
