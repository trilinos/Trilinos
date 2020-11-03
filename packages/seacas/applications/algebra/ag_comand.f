C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE COMAND (A, INLINE, INTYP, CFIELD, IFIELD, RFIELD,
     &   NAMECO, BLKTYP, NAMES, TIMES, IPTIMS,
     &   IDELB, VISELB, SELELB,
     &   QAREC, INFREC, RETVRB, MERR)
C=======================================================================
C   --*** COMAND *** (ALGEBRA) Process command
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --COMAND processes an ALGEBRA command.  The commands are:
C   --
C   --   TITLE                     Change database title
C   --
C   --   SAVE     {N,E,G,A,var,...} Save input database variable
C   --   DELETE   {var,var,...}    Delete assigned variable
C   --   ALIAS    {name,var,n}     Assign alias to vector/tensor variables
C   --
C   --   TMIN     {TMIN}           Set minimum selected time
C   --   TMAX     {TMAX}           Set maximum selected time
C   --   DELTIME  {DELT}           Set selected time interval
C   --   NINTV    {NINTV}          Change selected time interval
C   --   ZINTV    {NINTV}          Change selected time interval(zero interval)
C   --   ALLTIMES                  Select all times
C   --   TIMES    {t1,t2...}       Select specified times
C   --   STEPS    {n1,n2...}       Select specified steps
C   --
C   --   ZOOM     {x1,x2,y1,y2,z1,z2} Set zoomed mesh limits
C   --   VISIBLE  {id1,...}        Set element blocks to be written
C   --   BLOCKS   {id1,...}        Set selected element blocks
C   --   MATERIAL {id1,...}        Set selected element blocks
C   --
C   --   LIST     {option}         Display database information
C   --   SHOW     {option}         Display parameter setting
C   --   HELP     {option}         Help on ALGEBRA
C   --
C   --   END                       Exit to evaluate equations (also EXIT)
C   --   QUIT                      Exit but do not write database
C   --
C   --Parameters:
C   --   A         - IN - the dynamic memory base array
C   --   INLINE    - IN/OUT - the parsed input lines for the log file
C   --   INTYP     - IN - the input types from the free field reader
C   --   CFIELD    - IN - the character fields
C   --   IFIELD    - IN - the integer fields
C   --   RFIELD    - IN - the real fields
C   --   NAMECO    - IN - the coordinate names
C   --   BLKTYP    - IN - the element block names
C   --   NAMES     - IN - the global, nodal, and element variable names
C   --   TIMES     - IN - the database time steps
C   --   IPTIMS    - IN/OUT - the selected times steps
C   --   IDELB     - IN - the element block IDs
C   --   VISELB(i) - IN/OUT - true iff element block i is to be written
C   --   SELELB(i) - IN/OUT - true iff element block i is selected
C   --   QAREC     - IN - the QA records containing:
C   --              (1) - the analysis code name
C   --              (2) - the analysis code QA descriptor
C   --              (3) - the analysis date
C   --              (4) - the analysis time
C   --   INFREC    - IN - the information records
C   --   RETVRB    - OUT - returned action verb (END, QUIT, PRINT, LOG)
C   --
C   --Common Variables:
C   --   Uses NDIM, NSTEPS of /DBNUMS/
C   --   Sets TITLEO of /DBTITL/
C   --   Sets NPTIMS, TMIN, TMAX, DELT, NINTV, WHONLY of /TIMES/
C   --   Sets ISZOOM, ZMLIM of /ZOOM/
C   --   Uses FNCNAM of /FNCTBL/

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbtitl.blk'
      include 'ag_times.blk'
      include 'ag_zoom.blk'
      include 'ag_filter.blk'
      include 'ag_remove.blk'
      include 'ag_fnctbc.blk'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*(mxstln) cdebug

      LOGICAL FFEXST, FFNUMB, MATSTR

      DIMENSION A(*)
      CHARACTER*(*) INLINE(*)
      INTEGER       INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER       IFIELD(*)
      REAL          RFIELD(*)
      CHARACTER*(namlen) NAMECO(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(namlen) NAMES(*)
      REAL TIMES(*)
      INTEGER IPTIMS(*)
      INTEGER IDELB(*)
      LOGICAL VISELB(NELBLK)
      LOGICAL SELELB(NELBLK)
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(MXLNLN) INFREC(*)
      CHARACTER*(*) RETVRB
      CHARACTER*1 OPTION

      CHARACTER*(maxnam) VERB, WORD
      CHARACTER*(mxstln) HLPTYP
      LOGICAL ISON
      LOGICAL HELP
      LOGICAL XLIM, YLIM, ZLIM
      INTEGER MERR

      CHARACTER*(mxstln) CMDTBL(26)
      SAVE CMDTBL
C      --CMDTBL - the commands table

      CHARACTER*(mxstln) HLPTBL(4)
      SAVE HLPTBL
C      --HLPTBL - the HELP type table

      MERR = 0

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1  'TITLE                           ',
     2  'SAVE                            ',
     3  'DELETE                          ',
     4  'ALIAS                           ',
     5  'TMIN                            ',
     6  'TMAX                            ',
     7  'DELTIME                         ',
     8  'NINTV                           ',
     9  'ZINTV                           ',
     *  'ALLTIMES                        ',
     1  'TIMES                           ',
     2  'STEPS                           ',
     3  'ZOOM                            ',
     4  'FILTER                          ',
     5  'VISIBLE                         ',
     6  'BLOCKS                          ',
     7  'MATERIAL                        ',
     8  'LOG                             ',
     9  'LIST                            ',
     *  'SHOW                            ',
     1  'HELP                            ',
     2  'END                             ',
     3  'EXIT                            ',
     4  'QUIT                            ',
     4  'REMOVE                          ',
     *  '                                ' /

C   --HELP type table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA HLPTBL /
     &  'RULES                           ',
     *  'COMMANDS                        ',
     *  'FUNCTION                        ',
     &  '                                ' /

      RETVRB = ' '

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF (VERB .EQ. ' ') VERB = WORD

      IF (VERB .EQ. 'TITLE') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL GETINP (0, 0, 'TITLE> ', TITLEO, IOSTAT)
         INLINE(2) = TITLEO

      ELSE IF (VERB .EQ. 'ALIAS') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL ALICMD (INLINE, INTYP(IFLD), CFIELD(IFLD), IFIELD(IFLD),
     &      NAMES, *150)

      ELSE IF (VERB .EQ. 'SAVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL SAVCMD (INLINE, INTYP(IFLD), CFIELD(IFLD), NAMES, *150)

      ELSE IF (VERB .EQ. 'DELETE') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL DELCMD (INLINE, INTYP(IFLD), CFIELD(IFLD), *150)

      ELSE IF ((VERB .EQ. 'TMIN') .OR. (VERB .EQ. 'TMAX')
     &   .OR. (VERB .EQ. 'DELTIME')
     &   .OR. (VERB .EQ. 'NINTV') .OR. (VERB .EQ. 'ZINTV')
     &   .OR. (VERB .EQ. 'ALLTIMES')
     &   .OR. (VERB .EQ. 'TIMES') .OR. (VERB .EQ. 'STEPS')) THEN
         CALL CMDTIM (INLINE, VERB, IFLD, INTYP, CFIELD, IFIELD,
     &                RFIELD, NSTEPS, TIMES, TMIN, TMAX, DELT,
     &                NINTV, NPTIMS, IPTIMS)
      ELSE IF (VERB .EQ. 'FILTER') THEN
         CALL FFADDC (VERB, INLINE(1))
C ... Get name of element variable to use for filtering
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         if (.not. matstr(word, 'ELEMENTS', 1)) THEN
           CALL PRTERR ('CMDERR',
     &       'Expected "ELEMENTS".  Syntax: FILTER ELEMENTS {elem_var}')
           VERB = ' '
           GOTO 150
         ENDIF
         CALL FFADDC (word, INLINE(1))
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL DBVIX('E', 1, IEV)
         IDXFLT = LOCSTR (WORD, NVAREL, NAMES(IEV))
         if (idxflt .eq. 0) then
           CALL PRTERR ('CMDERR',
     &       'Variable name specified does not exist on database.')
           VERB = ' '
           GOTO 150
         end if
         CALL FFADDC (NAMES(idxflt+iev-1), INLINE(1))

C ... Check for comparison -- uses fortran LT, LE, EQ, NE, GT, GE
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         if (.not. matstr(word, 'LT', 2) .and.
     *       .not. matstr(word, 'LE', 2) .and.
     *       .not. matstr(word, 'EQ', 2) .and.
     *       .not. matstr(word, 'NE', 2) .and.
     *       .not. matstr(word, 'GT', 2) .and.
     *       .not. matstr(word, 'GE', 2)) then
           CALL PRTERR ('CMDERR',
     &       'Expected "LT", "LE", "EQ", "NE", "GT", or "GE".')
           VERB = ' '
           GOTO 150
         end if
         CALL FFADDC (word, INLINE(1))
         if (matstr(word, 'LT', 2)) cmpflt = 1
         if (matstr(word, 'LE', 2)) cmpflt = 2
         if (matstr(word, 'EQ', 2)) cmpflt = 3
         if (matstr(word, 'NE', 2)) cmpflt = 4
         if (matstr(word, 'GT', 2)) cmpflt = 5
         if (matstr(word, 'GE', 2)) cmpflt = 6

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &     'filter value', 0.0, VALFLT, *150)
         CALL FFADDR (VALFLT, INLINE(1))

         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         if (.not. matstr(word, 'TIME', 1)) THEN
           CALL PRTERR ('CMDERR',
     &       'Expected "TIME". Syntax:')
           CALL PRTERR ('CMDSPEC',
     *       'FILTER ELEMENTS {elem_var} {CMP} {value} TIME {time}')
           VERB = ' '
           GOTO 150
         ENDIF
         CALL FFADDC (word, INLINE(1))

         CALL FFREAL (IFLD, INTYP, RFIELD,
     &     'database time', 0.0, TIMFLT, *150)
C ... Check that specified time is in range of times on database. Assume TIMES() is sorted.
         if (timflt .lt. times(1)) then
           call prterr('WARNING',
     *       'Specified time is less than minimum database time.')
           call prterr('CMDSPEC',
     *       'Setting filter time to minimum database time.')
           timflt = times(1)
         else if (timflt .gt. times(nsteps)) then
           call prterr('WARNING',
     *       'Specified time is greater than maximum database time.')
           call prterr('CMDSPEC',
     *       'Setting filter time to maximum database time.')
           timflt = times(nsteps)
         end if
         CALL FFADDR (TIMFLT, INLINE(1))
         ISFILTER = .TRUE.

      ELSE IF (VERB .EQ. 'REMOVE') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         if (.not. matstr(word, 'ELEMENTS', 1)) THEN
           CALL PRTERR ('CMDERR',
     &       'Expected "ELEMENTS".  Syntax: FILTER ELEMENTS {elem_var}')
           VERB = ' '
           GOTO 150
         ENDIF
         CALL FFADDC (word, INLINE(1))
C ... See if LOCAL or GLOBAL or nothing specified for id space...
         IF (INTYP(IFLD) .EQ. 0) THEN
           call ffchar(IFLD, INTYP, CFIELD, ' ', WORD)
           if (matstr(word, 'GLOBAL', 1)) then
             idsglobal = .true.
           else if (matstr(word, 'LOCAL', 1)) then
             idsglobal = .false.
           end if
         else
C ... Default is local ids or last set value.
           if (irmcnt .eq. 0) then
             idsglobal = .false.
           end if
         END IF
C ... Gather ids of elements to delete...
 99      continue
         IF (FFEXST (IFLD, INTYP)) THEN
           CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'element ID', 0, ID, *129)
           IRMCNT = IRMCNT + 1
           if (IRMCNT .gt. 1024) then
             CALL PRTERR ('CMDERR',
     &         'Can only specify a maximum of 1024 elements to delete')
             goto 129
           end if
           idsrem(irmcnt) = id
           goto 99
         END IF
 129     continue
         ISREMOVE = .TRUE.

      ELSE IF (VERB .EQ. 'ZOOM') THEN
         CALL FFADDC (VERB, INLINE(1))
         CALL FFNEED (IFLD, INTYP, 'R', MIN(3,NDIM)*2,
     &      'zoom mesh limits for all dimensions', *150)
         DO 100 I = 1, MIN(3,NDIM)*2
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &         'zoom mesh limit', 0.0, ZMLIM(I), *150)
            CALL FFADDR (ZMLIM(I), INLINE(1))
  100    CONTINUE
         CALL FFCHAR (IFLD, INTYP, CFIELD, 'INSIDE', WORD)
         IF (MATSTR (WORD, 'INSIDE', 1)) THEN
           CALL FFADDC ('INSIDE', INLINE(1))
           ZOOMIN = .TRUE.
         ELSE IF (MATSTR (WORD, 'OUTSIDE', 1)) THEN
           CALL FFADDC ('OUTSIDE', INLINE(1))
           ZOOMIN = .FALSE.
         END IF
C        Set logical variables which determine if zoom mesh limits
C        are of the correct format (minlimit < maxlimit)
         XLIM = (ZMLIM(1) .LE. ZMLIM(2))
         YLIM = (ZMLIM(3) .LE. ZMLIM(4))
         IF (NDIM .EQ. 3) THEN
            ZLIM = (ZMLIM(5) .LE. ZMLIM(6))
         ELSE
            ZLIM = .TRUE.
         ENDIF
         IF ((.NOT. XLIM) .OR. (.NOT. YLIM) .OR.
     &       (.NOT. ZLIM)) THEN
C           Zoom limits are incorrect
            CALL PRTERR ('CMDERR',
     &      'Zoom limits are incorrect. (minlimit > maxlimit)')
            INLINE(1) = ' '
         ELSE
            ISZOOM = .TRUE.
         END IF
      ELSE IF (VERB .EQ. 'VISIBLE') THEN
         CALL FFADDC (VERB, INLINE(1))

         IF (.NOT. FFEXST (IFLD, INTYP)) THEN

C         --Select all element blocks if no parameters

            CALL INILOG (NELBLK, .TRUE., VISELB)
            OPTION = '+'

         ELSE IF (.NOT. FFNUMB (IFLD, INTYP)) THEN

C         --Strip off ADD or DELETE option
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ADD', 1)) THEN
               CALL FFADDC ('ADD', INLINE(1))
               OPTION = '+'
            ELSE IF (MATSTR (WORD, 'DELETE', 1)) THEN
               CALL FFADDC ('DELETE', INLINE(1))
               OPTION = '-'
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "ADD" or "DELETE" or element block ID')
               GOTO 160
            END IF

         ELSE

C         --De-select all element blocks so only listed blocks are selected
            CALL INILOG (NELBLK, .FALSE., VISELB)
            OPTION = '+'
         END IF

  110    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'element block ID', 0, ID, *120)
            IELB = LOCINT (ID, NELBLK, IDELB)
            IF (IELB .GT. 0) THEN
               CALL FFADDI (ID, INLINE(1))
               VISELB(IELB) = (OPTION .EQ. '+')
            ELSE
               CALL INTSTR (1, 0, ID, WORD, LSTR)
               CALL PRTERR ('CMDERR', 'Element block ID ' //
     &            WORD(:LSTR) // ' does not exist, ignored')
            END IF
  120       CONTINUE
            GOTO 110
         END IF

      ELSE IF ((VERB .EQ. 'BLOCKS') .OR. (VERB .EQ. 'MATERIAL')) THEN
         CALL FFADDC (VERB, INLINE(1))
         IF (VERB .EQ. 'MATERIAL') VERB = 'BLOCKS'

         IF (.NOT. FFEXST (IFLD, INTYP)) THEN

C         --Select all element blocks if no parameters
            CALL INILOG (NELBLK, .TRUE., SELELB)
            OPTION = '+'

         ELSE IF (.NOT. FFNUMB (IFLD, INTYP)) THEN

C         --Strip off ADD or DELETE option
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'ADD', 1)) THEN
               CALL FFADDC ('ADD', INLINE(1))
               OPTION = '+'
            ELSE IF (MATSTR (WORD, 'DELETE', 1)) THEN
               CALL FFADDC ('DELETE', INLINE(1))
               OPTION = '-'
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "ADD" or "DELETE" or element block ID')
               GOTO 160
            END IF

         ELSE

C         --De-select all element blocks so only listed blocks are selected
            CALL INILOG (NELBLK, .FALSE., SELELB)
            OPTION = '+'
         END IF

  130    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &         'element block ID', 0, ID, *140)
            IELB = LOCINT (ID, NELBLK, IDELB)
            IF (IELB .GT. 0) THEN
               CALL FFADDI (ID, INLINE(1))
               SELELB(IELB) = (OPTION .EQ. '+')
            ELSE
               CALL INTSTR (1, 0, ID, WORD, LSTR)
               CALL PRTERR ('CMDERR', 'Element block ID ' //
     &            WORD(:LSTR) // ' does not exist, ignored')
            END IF
  140       CONTINUE
            GOTO 130
         END IF

         RETVRB = 'BLOCKS'

      ELSE IF ((VERB .EQ. 'LOG') .or. (verb .eq. 'SAVELOG')) THEN
         if (verb .eq. 'SAVELOG') then
            call prterr ('CMDSPEC', 'Please use the LOG command')
            verb = 'LOG'
         end if
         VERB = ' '
         RETVRB = 'LOG'

      ELSE IF (VERB .EQ. 'LIST') THEN
         VERB = ' '
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL DBLIST (WORD, A, NAMECO, BLKTYP, NAMES,
     &                TIMES, IDELB, QAREC, INFREC, MERR)
         IF (MERR .EQ. 1) RETURN

      ELSE IF (VERB .EQ. 'SHOW') THEN
         VERB = ' '
         CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
         CALL SHOW (WORD, CMDTBL, NAMECO, BLKTYP, NAMES,
     &        TIMES, IPTIMS, IDELB, VISELB, SELELB)

         IF (WORD .EQ. ' ') RETVRB = 'PRINT'

      ELSE IF (VERB .EQ. 'HELP') THEN
         VERB = ' '
         CALL ABRSTR (HLPTYP, WORD, HLPTBL)
         ISON = HELP ('ALGEBRA', HLPTYP, ' ')
         IF (.NOT. ISON) THEN
            IF (HLPTYP .EQ. 'RULES') THEN
               CONTINUE
            ELSE IF (HLPTYP .EQ. 'COMMANDS') THEN
               CALL SHOCMD ('COMMANDS', CMDTBL)
            ELSE IF (HLPTYP .EQ. 'FUNCTION') THEN
               CALL SHOCMD ('Available Functions:', FNCNAM)
            ELSE
               CALL SHOCMD ('HELP Options:', HLPTBL)
            END IF
         END IF

         RETVRB = 'PRINT'

      ELSE IF ((VERB .EQ. 'END') .OR. (VERB .EQ. 'EXIT')) THEN
         VERB = ' '
         CALL FFADDC ('END', INLINE(1))
         RETVRB = 'END'

      ELSE IF (VERB .EQ. 'QUIT') THEN
         VERB = ' '
         CALL FFADDC ('QUIT', INLINE(1))
         RETVRB = 'QUIT'

C   --This command allows selectable debugging
      else if (verb .eq. 'DEBUG') then
         verb = ' '
         call ffchar (ifld, intyp, cfield, ' ', cdebug)
         if (matstr (cdebug, 'EQUATION', 3)) then
            call prtdeb ('EQUATION', 0)
         else if (matstr (cdebug, 'VARIABLE', 3)) then
            call prtdeb ('VARIABLE', 0)
         else if (matstr (cdebug, 'ALL', 3)) then
            call prtdeb ('VARIABLE', 0)
            call prtdeb ('EQUATION', -1)
         end if

      ELSE
         CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &      // '" is an invalid command')
         VERB = ' '
         GOTO 150
      END IF

      GOTO 160

  150 CONTINUE
      INLINE(1) = ' '

  160 CONTINUE
      IF (VERB .NE. ' ') THEN
         CALL SHOW (VERB, CMDTBL, NAMECO, BLKTYP, NAMES,
     &      TIMES, IPTIMS, IDELB, VISELB, SELELB)
      END IF

      RETURN
      END
