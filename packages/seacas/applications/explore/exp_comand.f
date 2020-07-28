C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE COMAND (A, IA, EXODUS, DBNAME, QAREC, INFO,
     &     NAMECO, EBTYPE, EBNAME, ATNAME,
     $     NAMIGV, NAMINV, NAMIEV, NAMINS, NAMISS,
     &     NAMOGV, NAMONV, NAMOEV, NAMONS, NAMOSS,
     &     CORD, MAPEL, DBMAPEL, MAPND, DBMAPND, DOMAPN, DOMAPE,
     &     do_check, IDELB, NUMELB, LENE, NUMLNK, NUMATR, LINK, ATRIB,
     &     IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS, LTNNPS, FACNPS, NSNAME,
     &     IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTSESS, FACESS,
     $     SSNAME, ISEVOK, ISNSVOK, ISSSVOK, TIMES,
     $     VARGL, VARNP, VAREL,  VARNS,  VARSS,
     &     LISNP, NLISEL, LISEL, LISNPS, LISESS,
     &     LISGV, LISNV, LISEV, LISMV, LISSV)
C=======================================================================

C   --*** COMAND *** (EXPLORE) Input and process commands
C   --
C   --COMAND inputs and executes an user command.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   EXODUS - IN - true iff database is in the EXODUS database format
C   --   DBNAME - the database name
C   --   QAREC - IN - the QA records containing:
C   --      (1) - the analysis code name
C   --      (2) - the analysis code QA descriptor
C   --      (3) - the analysis date
C   --      (4) - the analysis time
C   --   INFO - IN - the information records
C   --   NAMECO - IN - the names of the coordinates
C   --   EBTYPE - IN - the names of the element block types
C   --   EBNAME - IN - the names of the element blocks
C   --   NAMIGV - IN - the names of the global variables as input
C   --   NAMINV - IN - the names of the nodal variables as input
C   --   NAMIEV - IN - the names of the element variables as input
C   --   NAMOGV - IN - the names of the global variables for comparison
C   --   NAMONV - IN - the names of the nodal variables for comparison
C   --   NAMOEV - IN - the names of the element variables for comparison
C   --   CORD - IN - the coordinates
C   --   MAPEL - IN - the element order map
C   --   MAPND - IN - the node order map
C   --   IDELB - IN - the element block ID for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   LENE - IN - the cumulative element count by element block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the connectivity array for all blocks
C   --   ATRIB - IN - the attribute array for all blocks
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   NDNPS - IN - the number of distribution factors for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   IXDNPS - IN - the index of the first dist factor for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --   NSNAME - IN - the names of the nodesets
C   --   IDESS - IN - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTSESS - IN - the element sides for all sets
C   --   FACESS - IN - the distribution factors for all sets
C   --   SSNAME - IN - the names of the sidesets
C   --   TIMES - IN - the times for all time steps
C   --   VARGL - SCRATCH - the global variables for current time step
C   --   VARNP - SCRATCH - the nodal variables for current time step
C   --   VAREL - SCRATCH - the element variables for current time step
C   --   LISNP - SCRATCH - the indices of the selected nodes
C   --   NLISEL - IN/OUT - the number of selected elements for each block
C   --   LISEL - IN/OUT - the indices of the selected elements (by block)
C   --   LISNPS - SCRATCH - the indices of the selected nodal point sets
C   --   LISESS - SCRATCH - the indices of the selected element side sets
C   --   LISGV - SCRATCH - the indices of the selected global variables
C   --   LISNV - SCRATCH - the indices of the selected nodal variables
C   --   LISEV - SCRATCH - the indices of the selected element variables
C   --
C   --Common Variables:
C   --   Uses NDIM, NUMNP, NUMEL, NELBLK,
C   --      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, LESSDF,
C   --      NQAREC, NINFO, NVARGL, NVARNP, NVAREL, NSTEPS of /DBNUMS/
C   --   Uses NOUT, NCRT, NPRT, ANYPRT of /OUTFIL/

      PARAMETER (MAXFLD = 80)
      include 'exodusII.inc'
      INCLUDE 'exp_progqa.blk'
      INCLUDE 'exp_dbtitl.blk'
      INCLUDE 'exp_dbase.blk'
      INCLUDE 'exp_dbnums.blk'
      INCLUDE 'exp_outfil.blk'
      include 'exp_errcnt.blk'

      DIMENSION A(*)
      INTEGER IA(*)
      LOGICAL EXODUS, FFMATC, DO_CHECK
      LOGICAL DOMAPN, DOMAPE, DOBLK, DOELE
      CHARACTER*(*) DBNAME
      CHARACTER*(MXSTLN) QAREC(4,*)
      CHARACTER*(MXLNLN) INFO(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(NAMLEN) NAMECO(*)
      CHARACTER*(NAMLEN) EBNAME(*), ATNAME(*), NSNAME(*), SSNAME(*)
      CHARACTER*(NAMLEN) NAMIGV(*), NAMINV(*), NAMIEV(*),
     $     NAMINS(*), NAMISS(*)
      CHARACTER*(NAMLEN) NAMOGV(*), NAMONV(*), NAMOEV(*),
     $     NAMONS(*), NAMOSS(*)
      REAL CORD(*)
      INTEGER MAPEL(*)
      INTEGER DBMAPEL(*)
      INTEGER MAPND(*)
      INTEGER DBMAPND(*)
      INTEGER IDELB(*), NUMELB(*)
      INTEGER LENE(0:*)
      INTEGER NUMLNK(*), NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      INTEGER IDNPS(*), NNNPS(*), IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      INTEGER IDESS(*), NEESS(*), NNESS(*), IXEESS(*), IXNESS(*)
      INTEGER LTEESS(*), LTSESS(*)
      REAL FACESS(*)
      INTEGER ISEVOK(*), ISNSVOK(*), ISSSVOK(*)
      REAL TIMES(*)
      REAL VARGL(*), VARNP(*), VAREL(*), VARNS(*), VARSS(*)
      INTEGER LISNP(0:*)
      INTEGER NLISEL(0:*), LISEL(0:*)
      INTEGER LISNPS(0:*), LISESS(0:*)
      INTEGER LISGV(0:*), LISNV(0:*), LISEV(0:*), LISMV(0:*), LISSV(0:*)

      LOGICAL FFEXST, MATSTR

      CHARACTER*(MXNAME) WORD, VERB, LISTYP
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*(MXNAME) CFIELD(MAXFLD)
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
      LOGICAL HELP
      LOGICAL ISON
      CHARACTER*(MXNAME) OPT
      CHARACTER*(MXNAME) MMNAME
      CHARACTER MMTYP
      CHARACTER*80 DUMLIN
      character*2048 OUTPUT, LCOUTPUT

      CHARACTER*(MXSTLN) CMDTBL(15), SELTBL(15), LISTBL(38)
      SAVE CMDTBL, SELTBL, LISTBL, KINVC, KINVS
C      --CMDTBL - the valid commands table
C      --SELTBL - the valid SELECT options table
C      --LISTBL - the valid LIST/PRINT options table

C   --Command table follows.  Remember to change the dimensioned size when
C   --changing the table.
      DATA CMDTBL /
     1  'SELECT  ', 'LIST    ', 'PRINT   ', 'LIMITS  ',
     2  'MINMAX  ', 'CHECK   ', 'HELP    ', 'MAP     ',
     3  'EXIT    ', 'MAXERRS ', 'END     ', 'QUIT    ',
     4  'PRECISION','OUTPUT  ',
     5  '        ' /
      DATA SELTBL /
     1  'NODES   ', 'ELEMENTS', 'BLOCKS  ', 'MATERIAL',
     2  'NSETS   ', 'SSETS   ',
     3  'READ    ', 'STEP    ', 'TIME    ',
     4  'GVARS   ', 'NVARS   ', 'EVARS   ', 'NSVARS  ', 'SSVARS  ',
     5  '        ' /
      DATA LISTBL /
     1  'TITLE   ', 'VARS    ',
     2  'COORDINA', 'MAP     ', 'NMAP    ', 'NODEMAP ',
     3  'BLOCKS  ', 'MATERIAL', 'LINK    ', 'CONNECTI', 'ATTRIBUT',
     4  'NSETS   ', 'NNODES  ', 'NFACTORS', 'INVCON  ',
     5  'SSETS   ', 'SELEMS  ', 'SFACES  ', 'SFACTORS',
     6  'NAMES   ', 'TIMES   ', 'STEPS   ',
     7  'GVARS   ', 'GLOBALS ', 'GLOBNODE', 'GLOBELEM',
     8  'NVARS   ', 'NODALS  ', 'EVARS   ', 'ELEMENTS',
     9  'VERSION ', 'COMMANDS', 'NSVARS  ', 'NODESETVARS',
     $  'SSVARS  ', 'SIDESETVARS', 'FRAMES  ', '        ' /

      DATA KINVC,KINVS /0,0/
C   --Initialize

      OUTPUT = "explore.o"
      MAXERRS = 10

      LENE(0) = 0
      DO 100 I = 1, NELBLK
        LENE(I) = LENE(I-1) + NUMELB(I)
 100  CONTINUE

      LISNP(0) = NUMNP
      DO 110 I = 1, NUMNP
        LISNP(I) = I
 110  CONTINUE
      NLISEL(0) = NELBLK
      DO 120 I = 1, NELBLK
        NLISEL(I) = LENE(I) - LENE(I-1)
 120  CONTINUE
      LISEL(0) = NUMEL
      DO 130 I = 1, NUMEL
        LISEL(I) = I
 130  CONTINUE
      LISNPS(0) = NUMNPS
      DO 140 I = 1, NUMNPS
        LISNPS(I) = I
 140  CONTINUE
      LISESS(0) = NUMESS
      DO 150 I = 1, NUMESS
        LISESS(I) = I
 150  CONTINUE
      IF (EXODUS) THEN
        LISGV(0) = NVARGL
        DO 170 I = 1, NVARGL
          LISGV(I) = I
 170    CONTINUE
        LISNV(0) = NVARNP
        DO 180 I = 1, NVARNP
          LISNV(I) = I
 180    CONTINUE
        LISEV(0) = NVAREL
        DO 190 I = 1, NVAREL
          LISEV(I) = I
 190    CONTINUE
        LISMV(0) = NVARNS
        DO 191 I = 1, NVARNS
          LISMV(I) = I
 191   CONTINUE
        LISSV(0) = NVARSS
        DO 192 I = 1, NVARSS
          LISSV(I) = I
 192   CONTINUE
      END IF

      mmname = ' '
      mmvar  = 0

C   --Read first time step variables

      IF (EXODUS) THEN
        NCSTEP = 999
        NSTEPNS = -1
        NSTEPSS = -1
        NSTEP = 1
        CALL TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &    TIME, VARGL, VARNP, VAREL)
      END IF

      WRITE (*, *)
      CALL PRTERR ('CMDREQ',
     & 'Use "precision low|normal|high|#" to control" output precision')

      if (domape .and. domapn) then
        call PRTERR('CMDREQ',
     *    'Nodes and Elements using Global Ids')
      else if (domape) then
        call PRTERR('CMDREQ',
     *    'Elements use Global Ids, Node Ids are Local')
      else if (domapn) then
        call PRTERR('CMDREQ',
     *    'Element use Local Ids, Node Ids are Global')
      else
        call PRTERR('CMDREQ',
     *    'Nodes and Elements using Local Ids')
      end if

 200  CONTINUE

C   --Read command line

      if (do_check) then
        call check(a, ia, exodus, idelb, ebtype, numelb, isevok, numlnk,
     *    numatr, link, atrib, atname, mapnd, dbmapnd, mapel, dbmapel,
     *    idnps, nnnps, ixnnps,
     *    ltnnps, facnps, idess, neess, nness, ixeess, ixness, lteess,
     *    ltsess, facess, vargl, varnp, varel)
        return
      endif
      WRITE (*, *)
      CALL FREFLD (0, 0, 'EXPLORE> ', MAXFLD,
     &  IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      IF (IOSTAT .LT. 0) GOTO 280
      IF (NUMFLD .EQ. 0) GOTO 200
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      IFLD = 1
      CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
      CALL ABRSTR (VERB, WORD, CMDTBL)
      IF ((VERB .EQ. ' ') .AND. (WORD .NE. ' ')) THEN
        VERB = '*'
        IFLD = IFLD - 1
      END IF

      IF (VERB .EQ. 'SELECT') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL ABRSTR (LISTYP, WORD, SELTBL)
        IF (LISTYP .EQ. ' ') LISTYP = WORD

      ELSE IF ((VERB .EQ. 'LIST') .OR. (VERB .EQ. 'PRINT')) THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        IF (VERB .EQ. 'PRINT') THEN
          IF ((WORD .EQ. ' ')
     &      .OR. (WORD .EQ. 'ON') .OR. (WORD .EQ. 'OFF')) THEN
            CALL PRTERR ('CMDERR',
     &        'Please use the new PRINT command'
     &        // ' (e.g., PRINT NAMES)')
            GOTO 270
          END IF
        END IF
        CALL ABRSTR (LISTYP, WORD, LISTBL)
        IF (LISTYP .EQ. ' ') LISTYP = WORD

      ELSE IF (VERB .EQ. '*') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
        CALL ABRSTR (LISTYP, WORD, LISTBL)
        IF (LISTYP .EQ. ' ') THEN
          CALL ABRSTR (LISTYP, WORD, SELTBL)
          IF (LISTYP .NE. ' ') THEN
            CALL PRTERR ('CMDREQ', 'Please use the SELECT command')
            VERB = 'SELECT'
          ELSE
            LISTYP = WORD
          END IF
        END IF
      END IF

      IF (VERB .EQ. 'PRINT') THEN
        IF (NPRT .LE. 0) THEN
          CALL PRTERR ('CMDERR', 'Print file cannot be opened')
          GOTO 270
        END IF

        IF (.NOT. ANYPRT) THEN

C         --Set up the print file

          OPEN (unit=nprt, file=OUTPUT, IOSTAT=IERR)
          IF (IERR .NE. 0) THEN
            CALL PRTERR ('CMDERR', 'Print file cannot be opened')
            NPRT = -1
            GOTO 270
          END IF

          CALL BANNER (NPRT, QAINFO,
     &      ' ', ' ', ' ')
          CALL PRINIT ('N', NPRT, DBNAME, TITLE,
     &      NDIM, NUMNP, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF,
     &      NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)

          ANYPRT = .TRUE.
        END IF

        NOUT = NPRT
      ELSE
        NOUT = NCRT
      END IF

      DUMLIN = ' '

C *** GENESIS Print Commands ***

      IF ((VERB .EQ. 'SELECT') .OR.
     &  (FFEXST (IFLD, INTYP) .AND. ((VERB .EQ. '*')
     &  .OR. (VERB .EQ. 'LIST') .OR. (VERB .EQ. 'PRINT')))) THEN

        IF (LISTYP .EQ. 'NODES') THEN
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)

          if (FFMATC (IFLD, INTYP, CFIELD, 'NSET', 4) .OR.
     *      FFMATC (IFLD, INTYP, CFIELD, 'NODESET', 7)) THEN
            CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &        'nodal point set ID',
     &        NUMNPS, IDNPS, LISNPS(0), LISNPS(1), *270)
            if (lisnps(0) .gt. 0) then
              call selset(lisnp(0), lisnp(1),
     *          numnps, lisnps, lnpsnl,
     *          idnps, nnnps, ixnnps, ltnnps, "nodes")
            end if

          else if (FFMATC (IFLD, INTYP, CFIELD, 'BLOCK', 3) .OR.
     *        FFMATC (IFLD, INTYP, CFIELD, 'MATERIAL', 3)) THEN
            CALL MDRSRV ('SCRSEL', KLELB, 1+NELBLK)
            CALL MDRSRV ('SCR',    KSCR,  NUMNP)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 280

            CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &        'element block ID',
     &        NELBLK, IDELB, IA(KLELB), IA(KLELB+1), *205)
 205        continue
            if (IA(KLELB) .gt. 0) then
              call selblk(lisnp(0), lisnp(1),
     *          nelblk, IA(KLELB), numelb, numlnk, link,
     *          A(KSCR), numnp, ebtype)
            end if
            CALL MDDEL ('SCRSEL')
            CALL MDDEL ('SCR')
          else
            CALL RMIXINT (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &        'node number', NUMNP, LISNP(0), LISNP(1), MAPND, *270)
          end if

        ELSE IF ((VERB .NE. 'SELECT')
     &      .AND. (LISTYP .EQ. 'COORDINA')) THEN
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)
          CALL CKNONE (NDIM, .FALSE., 'coordinates', *270)

          CALL RMIXINT (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &      'node number', NUMNP, LISNP(0), LISNP(1), MAPND, *270)

          ELSE IF ((LISTYP .EQ. 'ELEMENTS')
     &      .OR. ((VERB .NE. 'SELECT')
     &      .AND. ((LISTYP .EQ. 'LINK') .OR. (LISTYP .EQ. 'CONNECTI')
     &      .OR. (LISTYP .EQ. 'ATTRIBUT')))) THEN
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKNONE (NUMEL, .FALSE., 'elements', *270)

          CALL MDRSRV ('SCRSEL', KLEL, 1+NUMEL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 280

          if (FFMATC (IFLD, INTYP, CFIELD, 'SSET', 4) .OR.
     *      FFMATC (IFLD, INTYP, CFIELD, 'SIDESET', 7)) THEN
            CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &        'side set ID',
     &        NUMESS, IDESS, LISESS(0), LISESS(1), *270)
            if (lisess(0) .gt. 0) then
              call selset(IA(KLEL), IA(KLEL+1),
     *          numess, lisess, lessel,
     *          idess, neess, ixeess, lteess, "elements")
            end if
          else
            CALL RMIXINT (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &        'element number', NUMEL, IA(KLEL), IA(KLEL+1), MAPEL,
     *        *270)
          end if

          CALL DBSBEL (NELBLK, NUMEL, LENE, A(KLEL), NLISEL, LISEL)

          CALL MDDEL ('SCRSEL')

        ELSE IF ((LISTYP .EQ. 'BLOCKS')
     &      .OR. (LISTYP .EQ. 'MATERIAL')) THEN
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKNONE (NELBLK, .FALSE., 'element blocks', *270)

          CALL MDRSRV ('SCRSEL', KLELB, 1+NELBLK)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 280

          CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &      'element block ID',
     &      NELBLK, IDELB, IA(KLELB), IA(KLELB+1), *220)
 220      CONTINUE

          CALL DBSELB (NELBLK, NUMEL, LENE, A(KLELB),
     &      NLISEL, LISEL)

          CALL MDDEL ('SCRSEL')

        ELSE IF (LISTYP .EQ. 'NSETS') THEN
          NSTEPNS = -1
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKNONE (NUMNPS, .FALSE., 'nodal point sets', *270)

          CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &      'nodal point set ID',
     &      NUMNPS, IDNPS, LISNPS(0), LISNPS(1), *270)

        ELSE IF (LISTYP .EQ. 'SSETS') THEN
          NSTEPSS = -1
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKNONE (NUMESS, .FALSE., 'element side sets', *270)

          CALL RIXID (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &      'element side set ID',
     &      NUMESS, IDESS, LISESS(0), LISESS(1), *270)

C *** EXODUS Movement Commands ***

        ELSE IF (((VERB .EQ. 'SELECT') .OR. (VERB .EQ. '*'))
     &      .AND. ((LISTYP .EQ. 'READ') .OR. (LISTYP .EQ. 'STEP')
     &      .OR. (LISTYP .EQ. 'TIME'))) THEN
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NSTEPS, .FALSE., 'time steps', *270)

          IF (LISTYP .EQ. 'READ') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &        'number of steps to read', 1, N, *230)
            NSTEP = MAX (1, NCSTEP + N)
          ELSE IF (LISTYP .EQ. 'STEP') THEN
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &        'step number', NCSTEP, N, *230)
            NSTEP = N
          ELSE IF (LISTYP .EQ. 'TIME') THEN
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'time step time', TIME, T, *230)
            NSTEP = LOCREA (T, NSTEPS, TIMES)
          END IF

          CALL TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK, TIME,
     &      VARGL, VARNP, VAREL)

 230      CONTINUE
          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

C *** EXODUS Print Commands ***

        ELSE IF (LISTYP .EQ. 'GVARS') THEN
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVARGL, .FALSE., 'global variables', *270)

          CALL RIXWRD (DUMLIN, IFLD, INTYP, CFIELD,
     &      'global variable name', NVARGL, NAMOGV,
     &      LISGV(0), LISGV(1), *270)

        ELSE IF (LISTYP .EQ. 'NVARS') THEN
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVARNP, .FALSE., 'nodal variables', *270)

          CALL RIXWRD (DUMLIN, IFLD, INTYP, CFIELD,
     &      'nodal variable name', NVARNP, NAMONV,
     &      LISNV(0), LISNV(1), *270)

        ELSE IF ((VERB .NE. 'SELECT')
     &      .AND. (LISTYP .EQ. 'NODALS')) THEN
          IF (VERB .EQ. '*') VERB = 'LIST'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)

          CALL RMIXINT (DUMLIN, IFLD, INTYP, CFIELD, IFIELD,
     &      'node number', NUMNP, LISNP(0), LISNP(1), MAPND, *270)

        ELSE IF (LISTYP .EQ. 'EVARS') THEN
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVAREL, .FALSE., 'element variables', *270)

          CALL RIXWRD (DUMLIN, IFLD, INTYP, CFIELD,
     &      'element variable name', NVAREL, NAMOEV,
     &      LISEV(0), LISEV(1), *270)

        ELSE IF (LISTYP .EQ. 'NSVARS') THEN
          NSTEPNS = -1
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVARNS, .FALSE., 'nodeset variables', *270)

          CALL RIXWRD (DUMLIN, IFLD, INTYP, CFIELD,
     &      'nodeset variable name', NVARNS, NAMONS,
     &      LISMV(0), LISMV(1), *270)

        ELSE IF (LISTYP .EQ. 'SSVARS') THEN
          NSTEPSS = -1
          IF (VERB .EQ. '*') VERB = 'SELECT'
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVARSS, .FALSE., 'sideset variables', *270)

          CALL RIXWRD (DUMLIN, IFLD, INTYP, CFIELD,
     &      'sideset variable name', NVARSS, NAMOSS,
     &      LISSV(0), LISSV(1), *270)

        ELSE IF (VERB .EQ. 'SELECT') THEN
          CALL SHOCMD ('SELECT Options:', SELTBL)
        END IF
      END IF

      IF (VERB .EQ. '*') VERB = 'LIST'

      IF (VERB .EQ. 'SELECT') THEN
        CONTINUE

      ELSE IF ((VERB .EQ. 'LIST') .OR. (VERB .EQ. 'PRINT')) THEN

        IF (listyp .eq. 'GLOBNODE') THEN
 233      continue
          if (ffexst(ifld, intyp)) THEN
            call ffintg(ifld, intyp, ifield, 'Global ID',
     *        0, IDGLO, *270)
            call selmap(nout, 'node', idglo, numnp, mapnd)
            go to 233
          end if

        ELSE IF (listyp .eq. 'GLOBELEM') THEN
 234      continue
          if (ffexst(ifld, intyp)) then
            call ffintg(ifld, intyp, ifield, 'Global ID',
     *        0, IDGLO, *270)
            call selmap(nout, 'element', idglo, numel, mapel)
            goto 234
          end if

        ELSE IF ((LISTYP .EQ. 'TITLE') .OR. (LISTYP .EQ. 'VARS')) THEN
          IF (EXODUS) THEN
            CALL PRINIT ('NTSICV', NOUT, DBNAME, TITLE,
     &        NDIM, NUMNP, NUMEL, NELBLK,
     &        NUMNPS, LNPSNL, LNPSDF,
     &        NUMESS, LESSEL, LESSNL, LESSDF,
     &        NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
          ELSE
            CALL PRINIT ('NTSIC', NOUT, DBNAME, TITLE,
     &        NDIM, NUMNP, NUMEL, NELBLK,
     &        NUMNPS, LNPSNL, LNPSDF,
     &        NUMESS, LESSEL, LESSNL, LESSDF,
     &        NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
          END IF

        ELSE IF (LISTYP .EQ. 'VERSION') THEN
          CALL PRVERS (NDB, NOUT)

        ELSE IF (LISTYP .EQ. 'COORDINA') THEN
          CALL CKNONE (NDIM, .FALSE., 'coordinates', *270)
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)
          CALL CKNONE (LISNP(0), .TRUE., 'nodes', *270)

          CALL PRXYZ ('*', NOUT, NDIM, NAMECO, NUMNP, LISNP, CORD,
     *      MAPND, DOMAPN)

        ELSE IF (LISTYP .EQ. 'MAP') THEN
          CALL CKNONE (NUMEL, .FALSE., 'elements', *270)

          CALL PRMAP ('*', NOUT, 'Element', NUMEL, DBMAPEL)

        ELSE IF (LISTYP .EQ. 'NMAP' .OR. LISTYP .EQ. 'NODEMAP') THEN
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)

          CALL PRMAP ('*', NOUT, 'Node', NUMNP, DBMAPND)

        ELSE IF ((LISTYP .EQ. 'BLOCKS') .OR. (LISTYP .EQ. 'MATERIAL')
     &      .OR. (LISTYP .EQ. 'LINK') .OR. (LISTYP .EQ. 'CONNECTI')
     &      .OR. (LISTYP .EQ. 'ATTRIBUT')) THEN
          CALL CKNONE (NELBLK, .FALSE., 'element blocks', *270)
          CALL CKNONE (NLISEL(0), .TRUE., 'element blocks', *270)

          IF ((LISTYP .EQ. 'BLOCKS')
     &      .OR. (LISTYP .EQ. 'MATERIAL')) THEN
            IF (EXODUS) THEN
              OPT = 'NV'
            ELSE
              OPT = 'N'
            END IF
          ELSE IF ((LISTYP .EQ. 'LINK')
     &        .OR. (LISTYP .EQ. 'CONNECTI')) THEN
            OPT = 'C'
          ELSE IF (LISTYP .EQ. 'ATTRIBUT') THEN
            OPT = 'A'
          END IF

          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDRSRV ('XLISEV', KXLSEV, NVAREL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 280
          END IF

          CALL PRELB (OPT, NOUT, NELBLK, NLISEL, LISEL,
     &      IDELB, LENE, NUMLNK, NUMATR, LINK, ATRIB,
     &      EBTYPE, EBNAME, NVAREL, NAMIEV, ISEVOK, A(KXLSEV),
     *      ATNAME, MAPND, DOMAPN, MAPEL, DOMAPE)

          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDDEL ('XLISEV')
          END IF

        ELSE IF ((LISTYP .EQ. 'NSETS')
     &      .OR. (LISTYP .EQ. 'NNODES')
     &      .OR. (LISTYP .EQ. 'NFACTORS')) THEN
          CALL CKNONE (NUMNPS, .FALSE., 'nodal point sets', *270)
          CALL CKNONE (LISNPS(0), .TRUE., 'nodal point sets', *270)

          IF (LISTYP .EQ. 'NSETS') THEN
             if (exodus) then
                opt = 'V'
             else
                OPT = ' '
             end if
          ELSE IF (LISTYP .EQ. 'NNODES') THEN
            OPT = 'N'
          ELSE IF (LISTYP .EQ. 'NFACTORS') THEN
            OPT = 'F'
          END IF
          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDRSRV ('XLISNV', KXLSNV, NVARNS)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 280
          END IF
          CALL PRNPS (OPT, NOUT, NUMNPS, LISNPS, LNPSNL,
     &         IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS, LTNNPS, FACNPS,
     *         NSNAME, nvarns, namins, isnsvok, a(kxlsnv),
     $         MAPND, DOMAPN)
          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDDEL ('XLISNV')
          END IF

        ELSE IF ((LISTYP .EQ. 'SSETS')
     &      .OR. (LISTYP .EQ. 'SELEMS')
     &      .OR. (LISTYP .EQ. 'SFACES')
     &      .OR. (LISTYP .EQ. 'SFACTORS')) THEN
          CALL CKNONE (NUMESS, .FALSE., 'element side sets', *270)
          CALL CKNONE (LISESS(0), .TRUE., 'element side sets', *270)

          IF (LISTYP .EQ. 'SSETS') THEN
             if (exodus) then
                opt = 'V'
             else
                OPT = ' '
             end if
          ELSE IF (LISTYP .EQ. 'SELEMS') THEN
            OPT = 'E'
          ELSE IF (LISTYP .EQ. 'SFACES') THEN
            OPT = 'N'
          ELSE IF (LISTYP .EQ. 'SFACTORS') THEN
            OPT = 'F'
          END IF
          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDRSRV ('XLISSV', KXLSSV, NVARSS)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 280
          END IF
          IF (INDEX (OPT, 'F') .GT. 0) then
            CALL MDRSRV ('NDFSID', KNDFSID, LESSEL)
            CALL MDRSRV ('NODSID', KNODSID, LESSNL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 280
          END IF

          CALL PRESS (OPT, NOUT, NUMESS, LISESS, LESSEL, LESSNL,
     &         IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &         LTEESS, LTSESS, FACESS, SSNAME,
     $         nvarss, namiss, isssvok, a(kxlssv),
     &         a(kndfsid), a(knodsid),MAPEL, MAPND, DOMAPE, DOMAPN)
          IF (INDEX (OPT, 'V') .GT. 0) THEN
            CALL MDDEL ('XLISSV')
          END IF
          IF (INDEX (OPT, 'F') .GT. 0) THEN
            CALL MDDEL ('NDFSID')
            CALL MDDEL ('NODSID')
          END IF

        ELSE IF (LISTYP .EQ. 'INVCON') THEN
           CALL FFCHAR (IFLD, INTYP, CFIELD,' ', WORD)
           doblk = .false.
           doele = .false.
           IF (MATSTR(WORD, 'ELEMENTS', 1)) THEN
              DOELE = .true.
           ELSE IF (MATSTR(WORD, 'BLOCKS', 1)) THEN
              DOBLK = .true.
           ELSE
              DOELE = .true.
              DOBLK = .true.
           END IF
           if (kinvc .eq. 0) then
              CALL MDRSRV ('INVCON', KINVC, NELBLK*NUMNP)
              CALL MDRSRV ('INVSCR', KINVS, 2*NELBLK)
              CALL MDRSRV ('NODMAP', KNDMP, 2*NUMNP)
              CALL MDSTAT (NERR, MEM)
              IF (NERR .GT. 0) GOTO 280
           end if

           call invcon(ia, nelblk, idelb, numelb, numlnk, link, numnp,
     *          ia(kinvc), ia(kinvs), ia(kndmp), lisnp, NOUT, MAPND,
     $          MAPEL, DOMAPN, DOMAPE, DOBLK, DOELE, ebtype)

        ELSE IF (LISTYP .EQ. 'QA') THEN

          CALL PRQA ('*', NOUT, NQAREC, QAREC, NINFO, INFO)

C *** EXODUS Print Commands ***

        ELSE IF (LISTYP .EQ. 'NAMES') THEN
          CALL CKEXOD (EXODUS, *270)

          CALL PRNAME ('*', NOUT,
     &         NAMIGV, NAMINV, NAMIEV, NAMINS, NAMISS)

        ELSE IF ((LISTYP .EQ. 'TIMES')
     &      .OR. (LISTYP .EQ. 'ALLTIMES')) THEN
          IF (LISTYP .EQ. 'ALLTIMES') CALL PRTERR ('CMDREQ',
     &      'Please use the LIST TIMES command')
          CALL CKEXOD (EXODUS, *270)

          CALL PRTIMS ('NT', NOUT, NSTEPS, TIMES)

        ELSE IF (LISTYP .EQ. 'STEPS') THEN
          CALL CKEXOD (EXODUS, *270)

          CALL PRTIMS ('NM', NOUT, NSTEPS, TIMES)

        ELSE IF ((LISTYP .EQ. 'GVARS')
     &      .OR. (LISTYP .EQ. 'GLOBALS')) THEN
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NVARGL, .FALSE., 'global variables', *270)
          CALL CKNONE (LISGV(0), .TRUE., 'global variables', *270)

          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

          CALL PRGLOB ('*', NOUT, NVARGL, LISGV, NAMIGV, VARGL)

        ELSE IF ((LISTYP .EQ. 'NVARS')
     &      .OR. (LISTYP .EQ. 'NODALS')) THEN
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NUMNP, .FALSE., 'nodes', *270)
          CALL CKNONE (NVARNP, .FALSE., 'nodal variables', *270)
          CALL CKNONE (LISNV(0), .TRUE., 'nodal variables', *270)

          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

          CALL PRNODE ('*', NOUT, NUMNP, LISNP, NVARNP, LISNV, NAMINV,
     &      VARNP, MAPND, DOMAPN)

        ELSE IF ((LISTYP .EQ. 'EVARS')
     &      .OR. (LISTYP .EQ. 'ELEMENTS')) THEN
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NUMEL, .FALSE., 'elements', *270)
          CALL CKNONE (NVAREL, .FALSE., 'element variables', *270)
          CALL CKNONE (LISEV(0), .TRUE., 'element variables', *270)

          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

          CALL PRELEM ('*', NOUT, NELBLK, NUMEL, NLISEL, LISEL, LENE,
     &      NVAREL, LISEV, NAMIEV, ISEVOK, VAREL, max(1,nvarel),
     *      MAPEL, DOMAPE)

        ELSE IF ((LISTYP .EQ. 'NSVARS')
     &      .OR. (LISTYP .EQ. 'NODESETVARS')) THEN
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NUMNPS, .FALSE., 'nodesets', *270)
          CALL CKNONE (NVARNS, .FALSE., 'nodeset variables', *270)
          CALL CKNONE (LISMV(0), .TRUE., 'nodeset variables', *270)

          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

          CALL PRNSV (NOUT, NCSTEP, NUMNPS, LISNPS, LNPSNL,
     &         IDNPS, NNNPS, IXNNPS, LTNNPS, NSNAME,
     $         NVARNS, LISMV(0), NAMINS, ISNSVOK, VARNS, max(1,nvarns),
     *         MAPND, DOMAPN)

        ELSE IF ((LISTYP .EQ. 'SSVARS')
     &      .OR. (LISTYP .EQ. 'SIDESETVARS')) THEN
          CALL CKEXOD (EXODUS, *270)
          CALL CKNONE (NUMESS, .FALSE., 'sidesets', *270)
          CALL CKNONE (NVARSS, .FALSE., 'sideset variables', *270)
          CALL CKNONE (LISSV(0), .TRUE., 'sideset variables', *270)

          CALL PRSTEP ('*', NOUT, TIME, NCSTEP, NSTEPS)

          CALL PRSSV (NOUT, NCSTEP, NUMESS, LISESS, LESSEL,
     &         IDESS, NEESS, IXEESS, LTEESS, LTSESS, SSNAME,
     $         NVARSS, LISSV(0), NAMISS, ISSSVOK, VARSS, max(1,nvarss),
     *         MAPEL, DOMAPE)

C ... Coordinate Frames
        ELSE IF (LISTYP .EQ. 'FRAMES') THEN
          CALL PRFRM (NOUT)
        ELSE IF (LISTYP .EQ. 'COMMANDS') THEN
          CALL SHOCMD ('COMMAND Options:', CMDTBL)
        ELSE
          CALL SHOCMD ('LIST/PRINT Options:', LISTBL)
        END IF

C *** Miscellaneous Commands ***

      ELSE IF (VERB .EQ. 'LIMITS') THEN
        call limits(ndim, numnp, cord)

      ELSE IF (VERB .EQ. 'MAP') THEN
        call prterr('CMDERR',
     *    'The map nodes|elements|both command is no longer'
     *    // ' supported.  The mapping of nodes/elements is'
     *    // ' selected at program startup with the -map '
     *    // ' or -nomap option.')
c$$$        CALL FFCHAR (IFLD, INTYP, CFIELD,'BOTH', WORD)
c$$$        CALL FFONOF (IFLD, INTYP, CFIELD, MAPTMP, *270)
c$$$
c$$$        IF (MATSTR(WORD, 'ELEMENTS', 1)) THEN
c$$$          DOMAPE = MAPTMP
c$$$        ELSE IF (MATSTR(WORD, 'NODES', 1)) THEN
c$$$          DOMAPN = MAPTMP
c$$$        ELSE IF (MATSTR(WORD, 'BOTH', 1)) THEN
c$$$          DOMAPE = MAPTMP
c$$$          DOMAPN = MAPTMP
c$$$        END IF

      ELSE IF (VERB .EQ. 'MAXERRS') THEN
         CALL FFINTG (IFLD, INTYP, IFIELD,
     &        'maximum errors to print; 0 for all', 10, MAXERRS, *235)
 235     continue

      ELSE IF (VERB .EQ. 'OUTPUT') THEN
        CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', OUTPUT)
C ... Filename will be converted to lowercase -- FREFLD returns
C     everything in uppercase Assume user wanted lowercase; if they
C     didn't, need to rewrite frefld to return mixed case.
        call lowstr(lcoutput,output)
        if (anyprt) then
          close (nprt)
        end if
        OPEN (unit=nprt, file=lcoutput, IOSTAT=IERR)
        IF (IERR .NE. 0) THEN
          CALL PRTERR ('CMDERR', 'Print file cannot be opened')
          NPRT = -1
          GOTO 270
        END IF
        ANYPRT = .TRUE.

      ELSE IF (VERB .EQ. 'CHECK') THEN

        call check(a, ia, exodus, idelb, ebtype, numelb, isevok, numlnk,
     *    numatr, link, atrib, atname, mapnd, dbmapnd, mapel, dbmapel,
     *    idnps, nnnps, ixnnps,
     *    ltnnps, facnps, idess, neess, nness, ixeess, ixness, lteess,
     *    ltsess, facess, vargl, varnp, varel)

      ELSE IF (VERB .EQ. 'MINMAX') THEN

        CALL RDMMAX (IFLD, INTYP, CFIELD,
     &    NAMIGV, NAMINV, NAMIEV,
     &    NAMOGV, NAMONV, NAMOEV,
     &    NCSTEP, MMSTEP, MMNAME, MMTYP, MMVAR, MMNUM, *270)

        IF (MMNUM .LT. 0) THEN
          CALL PRTERR ('CMDERR', 'All minimums and maximums'
     &      // ' for this variable have been displayed')
          GOTO 270
        END IF

        CALL PRMMAX (NOUT, MMSTEP, MMNAME, MMTYP, MMVAR, MMNUM,
     &    XMIN, XMAX, NUMELB, IDELB, ISEVOK, VARGL, VARNP, VAREL)

        IF (MMSTEP .EQ. 0) THEN
          NSTEP = 1
          CALL TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &      TIME, VARGL, VARNP, VAREL)
        END IF

C *** Miscellaneous Commands ***

      ELSE IF (VERB .EQ. 'HELP') THEN
        ISON = HELP ('EXPLORE', 'COMMANDS', CFIELD(IFLD))
        IF (.NOT. ISON)
     &    CALL SHOCMD ('COMMANDS', CMDTBL)

      ELSE IF (VERB .EQ. 'PRECISION') THEN
         if (FFMATC (IFLD, INTYP, CFIELD, 'HIGH', 1)) THEN
            IPREC = 9
         else if (FFMATC (IFLD, INTYP, CFIELD, 'LOW', 1) .OR.
     $           FFMATC (IFLD, INTYP, CFIELD, 'NORMAL', 1)) THEN
            IPREC = 4
         else if (ffexst(ifld, intyp)) THEN
            call ffintg(ifld, intyp, ifield, 'Precision',
     *        4, IPREC, *270)
         else
            CALL PRTERR ('CMDERR',
     &        'Syntax: PRECISION high|low|normal|0..9')
            GOTO 270
         end if
         if (IPREC .gt. 16) then
            CALL PRTERR ('CMDERR',
     &        'Maximum precision is 16.')
            IPREC = 16
         end if
         if (IPREC .le. 0) then
            CALL PRTERR ('CMDERR',
     &        'Precision must be positive.')
            GOTO 270
         end if
         CALL SETPRC(IPREC,1)

      ELSE IF ((VERB .EQ. 'EXIT') .OR. (VERB .EQ. 'END')
     &    .OR. (VERB .EQ. 'QUIT')) THEN
        CALL SCNEOF
        GOTO 280

      ELSE
        CALL PRTERR ('CMDERR', '"' // VERB(:LENSTR(VERB))
     &    // '" is an invalid command')
      END IF

 270  CONTINUE
      GOTO 200

 280  CONTINUE
      RETURN
      END
