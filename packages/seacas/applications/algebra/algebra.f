C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      PROGRAM ALGEBRA2
C=======================================================================

C   --This version of ALGEBRA will read and write EXODUSIIV2 database
C   --format files.  Many changes have occurred since the first version
C   --of ALGEBRA.  The original database files, genesis and exodusI
C   --were sequential access files.  EXODUSIIV2 uses a random access
C   --file format. Previous versions of ALGEBRA would have to read the
C   --input database more than once in order to get the file pointer
C   --to the desired data.  With random access files we are able to
C   --select what we want to read or write at anytime.

C                         *** ALGEBRA 2.02 ***
C   --*** ALGEBRA *** (ALGEBRA) Algebraic Database Manipulation Program
C   --
C   --The ALGEBRA program allows the user to process data from a finite
C   --element analysis before it is plotted. The finite element output
C   --data is in the form of variable values (stress, strain, and
C   --velocity components, etc.) in an EXODUS database.  The ALGEBRA program
C   --evaluates user-supplied functions of the data and writes the results
C   --to an output EXODUS database which can be read by plot programs.
C   --
      include 'exodusII.inc'
      include 'ag_namlen.blk'

C     max_num_equation, max_parsed_entries
      include 'ag_numeqn.blk'
C     QAINFO array with program information
      include 'ag_progqa.blk'
C     input equation info
      include 'ag_ent.blk'
C     input variable info
      include 'ag_var.blk'
C     aliases
      include 'ag_alias.blk'
C     I/O file names
      include 'ag_dbase.blk'
C     I/O database titles
      include 'ag_dbtitl.blk'
C     num_of: nodes,elements,coordinated/node,element_blks
C     num_of: history,global,nodal,element variables, num_of_time_steps
      include 'ag_dbnums.blk'
C     node set/side sets number,length,dist factor
      include 'ag_dbnumg.blk'
C     database type, num_of qa and info records
      include 'ag_dbnumq.blk'
C     num_of: nodes,elements,element blks,node_sets,side_sets in zoom mesh
C     num_of: output history,global,nodal,element variables
      include 'ag_dbout.blk'
C     time_index, I_index/O_index for history,global,nodal and element vars
      include 'ag_dbxvar.blk'
C     time variables
      include 'ag_times.blk'
C     zoom info
      include 'ag_zoom.blk'
C     function variables
      include 'ag_fnctbc.blk'
C     equation line error messages
      include 'ag_eqnlns.blk'
C     floating point byte size
      include 'ag_dbws.blk'
      include 'argparse.inc'

C     Input/Output File Arguments
C     CPUWS - The word size in bytes of the floating point variables
C             used in the application program
C     IOWS  - The word size in bytes of the floating point data as they
C             are stored in the EXODUS II file
C     IERR  - error code
C     VERS  - version number of database
C     NERR  - error flag for dynamic memory errors
C     MERR  - memory error flag
C     IOERR - Input/output error flag
C      INTEGER CPUWS, IOWS, IERR, NERR, CERR, MERR
      INTEGER IERR, NERR, CERR, MERR
      REAL VERS
      LOGICAL MEMBUG
      CHARACTER*2048 FILNAM, SCRATCH

      CHARACTER*8 STR8

C     A - the dynamic memory base array
      DIMENSION A(1)
      CHARACTER*1 C(1)

C     Logical variables that indicate whether an input/output
C     file is open.  These variable are use during exit of the
C     program in hopes to cut down on the multitude of goto labels
C     while exiting the program
      LOGICAL INOPEN, OTOPEN

C     Executable code in qainfo.blk
      include 'ag_qainfo.blk'

      INOPEN = .FALSE.
      OTOPEN = .FALSE.
      MEMBUG = .FALSE.
      MERR   = 0
      IOERR  = 0

      CALL STRTUP (QAINFO)
      CALL BANNR2 (80, QAINFO(1), 0)
      CALL BANNER (0, QAINFO,
     &   'AN ALGEBRAIC MANIPULATION PROGRAM',
     &   'FOR POST-PROCESSING OF FINITE ELEMENT ANALYSES',
     &   'EXODUS II VERSION')
      CALL CPYRGT (0, '2008')
      WRITE(*,*)

C     Set compute and file float size to default values
C     values located in dbws.blk
      CPUWS=0
      IOWS=0

      IF (MEMBUG) THEN
         CALL MLIST()
      END IF

C     Initialize the Memory Manager
      CALL MDINIT (A)
      CALL MDFILL(0)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
C        Report dynamic memory error
         CALL MEMERR
C        Goto WRAPUP call and exit program
         GOTO 130
      END IF

C     Input and Output File ID's
C     See Algebra script for more details
C     dbase.blk: COMMON /DBASE/ NDBIN, NDBOUT
      NDBIN = 11
      NDBOUT = 12

C     Open the log file - temporary file unless the user decides to save it
      NLOG = 99
      CALL OPNLOG (NLOG)

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "algebra file_in file_out"')
        GOTO 130
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "algebra file_in file_out"')
        GOTO 130
      end if

C     Open the input database; Exit on error
      CALL get_argument(1,FILNAM, LFIL)
      ndbin = exopen(filnam(:lfil), EXREAD, cpuws, iows,
     &       vers, ierr)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//FILNAM(:LFIL)//'" does not exist.'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 130
      END IF

      call exinq(ndbin, EXDBMXUSNM, namlen, rdum, cdum, ierr)
      call exmxnm(ndbin, namlen, ierr)

      INOPEN = .TRUE.

C     Read the initial parameters from the database
C     ndbin  - file ID (input)
C     title  - title of input file dbtitl.blk
C     ndim   - number of coordinates per node dbnums.blk
C     numnp  - number of nodes dbnums.blk
C     numel  - number of elements dbnums.blk
C     nelblk - number of element blocks dbnums.blk
C     numnps - number of node sets dbnumg.blk
C     numess - number of side sets dbnumg.blk
C     ioerr   - error code
C     dbase.blk, dbtitl.blk, dbnums.blk, dbnumg.blk
      CALL EXGINI(NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &            NUMNPS, NUMESS, IERR)

C     if number of node sets > 0, gather node set information dbnumg.blk
      if (numnps .gt. 0) then
C        EXNSNL - request the length of the node set node list.
         CALL EXINQ (NDBIN, EXNSNL, LNPSNL, RDUM, CDUM, IERR)
C        EXNSDF - request the length of the node sets distribution factors list
         CALL EXINQ (NDBIN, EXNSDF, LNPSDF, RDUM, CDUM, IERR)
      else
         lnpsnl = 0
         lnpsdf = 0
      end if

C     if number of side sets > 0, gather side set information dbnumg.blk
      if (numess .gt. 0) then
C        EXSSNL request length of the side set node list
         CALL EXINQ (NDBIN, EXSSNL, LESSNL, RDUM, CDUM, IERR)
C        EXSSEL request the length of the side sets element list.
         CALL EXINQ (NDBIN, EXSSEL, LESSEL, RDUM, CDUM, IERR)
C        EXSSDF request the length of the side set distribution factors list
         CALL EXINQ (NDBIN, EXSSDF, LESSDF, RDUM, CDUM, IERR)
      else
         lessnl = 0
         lessel = 0
         lessdf = 0
      end if

C     Title of output db = title of input db dbtitl.blk
      TITLEO = TITLE

C     Print the database filename, title, and db initial variables
C     Cannot print global, element, or nodal info.  Data has not
C     been read from the input database yet.
      CALL DBPINI ('TIS', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &             NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL,
     &             LESSDF, IDUM, IDUM, IDUM, ' ')

cccc     DBLIST uses MDFIND to locate NUMELB, NUMLNK, NUMATR
cccc     Scan the connectivity to get the element block IDs
cccc     Read element block connectivity

C     call MxRSRV(array_name, ret_array_index, array_size)
C     IDELB: Element block ID's for each block
      CALL MDRSRV ('IDELB', KIDELB, NELBLK)
C     NUMELB: Number of elements in each block
      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
C     NUMLNK: Number of nodes per element in each block
      CALL MDRSRV ('NUMLNK', KNLNK, NELBLK)
C     NUMATR: Number of attributes in each block
      CALL MDRSRV ('NUMATR', KNATR, NELBLK)
C     VISELB(i) = True iff element block i is to be written
      CALL MDRSRV ('VISELB', KVISEB, NELBLK)
C     BLKTYP: Element type
      CALL MCRSRV ('BLKTYP', KNMLB, NELBLK * MXSTLN)
C     Reserving space for the element block connectivity arrays
      CALL MDRSRV ('LINK', KLINK, 0)
C     Reserving space for the element block attributes
      CALL MDRSRV ('ATRIB', KATRIB, 0)
C     Check for dynamic memory errors
      CALL MDSTAT (NERR, MEM)
      CALL MCSTAT (CERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0)) THEN
         CALL MEMERR
         GOTO 130
      END IF

C     DBIELB reads and returns the following:
C     1. Element block ID's
C     2. Element type in each element block
C     3. Number of elements in each element block
C     4. Number of nodes per element in each element block
C     5. Number of attributes per element in each element block
C     6. If Option = A or * -> element block attributes
C     7. If Option = C or * -> element block connectivity
C     8. If Option = H or * -> read all header information
C     9. If Option = I or * -> read the element block ID's
C     Note: Element block ID's must be read for all options
C     DBIELB reads the element block information from the database.
C     An error message is displayed if the end of file is read.
      CALL DBIELB (NDBIN, '*', 1, NELBLK, A(KIDELB), A(KNELB),
     &     A(KNLNK), A(KNATR), C(KNMLB), A, IELNK, IEATR, MERR)
C     Exit program on memory error
      IF (MERR .EQ. 1) GO TO 130

C     QA and Information record number stored in dbnumq.blk
C     Request the number of QA records.  Return the value
C     as an integer in nqarec
      CALL EXINQ (NDBIN, EXQA, NQAREC, RDUM, CDUM, IERR)
C     Request the number of information records.  Return the
C     value as an integer in ninfo
      CALL EXINQ (NDBIN, EXINFO, NINFO, RDUM, CDUM, IERR)

C     Reserve contiguous block space for arrays QA and information records
C     Reserve nqarec for input file + 1 for the current run
      CALL MCGET(((nqarec+1)*4*MXSTLN) + (ninfo*MXLNLN))
C     Reserve space to read the QA and information records
      call MCRSRV('QAREC', kqarec, (nqarec+1) * 4 * MXSTLN)
      call MCRSRV('INFREC', kinfo, ninfo * MXLNLN)
      CALL MCSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         GOTO 130
      END IF

C     Read the QA and Information Records
      CALL DBIQA(NDBIN, '*', NQAREC, C(KQAREC), NINFO, C(KINFO))

C     Read the number of global, node, and element variables
C     Read number of global variables = nvargl
      call exgvp(ndbin, 'G', nvargl, ierr)

C     Read number of nodal variables = nvarnp
      call exgvp(ndbin, 'N', nvarnp, ierr)

C     Read number of element variables = nvarel
      call exgvp(ndbin, 'E', nvarel, ierr)

C     Reserve memory for coordinate array names
      CALL MCRSRV ('NAMECO', KNACOR, namlen*NDIM)
C     Reserve memory for global, node, and element variable
      CALL MCRSRV ('NAMES', KNAMES, namlen*(NVARGL+NVARNP+NVAREL))
C     Reserve memory for element variable truth table
      CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
C     Temporary storage
      CALL MDRSRV ('ITMP', KITMP,  NELBLK * NVAREL)
      CALL MCSTAT (CERR, MEM)
      CALL MDSTAT (NERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0))  THEN
         CALL MEMERR
         MERR = 1
         GOTO 130
      END IF

C     Read the coordinate names
      CALL DBICON (NDBIN, NDIM, C(KNACOR))
      IF (IOERR .EQ. 1) GO TO 130

C     Read the names of the global, node, and element variables
      CALL DBINAM (NDBIN, C, C(KNAMES), NVARGL, NVARNP, NVAREL,
     &             KNAMGV, KNAMNV, KNAMEV, MAXVAR)

C     Print database variable names
C     call dbpnam (option, num_glb_var, num_nod_var, num_elem_var,
C                  glb_var_names, nod_var_names, elem_var_names)
      CALL DBPNAM ('*', NVARGL, NVARNP, NVAREL,
     &             C(KNAMES+NAMLEN*(KNAMGV-1)),
     &             C(KNAMES+NAMLEN*(KNAMNV-1)),
     &             C(KNAMES+NAMLEN*(KNAMEV-1)))

C     Read the element variable truth table
      CALL DBIVTT (NDBIN, A(KIEVOK), A(KITMP), NELBLK, NVAREL)

C     Delete temporary dynamic memory
      CALL MDDEL('ITMP')
      CALL MDSTAT (NERR, MEM)
      if (NERR .GT. 0) then
         CALL MEMERR
         merr = 1
         GO TO 130
      end if

C     Request the number of time steps from the database
      CALL EXINQ (NDBIN, EXTIMS, NSTEPS, RDUM, CDUM, IERR)
C     Reserve memory to hold the time step values
C     Reserve memory for step times
      CALL MDRSRV ('TIMES', KTIMES, max(1,NSTEPS))
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        GO TO 130
      END IF
      if (nsteps .gt. 0) then
        CALL EXGATM(NDBIN, A(KTIMES), IERR)
      else
C        Add a dummy step at time 0 in case user is trying to add a new step
        A(KTIMES) = 0.0
      END IF

C     Displays the number of time steps and the minimum and
C     maximum time on the database
C     'NM'   Print the number of time steps and the minimum and
C            maximum time step times
C     NSTEPS The number of time steps
C     A(KTIMES) The database time steps
      CALL DBPTIM ('NM', NSTEPS, A(KTIMES))

C     Selected time steps
      CALL MDRSRV ('IPTIMS', KPTIMS, MAX(1,NSTEPS))
C     selected element blocks
      CALL MDRSRV ('SELELB', KSELEB, NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         MERR = 1
         GOTO 130
      END IF

C ... See if the input database has any nodeset or sideset variables.
C     If there are any, then output a warning message telling the user
C     that these are not supported in algebra and will be stripped from
C     output file.
      call exgvp(ndbin, 'M', nvarns, ierr)
      call exgvp(ndbin, 'S', nvarss, ierr)
      if (nvarns .gt. 0) then
         CALL PRTERR('WARNING',
     $        'The Input Database contains NODESET VARIABLES')
         CALL PRTERR('CMDSPEC',
     $        '      NODESET VARIABLES are not supported by algebra')
         CALL PRTERR('CMDSPEC',
     $        '      These will NOT be saved to the output database')
      end if
      if (nvarss .gt. 0) then
         CALL PRTERR('WARNING',
     $        'The Input Database contains SIDESET VARIABLES')
         CALL PRTERR('CMDSPEC',
     $        '      SIDESET VARIABLES are not supported by algebra')
         CALL PRTERR('CMDSPEC',
     $        '      These will NOT be saved to the output database')
      end if

C     Read in the equations
C     call rdeqns(DMarray, DCarray, coord_names, elem_names, db_var_names,
C     qa_rec, info_rec, db_time_step, select_time_step,
C     elem_blk_ids, TRUE_IFF_elem(i)_written, DMindex_ISEVOK,
C     max_stack_size, return_if_quit)

      CALL RDEQNS (A, C, C(KNACOR), C(KNMLB), C(KNAMES), C(KQAREC),
     &             C(KINFO),  A(KTIMES), A(KPTIMS), A(KSELEB),
     &             A(KIDELB), A(KVISEB), KIEVOK, MAXSTK, NLOG, MERR)
      IF (MERR .EQ. 1) GO TO 130

C     RDEQNS processes the input equations as follows:
C        o Reads the line and parses it into fields.
C        o Checks the equation for syntax.
C        o Stores the variable names from the equations.
C        o Puts the equation in postfix form.
C        o Adds the equation variables to a global list.
C     If any errors are found, that equation is ignored.
C     After all the equations have been read in, all the variables are
C     gathered into the /VAR../ arrays.  The assigned variables are
C     after the expression variables.

C     Sort the input variables (by ID)
C     SORTID gathers the variables of the specified type into
C     one area of an array and orders these variable
      ITIME = 1
      CALL SORTID ('T', NUMINP, ITIME, I)
      ICOBEG = ITIME + 1
      CALL SORTID ('C', NUMINP, ICOBEG, ICOEND)
      IGVBEG = ICOEND + 1
      CALL SORTID ('G', NUMINP, IGVBEG, IGVEND)
      INVBEG = IGVEND + 1
      CALL SORTID ('N', NUMINP, INVBEG, INVEND)
      IEVBEG = INVEND + 1
      CALL SORTID ('E', NUMINP, IEVBEG, IEVEND)

C     Set the ID to the order in which the LHS variables were defined
      NLHS = 1
      DO 100 I = MAXVAR, IXLHS, -1
         IDVAR(I) = NLHS
         NLHS = NLHS + 1
  100 CONTINUE

C     Sort deleted variables to end of entries
C     SORDEL sorts the deleted variables so that they appear at the start
C     of all non-deleted entries.
      CALL SORDEL (MAXVAR, IXLHS, I)
      JTIME = I + 1

C     Sort the LHS variables in the order defined
      CALL SORTID ('T', MAXVAR, JTIME, I)
      JGVBEG = I + 1
      CALL SORTID ('G', MAXVAR, JGVBEG, JGVEND)
      JNVBEG = JGVEND + 1
      CALL SORTID ('N', MAXVAR, JNVBEG, JNVEND)
      JEVBEG = JNVEND + 1
      CALL SORTID ('E', MAXVAR, JEVBEG, JEVEND)
      NVARGO = JGVEND - JGVBEG + 1
      NVARNO = JNVEND - JNVBEG + 1
      NVAREO = JEVEND - JEVBEG + 1

C     Link the variables with storage
C     Assign storage for variables
C     LNKSTO sets up the storage locations for all the variables.
C     Input and output variables of the same name share storage (unless
C     one is a history/global and the other is not).  Time and history/global
C     variables are all in the first storage location: time is in slot 1,
C     followed by the input history variables (if any), the input global
C     variables (if any), then the output only history/global variables.
      CALL LNKSTO (C(KNAMES+NAMLEN*(KNAMGV-1)), NUMSTO, LTMENT, IOERR)
      IF (IOERR .EQ. 1) GO TO 130

C     LNKFNC sets up the storage locations for the time functions that
C     need storage for results that must be saved over time steps.
      CALL LNKFNC (NUMSTO, *130)

C     Relink the equation variables with the sorted name array
      CALL LNKVAR (*130)

C     Check displacement variables
      IF ((NVARNP .GT. 0) .OR. (NVARNO .GT. 0)) THEN
C     CHKDIS finds the displacement variables.  The first two/three nodal
C     variables are displacement variables if and only if they begin with
C     'D' and end with the last character of the corresponding coordinate
C     name.
         CALL CHKDIS (NDIM, C(KNACOR), NVARNO, NAMVAR(JNVBEG),
     *    namlen, maxnam)
      END IF

C *************************************************************
C                   Open the Output Database
C *************************************************************
      CALL get_argument(2,FILNAM, LFIL)
      ndbout = excre(filnam(:lfil), EXCLOB, CPUWS, IOWS, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Problems creating database "'//FILNAM(:LFIL)//'".'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        goto 130
      END IF
      OTOPEN = .TRUE.
      call exmxnm(ndbout, namlen, ierr)

C     RWEVAL reads the unread information from the input database,
C     processes the data for zoom mesh, and write the output database.

      CALL RWEVAL (NDBIN, NDBOUT, A, A, C, NPTIMS,
     &             NUMSTO, LTMENT, MAXSTK, NWRIT, IOERR, MERR)
      IF ((IOERR .EQ. 1) .OR. (MERR .EQ. 1)) GO TO 130

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
      END IF

C     Converts a set of integer numbers into a consistent
C     set of strings of the same length (right justified)
      CALL INTSTR (1, 0, NWRIT, STR8, LSTR)
      WRITE (*, 10000) STR8(:LSTR)
10000 FORMAT (/, 4X, A,
     &   ' time steps have been written to the database')

  130 CONTINUE

      call mddel ('IDELB')
      CALL MDDEL ('NUMELB')
      CALL MDDEL ('VISELB')
      CALL MCDEL ('BLKTYP')
      CALL MDDEL ('SELELB')
      CALL MDDEL ('TIMES')
      CALL MDDEL ('IPTIMS')
      CALL MCDEL ('NAMECO')
      CALL MCDEL ('NAMES')
      CALL MDDEL ('ISEVOK')

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
      END IF

C     Exiting programs - close input/output files if they are open
C     or even if they opened and may have an error
      IF (INOPEN) CALL exclos(NDBIN, IERR)

      IF (OTOPEN) CALL exclos(NDBOUT, IERR)

C     Performs finishing details common to all programs.
C     Specifically, it gets and displays the CPU time used.

      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

      END
