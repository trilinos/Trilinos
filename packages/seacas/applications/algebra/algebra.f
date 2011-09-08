C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
C=======================================================================
      PROGRAM ALGEBRA2
C=======================================================================
C $Id: algebra2.f,v 1.16 2009/03/25 14:42:08 gdsjaar Exp $
C
C   --This version of ALGEBRA will read and write EXODUSIIV2 database
C   --format files.  Many changes have occurred since the first version
C   --of ALGEBRA.  The original database files, genesis and exodusI
C   --were sequential access files.  EXODUSIIV2 uses a random access
C   --file format. Previous verions of ALGEBRA would have to read the
C   --input database more than once in order to get the file pointer
C   --to the desired data.  With random access files we are able to
C   --select what we want to read or write at anytime.
C
C                         *** ALGEBRA 2.02 ***
C   --*** ALGEBRA *** (ALGEBRA) Algebraic Database Manipulation Program
C   --   Written by Amy Gilkey - revised 05/18/88
C   --   Modified by Greg Sjaardema - 4/17/91
C   --     o Initial unix version
C   --     o Removed need for slatec library - wrote princ3.f
C   --   Modified by Frank Mello - 4/01/92
C   --     o Added _ as valid character in variable names
C   --   Modified by Greg Sjaardema - 4/01/92
C   --     o Hopefully fixed problem with deleting element blocks
C   --   Modified by Christi Forsythe - 8/26/95
C   --     o Read and write ExodusIIV2 database format
C   --
C   --The ALGEBRA program allows the user to process data from a finite
C   --element analysis before it is plotted. The finite element output
C   --data is in the form of variable values (stress, strain, and
C   --velocity components, etc.) in an EXODUS database.  The ALGEBRA program
C   --evaluates user-supplied functions of the data and writes the results
C   --to an output EXODUS database which can be read by plot programs.
C   --
C   --Expected input:
C   --   o The equations and commands on the standard input device.
C   --   o The input EXODUSIIV2 database on unit 11.
C   --
C   --Output:
C   --   o A listing of the input database information and any errors
C   --     found on the standard output device.
C   --   o The output EXODUSIIV2 database on unit 12.

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                    ISSUED BY SANDIA LABORATORIES,                   *
*                      A PRIME CONTRACTOR TO THE                      *
*                  UNITED STATES DEPARTMENT OF ENERGY                 *
* * * * * * * * * * * * * *   N O T I C E   * * * * * * * * * * * * * *
* This program was prepared as an account of work sponsored by the    *
* United States Government.  Neither the United States nor the United *
* States Department of Energy nor any of their employees, nor any of  *
* their contractors, subcontractors, or their employees, makes any    *
* warranty, express or implied, or assumes any legal liability or     *
* responsibility for the accuracy, completeness or usefulness of any  *
* information, apparatus, product or process disclosed, or represents *
* that its use would not infringe privately owned rights.             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

C   --Developed at Sandia National Laboratories.
C   --
C   --Current author and code sponsor: Amy Gilkey
C   --
C   --Revision History:
C   --   8/95   Converted from ExodusIIv1 to ExodusIIv2 database format
C   --   11/87  Converted from SEACO to EXODUS database (Amy Gilkey)
C   --   11/84  New sponsor (Amy Gilkey)
C   --   05/82  Changed to CRAY (M.A. Richgels)
C   --   11/81  Altered to FORTRAN 77 (M.A. Richgels)
C   --   11/81  Changed to VAX interactive (M.A. Richgels)
C   --   06/80  Created (M.A. Richgels)
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   Some of the SUPES subroutines have been rewritten for this code
C   --   SUPES package (dynamic memory, FORTRAN extensions)
C   --   SLATEC mathematics library
C   --
C   --Runs on VAX VMS !#VAX
C#CTSSC   --Runs on CRAY CTSS

C   --Documentation: SAND86-0881, printed May 1986
C   --   "ALGEBRA - A Program That Algebraically Manipulates the Output
C   --   of a Finite Element Analysis"

      include 'exodusII.inc'
      include 'namlen.blk'
      
C     max_num_equation, max_parsed_entries
      include 'numeqn.blk'
C     QAINFO array with program information
      include 'progqa.blk'
C     input equation info
      include 'ent.blk'
C     input variable info
      include 'var.blk'
C     aliases
      include 'alias.blk'
C     I/O file names
      include 'dbase.blk'
C     I/O database titles
      include 'dbtitl.blk'
C     num_of: nodes,elements,coordinated/node,element_blks
C     num_of: history,global,nodal,element variables, num_of_time_steps
      include 'dbnums.blk'
C     node set/side sets number,length,dist factor
      include 'dbnumg.blk'
C     database type, num_of qa and info records
      include 'dbnumq.blk'
C     num_of: nodes,elements,element blks,node_sets,side_sets in zoom mesh
C     num_of: output history,global,nodal,element variables
      include 'dbout.blk'
C     time_index, I_index/O_index for history,global,nodal and element vars
      include 'dbxvar.blk'
C     time variables
      include 'times.blk'
C     zoom info
      include 'zoom.blk'
C     function variables
      include 'fnctbc.blk'
C     equation line error messages
      include 'eqnlns.blk'
C     floating point byte size
      include 'dbws.blk'
      include 'f2kcli.inc'

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
      CHARACTER*256 FILNAM, SCRATCH

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
      include 'qainfo.blk'

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
      NARG = COMMAND_ARGUMENT_COUNT()
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
      CALL GET_COMMAND_ARGUMENT(1,FILNAM, LFIL, ISTATUS)
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
     &             NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &             LESSDF, IDUM, IDUM, IDUM)

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
      IF (nsteps .gt. 0) THEN
C        Reserve memory for step times
         CALL MDRSRV ('TIMES', KTIMES, NSTEPS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            GO TO 130
         END IF
         CALL EXGATM(NDBIN, A(KTIMES), IERR)
       ELSE
         KTIMES=1
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
      CALL GET_COMMAND_ARGUMENT(2,FILNAM, LFIL, ISTATUS)
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

C     Exiting programs - close input/output files if they are open
C     or even if they opened and may have an error
      IF (INOPEN) CALL exclos(NDBIN, IERR)

      IF (OTOPEN) CALL exclos(NDBOUT, IERR)

C     Performs finishing details common to all programs.
C     Specifically, it gets and displays the CPU time used.

      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

      END


