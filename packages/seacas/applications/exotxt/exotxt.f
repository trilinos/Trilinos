C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      PROGRAM EXOTXT
C=======================================================================

C   --*** EXO2TXT *** Renamed for ExodusIIv2 database
C   --*** EXOTXT *** (EXOTXT) EXODUS to TEXT translator
C   --   Written by Amy Gilkey - revised 03/02/88
C   --   Modified for ExodusIIv2 database format
C   --
C   --EXOTXT reads either from an EXODUS database or from the user
C   --and writes a text file of the database contents.
C   --
C   --Expects the input database on unit 9, the output text file on unit 20.

      include 'exodusII.inc'
      INCLUDE 'argparse.inc'


C     Input/Output File Arguments
C     CPUWS - The word size in bytes of the floating point variables
C             used in the application program
C     IOWS  - The word size in bytes of the floating point data as they
C             are stored in the EXODUS II file
C     IERR  - error code
C     NERR  - error flag for dynamic memory errors - numeric
C     CERR  - error flag for dynamic memory errors - character
C     IOERR - Input/output error flag
      INTEGER CPUWS, IOWS, IERR, NERR, CERR, IOERR
C     VERS  - version number of database
      REAL VERS
C     QA program information
      CHARACTER*(MXSTLN) QAINFO(6)
C     Title of database/experiment
      CHARACTER*(MXLNLN) TITLE
C     Output string for number of time steps processed
      CHARACTER*8 STR8
C     Input filename
      CHARACTER*2048 FILNAM, SCRATCH
C     LDAT - length of date information
C     LREV - length of revision information to actual revision number
      INTEGER LDAT, LREV

C     A(1) dynamic memory base array for numeric data
C     C(1) dynamic memory base array for character data
      DIMENSION A(1)
      CHARACTER*1 C(1)

C     Set ExodusII error reporting level
C     Set level to EXABRT for production use
      call exopts(EXABRT, ierr)

C     Program Information
C.
      QAINFO(1) = 'exotxt                          '
      QAINFO(2) = '2012/11/07                      '
      QAINFO(3) = ' 1.22                           '
      QAINFO(4) = '                                '
      QAINFO(5) = '                                '
      QAINFO(6) = '                                '

      CALL STRTUP (QAINFO)

      CALL BANNER (0, QAINFO,
     &   'EXODUSII DATABASE TO TEXT FILE TRANSLATOR',
     &   ' ', ' ')

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      CALL MCSTAT (CERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0)) THEN
         CALL MEMERR
         GOTO 140
      END IF

C   --Open the input and output files

      NDB = 11
      NTXT = 20

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filename not specified.')
        CALL PRTERR ('FATAL', 'Syntax is: "exotxt db_file text_file"')
        GOTO 140
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL', 'Syntax is: "exotxt db_file text_file"')
        GOTO 140
      end if

C     Open the input database; Exit on error
      CALL get_argument(1,FILNAM, LFIL)
      CPUWS = 0
      IOWS  = 0
      NDB = EXOPEN(FILNAM(:LFIL), EXREAD, CPUWS, IOWS,
     &       VERS, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//FILNAM(:LFIL)//'" does not exist.'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 140
      END IF

      CALL get_argument(2,FILNAM, LFIL)
      open(unit=ntxt, file=filnam(:lfil), iostat=ierr)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Could not create "'//FILNAM(:LFIL)//'"'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 140
      END IF

C   --Read the initial variables

      CALL EXGINI (NDB, TITLE, NNDIM, NUMNP, NUMEL, NELBLK,
     &             NUMNPS, NUMESS, IERR)

      call exinq(ndb, EXDBMXUSNM, namlen, rdum, cdum, ierr)
      call exmxnm(ndb, namlen, ierr)

C     Request length of the concatenated node set node list
      CALL EXINQ (NDB, EXNSNL, LNPSNL, RNUM, CNUM, IERR)
C     Request length of the node set distribution factors list
      CALL EXINQ (NDB, EXNSDF, LNPSDF, RDUM, CDUM, IERR)
C     Request length of the concatenated side sets element list
      CALL EXINQ (NDB, EXSSEL, LESSEL, RNUM, CNUM, IERR)
C     Request length of the concatenated side sets node list
      CALL EXINQ (NDB, EXSSNL, LESSNL, RNUM, CNUM, IERR)
C     Request length of the side set distribution factors list
      CALL EXINQ (NDB, EXSSDF, LESSDF, RNUM, CNUM, IERR)

      CALL WRINIT (NTXT, VERS, TITLE, NNDIM, NUMNP, NUMEL, NELBLK,
     &             NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &             LESSDF,QAINFO, NAMLEN)
      CALL DBPINI ('TIS', NDB, TITLE, NNDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSDF,
     &   IDUM, IDUM, IDUM, ' ')

C   --Read the coordinates

      CALL MDRSRV ('XN', KXN, NUMNP)
      IF (NNDIM .GE. 2) CALL MDRSRV ('YN', KYN, NUMNP)
      IF (NNDIM .GE. 3) CALL MDRSRV ('ZN', KZN, NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         GOTO 140
      END IF

C     Reserve memory for coordinate array names
      CALL MCRSRV ('NAMECO', KNACOR, NAMLEN*NNDIM)

      call getxyz(ndb, c(knacor), a(kXN), A(KYN), A(KZN), IERR, namlen)
      CALL WRXYZ (NTXT, NNDIM, NUMNP, A(KXN), A(KYN), A(KZN),
     *  C(knacor), namlen)

      CALL MDDEL ('XN')
      IF (NNDIM .GE. 2) CALL MDDEL ('YN')
      IF (NNDIM .GE. 3) CALL MDDEL ('ZN')
      CALL MCDEL ('NAMECO')

C   --ExodusIIv2 Read node number map, element number map,
C   --and element order map
      CALL MDRSRV ('NPMAP', KNPMAP, NUMNP)
      CALL MDRSRV ('ELMAP', KELMAP, NUMEL)
      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)   
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         GOTO 140
      END IF

C     Read node number map
      CALL EXGNNM (NDB, A(KNPMAP), IERR)
C     Read element number map
      CALL EXGENM (NDB, A(KELMAP), IERR)
C     Read element order map
      CALL EXGMAP (NDB, A(KMAPEL), IERR)

      CALL WRMAP (NTXT, '*', NUMNP, NUMEL,
     &             A(KNPMAP), A(KELMAP), A(KMAPEL))

      CALL MDDEL ('NPMAP')
      CALL MDDEL ('ELMAP')
      CALL MDDEL ('MAPEL')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
         CALL MEMERR
         GOTO 140
      END IF

C   --Read/write the element blocks
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
C     NAMELB: Element type
      CALL MCRSRV ('NAMELB', KNMLB, NELBLK * MXSTLN)
C     Reserving space for the element block connectivity arrays
      CALL MDRSRV ('LINK', KLINK, 0)
C     Reserving space for the element block attributes
      CALL MDRSRV ('ATRIB', KATRIB, 0)
C     Check for dynamic memory errors
      CALL MDSTAT (NERR, MEM)
      CALL MCSTAT (CERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0)) THEN
         CALL MEMERR
         GOTO 140
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
      CALL DBIELB (NDB, '*', 1, NELBLK, A(KIDELB), A(KNELB),
     &     A(KNLNK), A(KNATR), C(KNMLB), A, IELNK, IEATR, IOERR)
C     Exit program on error reading element block information
      IF (IOERR .EQ. 1) GO TO 140
      CALL MDFIND ('LINK', KLINK, IELNK)
      CALL MDFIND ('ATRIB', KATRIB, IEATR)
      CALL WRELB (NTXT, NELBLK, A(KIDELB), A(KNELB), A(KNLNK),
     &     A(KNATR), C(KNMLB), A(KLINK), A(KATRIB))

      CALL MDDEL ('NUMLNK')
      CALL MDDEL ('NUMATR')
      CALL MDDEL ('LINK')
      CALL MDDEL ('VISELB')
      CALL MDDEL ('ATRIB')
      CALL MCDEL ('NAMELB')
      CALL MDSTAT (NERR,MEM)
      CALL MCSTAT (CERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0)) THEN
         CALL MEMERR
         GO TO 140
      END IF

C   --Read/write the node sets

C     if number of node set > 0
      IF (NUMNPS .GT. 0) THEN
C     IDNPS  - array containing the node set ID's for each node set
C     NNNPS  - array containing the number of nodes for each node set
C     NDNPS  - array containing number of dist. fact for each node set
C     IXNNPS - array containing indices into the LTNNPS array which
C              are the location of the 1st nodes for each set
C     IXDNPS - array containing indices into the FACNPS array which
C              are the location of the 1st dist factor for each set
C     LTNNPS - Returned array containing the nodes for all node sets
C              Internal node IDs
C     FACNPS - Returned array containing the distribtion factors
C              for all node sets
         CALL MDRSRV ('IDNPS',   KIDNS, NUMNPS)
         CALL MDRSRV ('NNNPS',   KNNNS, NUMNPS)
         CALL MDRSRV ('NDNPS',  KNDNPS, NUMNPS)
         CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
         CALL MDRSRV ('IXDNPS', KDISNS, NUMNPS)
         CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
         CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)
         CALL MDSTAT (NERR,MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            GO TO 140
         END IF
C        Read the concatenated node sets
         CALL EXGCNS (ndb, a(kidns), a(knnns), a(kndnps),
     &                a(kixnns), a(kdisns), a(kltnns),
     &                a(kfacns), IERR)
         CALL WRNPS (NTXT, NUMNPS, LNPSNL, LNPSDF, A(KIDNS),
     &               A(KNNNS), A(KNDNPS), A(KIXNNS), A(KDISNS),
     &               A(KLTNNS), A(KFACNS))

         CALL MDDEL ('IDNPS')
         CALL MDDEL ('NNNPS')
         CALL MDDEL ('NDNPS')
         CALL MDDEL ('IXNNPS')
         CALL MDDEL ('IXDNPS')
         CALL MDDEL ('LTNNPS')
         CALL MDDEL ('FACNPS')
         CALL MDSTAT (NERR,MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            GO TO 140
         END IF
      END IF


C   --Read/write the element side sets

      IF (NUMESS .GT. 0) THEN
C     IDESS  - array containing side set IDS
C     NEESS  - array containing the number of sides for each sets
C     KNDSS  - Returned array containing the number of dist 
C              factors for each set
C     IXEESS - returned array containing the indices into the 
C              LTEESS array which are the locations of the 1st
C              element of each set
C     IXNESS - Returned array containing the indices into the
C              FACESS array which are the locations of the 1st
C              distribution factor for each set.
C     LTEESS - Returned array containing the elements for all side
C              sets. Internal element IDS are used in this list
C     LTSESS - Returned array containing the sides for all side sets
C     FACESS - Returned aray containing dist factors for all side sets
C     NNESS  - the number of nodes for each side set
C     IXNESS - index into LTNESS - the 1st node for each side set
C     LTNESS - array of nodes for all side sets
C     LTNNN  - array of number of nodes for each side in a side set
         CALL MDRSRV ('IDESS' , KIDSS , NUMESS)
         CALL MDRSRV ('NEESS' , KNESS , NUMESS)
         CALL MDRSRV ('NDESS' , KNDSS , NUMESS)
         CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
         CALL MDRSRV ('IXDESS', KIDESS, NUMESS)
         CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
         CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
         CALL MDRSRV ('FACESS', KFACSS, LESSNL)
         CALL MDRSRV ('NNESS' , KNNSS,  NUMESS)
         CALL MDRSRV ('IXNESS', KIXNSS, NUMESS)
         CALL MDRSRV ('LTNESS', KLTNSS, LESSNL)
         CALL MDRSRV ('LTNNN' , KLTNNN, LESSEL)
         CALL MDSTAT(NERR, MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            GO TO 140
         END IF

C        Read concatenated side sets
         CALL EXGCSS (NDB, A(KIDSS), A(KNESS), A(KNDSS),
     &                A(KIXESS), A(KIDESS), A(KLTESS),
     &                A(KLTSSS), A(KFACSS), IOERR)

C        Convert sides to nodes
         CALL DBIGN (NDB, NUMESS, A(KIDSS), A(KNNSS), 
     &               A(KIXNSS), A(KLTNSS), A(KLTNNN), IOERR)
         IF (IOERR .EQ. 1) GO TO 140


         CALL WRESS (NTXT, NUMESS, LESSEL, LESSNL, LESSDF,
     &        A(KIDSS), A(KNESS), A(KNDSS), A(KIXESS), A(KIDESS),
     &        A(KLTESS), A(KLTSSS), A(KFACSS), A(KNNSS), A(KIXNSS),
     &        A(KLTNSS), A(KLTNNN))

         CALL MDDEL ('IDESS')
         CALL MDDEL ('NEESS')
         CALL MDDEL ('NDESS')
         CALL MDDEL ('IXEESS')
         CALL MDDEL ('IXDESS')
         CALL MDDEL ('LTEESS')
         CALL MDDEL ('LTSESS')
         CALL MDDEL ('FACESS')
         CALL MDDEL ('NNESS')
         CALL MDDEL ('IXNESS')
         CALL MDDEL ('LTNESS')
         CALL MDDEL ('LTNNN')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            GOTO 140
         END IF

      END IF

      write (ntxt, '(A)') '! Properties'
C     Element Block, Node Set, and Side Set Properties
C     number of element block properties
      call exinq(ndb, EXNEBP, numebp, rdum, cdum, ierr)
C     number of node set properties
      call exinq(ndb, EXNNSP, numnsp, rdum, cdum, ierr)
C     number of side set properties
      call exinq(ndb, EXNSSP, numssp, rdum, cdum, ierr)
C************************************************************************
C     Read element block properties
C************************************************************************
C     Reserve memory for element block properties names and values
      call mcrsrv ('EBPNAM', iebpn, numebp * namlen)
      call mdrsrv ('EBPVAL', iebpv, numebp * nelblk)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         ioerr = 1
      end if
      write (NTXT, '(I10,12x,A)') numebp, 
     &      '! Number of ELEMENT BLOCK Properties'
      call wrprop (ndb, NTXT, EXEBLK, numebp, nelblk,
     &             c(iebpn), a(iebpv), namlen)

      call mcdel ('EBPNAM')
      call mddel ('EBPVAL')

C************************************************************************
C     read/write node set properties
C************************************************************************
      call mcrsrv ('NSPNAM', inspn, numnsp * namlen)
      call mdrsrv ('NSPVAL', inspv, numnsp * numnps)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         ioerr = 1
      end if

      write (NTXT, '(I10,12x,A)') numnsp, 
     &      '! Number of NODE SET Properties'
      call wrprop (ndb, NTXT, EXNSET, numnsp, numnps,
     &             c(inspn), a(inspv), namlen)

      call mcdel ('NSPNAM')
      call mddel ('NSPVAL')

C************************************************************************
C     read/write side set properties
C************************************************************************

      call mcrsrv ('SSPNAM', isspn, numssp * namlen)
      call mdrsrv ('SSPVAL', isspv, numssp * numess)
      call mcstat (cerr, mem)
      call mdstat (nerr, mem)
      if ((nerr .ne. 0) .or. (cerr .ne. 0)) then
         call memerr()
         ioerr = 1
      end if

      write (NTXT, '(I10,12x,A)') numssp, 
     &      '! Number of SIDE SET Properties'
      call wrprop (ndb, NTXT, EXSSET, numssp, numess,
     &             c(isspn), a(isspv), namlen)

      call mcdel ('SSPNAM')
      call mddel ('SSPVAL')
      call mdstat (nerr, mem)
      if (nerr .ne. 0) then
         call memerr()
         ioerr = 1
      end if


C   --Read the QA records
C     QA and Information record number stored in dbnumq.blk
C     Request the number of QA records.  Return the value
C     as an integer in nqarec
      call exinq(ndb, EXQA, nqarec, rdum, cdum, ierr)
C     Request the number of information records.  Return the
C     value as an integer in ninfo
      call exinq(ndb, EXINFO, ninfo, rdum, cdum, ierr)

C     Reserve space to read the QA and information records
      CALL MCRSRV ('QAREC', KQAREC, (NQAREC+1)*4*MXSTLN)
      CALL MCRSRV ('INFREC', KINFO, NINFO*MXLNLN)
      CALL MCSTAT (CERR, MEM)
      IF (CERR .GT. 0) THEN
         CALL MEMERR
         GOTO 140
      END IF

      CALL DBIQA (NDB, '*', NQAREC, C(KQAREC), NINFO, C(KINFO))
      CALL WRQA (NTXT, NQAREC, C(KQAREC), NINFO, C(KINFO), QAINFO)

      CALL MCDEL ('QAREC')
      CALL MCDEL ('INFREC')
      CALL MCSTAT (CERR, MEM)
      IF (CERR .GT. 0) THEN
         CALL MEMERR
         GOTO 140
      END IF


C   --Read the database names

C     Read the number of global, node, and element variables
C     EXGVP reads the number of global, nodal, or element
C     variables stored in the database
C     Read number of global variables = nvargl
      call exgvp(ndb, 'G', nvargl, ierr)
C     Read number of nodal variables = nvarnp
      call exgvp(ndb, 'N', nvarnp, ierr)
C     Read number of element variables = nvarel
      call exgvp(ndb, 'E', nvarel, ierr)

C     Reserve memory for global, node, and element variable
      CALL MCRSRV ('NAMES', KNAMES, NAMLEN*(NVARGL+NVARNP+NVAREL))
C     Reserve memory for element variable truth table
      CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
      CALL MDRSRV ('ITMP', KITMP,  NELBLK * NVAREL)
      CALL MCSTAT (CERR, MEM)
      CALL MDSTAT (NERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0))  THEN
         CALL MEMERR
         GOTO 140
      END IF

C     Read the names of the global, node, and element variables
      CALL DBINAM (NDB, C, C(KNAMES), NVARGL, NVARNP, NVAREL,
     &             KNAMGV, KNAMNV, KNAMEV, IOERR, NAMLEN)
      IF (IOERR .EQ. 1) GO TO 140

      CALL DBPINI ('V', NDB, TITLE, NNDIM, NUMNP, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSDF,
     &      NVARGL, NVARNP, NVAREL, ' ')

C     Read the element variable truth table
      CALL DBIVTT (NDB, A(KIEVOK), A(KITMP), NELBLK, NVAREL)

      CALL WRNAME (NTXT, NNDIM, NELBLK, NVARGL, NVARNP, NVAREL,
     &  C(KNAMES+NAMLEN*(KNAMGV-1)), C(KNAMES+NAMLEN*(KNAMNV-1)),
     &  C(KNAMES+NAMLEN*(KNAMEV-1)), A(KIEVOK), NAMLEN)

C     Delete temporary dynamic memory
      CALL MCDEL('NAMES')
      CALL MDDEL('ITMP')
      CALL MCSTAT (CERR, MEM)
      CALL MDSTAT (NERR, MEM)
      if ((NERR .GT. 0) .OR. (CERR .GT. 0)) then
         CALL MEMERR
         GO TO 140
      end if

C   --Read the number of database time steps
      call exinq(ndb, EXTIMS, nsteps, rdum, cdum, ierr)

      CALL MDRSRV ('VARGL', KVARGL, NVARGL)
      CALL MDRSRV ('VARNP', KVARNP, NVARNP * NUMNP)
      CALL MDRSRV ('VAREL', KVAREL, NVAREL * NUMEL)
      CALL MDSTAT (NERR, MEM)
      if (NERR .GT. 0) then
         CALL MEMERR
         GO TO 140
      end if

C ... Zero out memory to account for variables not read from DB 
C     due to truth table...
      CALL INIREA(NVAREL*NUMEL, 0.0, A(KVAREL))
      CALL INIREA(NVARNP*NUMNP, 0.0, A(KVARNP))

      WRITE (*, *)
      WRITE (*, *)

      WRITE (*, 10020) NSTEPS
10020 FORMAT (' ', I8, ' time steps on the input database') 
      DO 110 ISTEP = 1, NSTEPS
         CALL DBISTE (NDB, '*', ISTEP, NELBLK, TIME, 
     &                NVARGL, NVARNP, NVAREL, NUMNP, 
     &                A(KIDELB), A(KNELB), A(KIEVOK),
     &                A(KVARGL), A(KVARNP), A(KVAREL), IOERR)
         IF (IOERR .EQ. 1) GO TO 140

         CALL WRSTEP (NTXT, ISTEP, NELBLK, TIME,
     &                NVARGL, NVARNP, NVAREL, NUMNP,
     &                A(KIDELB), A(KNELB), A(KIEVOK),
     &                A(KVARGL), A(KVARNP), A(KVAREL),
     &                MAX(1, NVARGL), MAX(1, NVARNP),
     &                MAX(1, NVAREL), IOERR)
         IF (IOERR .EQ. 1) GO TO 140

         WRITE (*, 10000) ISTEP
10000     FORMAT (I8, ' time steps processed') 
  110 CONTINUE

      WRITE (STR8, '(I8)', IOSTAT=K) NSTEPS
      CALL SQZSTR (STR8, LSTR)
      WRITE (*, 10010) STR8(:LSTR)
10010  FORMAT (/, 4X, A,
     &   ' time steps have been written to the text file')

C     Delete dynamic memory
      CALL MDDEL ('IDELB')
      CALL MDDEL ('NUMELB')
      CALL MDDEL ('ISEVOK')
      CALL MDDEL ('VARGL')
      CALL MDDEL ('VARNP')
      CALL MDDEL ('VAREL')
      CALL MDSTAT (NERR, MEM)
      if (NERR .GT. 0) CALL MEMERR

  140 CONTINUE

      CLOSE (NTXT, IOSTAT=K)
      call exclos(ndb, ierr)
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

      END


C ... Written as wrapper to get string lengths correct on coordinate
C     name array which is dynamically allocated      
      SUBROUTINE GETXYZ(NDB, NAMECO, X, Y, Z, ierr, namlen)
      character*(NAMLEN) NAMECO(*)
      real x(*), y(*), z(*)
      CALL EXGCOR(NDB, X, Y, Z, ierr)
      call exgcon(ndb, nameco, ierr)
      return
      end
 
