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
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++ Copyright 1988, Sandia Corporation. The United States Government
C+++ retains a limited license in this software as prescribed in AL 88-1
C+++ and AL 91-7. Export of this program may require a license from
C+++ the United States Government.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C=======================================================================
      PROGRAM GROPE
C=======================================================================

C                         *** GROPE 2.00 ***

C   --*** GROPE *** (GROPE) GENESIS/EXODUS database examination program
C   --
C   --GROPE is a post-processing program to examine the output of a
C   --finite element analysis, which is in the GENESIS or EXODUS database
C   --format.  GROPE allows the user to examine any values in the database.
C   --The display can be directed to the CRT or to a print file.
C   --
C   --Expected input:
C   --   o The commands on the standard input device.
C   --   o The GENESIS/EXODUS database on unit 11.
C   --
C   --Output:
C   --   o A listing of the input database information and any errors
C   --     found on the standard output device.
C   --   o A print file of requested information on unit 20.

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
C   --   01/88  Changed to SELECT/LIST/PRINT structure
C   --   10/87  Added EXODUS database (Amy Gilkey)
C   --   04/86  Created (Amy Gilkey)
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C   --
C   --Runs on VAX VMS !#VAX
C#CTSSC   --Runs on CRAY CTSS

C   --Documentation:
C   --   "User's Manual for GROPE"

      include 'exodusII.inc'
      INCLUDE 'progqa.blk'
      INCLUDE 'outfil.blk'
      INCLUDE 'dbase.blk'
      INCLUDE 'dbtitl.blk'
      INCLUDE 'dbnums.blk'
      INCLUDE 'f2kcli.inc'

C      --A - the dynamic numeric memory base array
      DIMENSION A(1)
      INTEGER IA(1)
      EQUIVALENCE (A(1), IA(1))
      CHARACTER*1 C(1)

      CHARACTER*132 DBNAME
      CHARACTER*132 SCRATCH
C      --DBNAME - the database name, needed because the database may be closed

      LOGICAL EXODUS
C      --EXODUS - true iff EXODUS file versus GENESIS file

      LOGICAL ISEOF
      CHARACTER*1 cdum

      INCLUDE 'qainfo.blk'
      CALL STRTUP (QAINFO)

C   --Set up the print file

      NCRT = -1
      NOUT = NCRT
      NPRT = 20
      ANYPRT = .FALSE.

C   --Print banner to CRT and print file

      CALL BANNER (0, QAINFO,
     &   'A GENESIS/EXODUS DATABASE EXAMINATION PROGRAM',
     &   ' ', ' ')

      call cpyrgt (0, "2008")

C   --Open the database

      NDB = 11

      CMPSIZ = 0
      IOWS   = 0
      DBNAME  = ' '

C .. Get filename from command line.  If not specified, emit error message
      NARG = COMMAND_ARGUMENT_COUNT()
      if (narg .eq. 0) then
        CALL PRTERR ('FATAL', 'Filename not specified.')
        CALL PRTERR ('FATAL', 'Syntax is: "grope filename"')
        GOTO 120
      else if (narg .gt. 1) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL', 'Syntax is: "grope filename"')
        GOTO 120
      end if

      CALL GET_COMMAND_ARGUMENT(1,DBNAME, LNAM, ISTATUS)
      NDB = exopen(dbname(:lnam), EXREAD, CMPSIZ, IOWS, vers, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//DBNAME(:LNAM)//'" does not exist.'
         CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
         GOTO 120
      END IF

C   --Initialize dynamic memory
      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Read the initial variables
      call exgini(ndb, title, ndim, numnp, numel, nelblk,
     *     numnps, numess, ierr)
      IF (IERR .NE. 0) GOTO 110
      if (numnps .gt. 0) then
         call exinq(ndb, EXNSNL, lnpsnl, rdum, cdum, ierr)
         IF (IERR .NE. 0) GOTO 110
         call exinq(ndb, EXNSDF, lnpsdf, rdum, cdum, ierr)
         IF (IERR .NE. 0) GOTO 110
      else
         lnpsnl = 0
         lnpsdf = 0
      end if
      if (numess .gt. 0) then
         lessnl = 0
         call exinq(ndb, EXSSEL, lessel, rdum, cdum, ierr)
         IF (IERR .NE. 0) GOTO 110
         call exinq(ndb, EXSSDF, lessdf, rdum, cdum, ierr)
         IF (IERR .NE. 0) GOTO 110
      else
         lessnl = 0
         lessel = 0
         lessdf = 0
      end if
      
      call exinq(ndb, EXDBMXUSNM, namlen, rdum, cdum, ierr)
      IF (IERR .NE. 0) GOTO 110
      call exmxnm(ndb, namlen, ierr)
      IF (IERR .NE. 0) GOTO 110
      
      CALL PRINIT ('NTISC', NOUT, DBNAME, TITLE,
     &     NDIM, NUMNP, NUMEL, NELBLK,
     &     NUMNPS, LNPSNL, lnpsdf, NUMESS, LESSEL, LESSNL, LESSDF,
     &     NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)

C ... See if there are any timesteps on the database (is EXODUS)
      call exinq (ndb, EXTIMS, NSTEPS, rdum, cdum, ierr)
      IF (IERR .NE. 0) GOTO 110
      EXODUS = (NSTEPS .gt. 0)

      if (EXODUS) THEN
        CALL EXGVP (NDB,"G",NVARGL,IERR)
        IF (IERR .NE. 0) GOTO 110
        CALL EXGVP (NDB,"E",NVAREL,IERR)
        IF (IERR .NE. 0) GOTO 110
        CALL EXGVP (NDB,"N",NVARNP,IERR)
        IF (IERR .NE. 0) GOTO 110
        CALL EXGVP (NDB,"M",NVARNS,IERR)
        IF (IERR .NE. 0) GOTO 110
        CALL EXGVP (NDB,"S",NVARSS,IERR)
        IF (IERR .NE. 0) GOTO 110

        CALL PRINIT ('V', NOUT, DBNAME, TITLE,
     &       NDIM, NUMNP, NUMEL, NELBLK,
     &       NUMNPS, LNPSNL, lnpsdf, NUMESS, LESSEL, LESSNL, LESSDF,
     &       NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)
      END IF
      CALL SETPRC(4,0)

C ... Read coordinate data
      CALL MDRSRV ('CORD', KCORD, NUMNP * NDIM)
      CALL MCRSRV ('NAMECO', KNMCO, NAMLEN*NDIM)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      CALL RDCORD (NDB, NDIM, NUMNP, A(KCORD), C(KNMCO), ISEOF, NAMLEN)

C ... Read element map
      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      CALL RDMAP (NDB, NUMEL, A(KMAPEL), ISEOF)

C ... Read node map
      CALL MDRSRV ('MAPNO', KMAPNO, NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      CALL RDNMAP (NDB, NUMNP, A(KMAPNO), ISEOF)

C ... Read element blocks
      CALL MDRSRV ('IDELB', KIDELB, NELBLK)
      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
      CALL MDRSRV ('NUMLNK', KNLNK, NELBLK)
      CALL MDRSRV ('NUMATR', KNATR, NELBLK)
      CALL MDRSRV ('LENE', KLENE, 1+NELBLK)
      CALL MCRSRV ('EBTYPE', KNMLB, MXSTLN*NELBLK)
      CALL MCRSRV ('EBNAME', KNMEB, NAMLEN*NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      CALL RDELB (NDB, NELBLK,
     &   A(KIDELB), A(KNELB), A(KNLNK), A(KNATR), 
     &   A, C, KLINK, KATRIB, KATRNM, ISEOF, C(KNMLB), C(KNMEB),
     &   NAMLEN)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C ... Read nodesets
      CALL MDRSRV ('IDNPS',  KIDNS,  NUMNPS)
      CALL MDRSRV ('NNNPS',  KNNNS,  NUMNPS)
      CALL MDRSRV ('NDNPS',  KNDNPS, NUMNPS)
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
      CALL MDRSRV ('IXDNPS', KIXDNS, NUMNPS)
      CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
      CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)
      CALL MCRSRV ('NSNAME', KNMNS,  NAMLEN*NUMNPS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      CALL RDNPS (NDB, NUMNPS, LNPSNL,
     &     A(KIDNS), A(KNNNS), A(KNDNPS), A(KIXNNS), A(KIXDNS),
     &     A(KLTNNS), A(KFACNS), C(KNMNS), ISEOF, NAMLEN)

C ... Read sidesets
      CALL MDRSRV ('IDESS',  KIDSS,  NUMESS)
      CALL MDRSRV ('NEESS',  KNESS,  NUMESS)
      CALL MDRSRV ('NDESS',  KNDSS,  NUMESS)
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
      CALL MDRSRV ('IXDESS', KIXNSS, NUMESS)
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
      CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
      CALL MDRSRV ('FACESS', KFACSS, LESSDF)
      CALL MCRSRV ('SSNAME', KNMSS,  NAMLEN*NUMESS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      CALL RDESS (NDB, NUMESS, LESSEL, LESSDF,
     &   A(KIDSS), A(KNESS), A(KNDSS), A(KIXESS), A(KIXNSS),
     &   A(KLTESS), A(KLTSSS), A(KFACSS), C(KNMSS), ISEOF,
     &   NAMLEN)

C ... Read QA Information
      CALL RDQA (NDB, NQAREC, NINFO, KQAREC, KINFO, C)

      if (exodus) then
C ...    Read variable names and truth table
         CALL RDNAME (A, C, NDB, KVNAMI, KVNAMO,
     &        IXGV, IXNV, IXEV, IXNS, IXSS, 
     &        KIEVOK, KNSVOK, KSSVOK)
         
C ... Read in the times for all the time steps from the database
         call mdrsrv('TIMES', KTIMES, NSTEPS)
         CALL RDTIMS (NDB, A(KTIMES))

         WRITE (*, *)
         IF (NSTEPS .GT. 0) THEN
            WRITE (*, 10000, IOSTAT=IDUM) NSTEPS
10000       FORMAT (1X, 'Number of time steps on the database =', I10)
         END IF
         
C ... Get the memory for the variables
         CALL MDRSRV ('VARGL', KVARGL, NVARGL)
         CALL MDRSRV ('VARNP', KVARNP, NVARNP * NUMNP)
         CALL MDRSRV ('VAREL', KVAREL, NVAREL * NUMEL)
         CALL MDRSRV ('VARNS', KVARNS, NVARNS * LNPSNL)
         CALL MDRSRV ('VARSS', KVARSS, NVARSS * LESSEL)
         
      ELSE
         NSTEPS = 0
         KVARGL = 1
         KVARNP = 1
         KVAREL = 1
         KVARNS = 1
         KVARSS = 1
      END IF

C   --Get the memory for the logical arrays

      CALL MDRSRV ('LISNP', KLISNP, 1+NUMNP)
      CALL MDRSRV ('NLISEL', KNLISE, 1+NELBLK)
      CALL MDRSRV ('LISEL', KLISEL, 1+NUMEL)
      CALL MDRSRV ('LISBEL', KLISBE, 1+NUMEL)
      CALL MDRSRV ('LISNPS', KLISNS, 1+NUMNPS)
      CALL MDRSRV ('LISESS', KLISSS, 1+NUMESS)
      IF (EXODUS) THEN
         CALL MDRSRV ('LISGV', KLISGV, 1+NVARGL)
         CALL MDRSRV ('LISNV', KLISNV, 1+NVARNP)
         CALL MDRSRV ('LISEV', KLISEV, 1+NVAREL)
         CALL MDRSRV ('LISMV', KLISMV, 1+NVARNS)
         CALL MDRSRV ('LISSV', KLISSV, 1+NVARSS)
      ELSE
         KLISGV = 1
         KLISNV = 1
         KLISEV = 1
         KLISMV = 1
         KLISSV = 1
      END IF

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

C   --Process commands

      CALL COMAND (A, IA, EXODUS, DBNAME, C(KQAREC), C(KINFO),
     &     C(KNMCO), C(KNMLB), C(KNMEB),C(KATRNM),
     &     C(KVNAMI+NAMLEN*(IXGV-1)), C(KVNAMI+NAMLEN*(IXNV-1)),
     &     C(KVNAMI+NAMLEN*(IXEV-1)), C(KVNAMI+NAMLEN*(IXNS-1)),
     $     C(KVNAMI+NAMLEN*(IXSS-1)),
     &     C(KVNAMO+NAMLEN*(IXGV-1)), C(KVNAMO+NAMLEN*(IXNV-1)),
     &     C(KVNAMO+NAMLEN*(IXEV-1)), C(KVNAMO+NAMLEN*(IXNS-1)),
     $     C(KVNAMO+NAMLEN*(IXSS-1)),
     &     A(KCORD), A(KMAPEL), A(KMAPNO),
     &     A(KIDELB), A(KNELB), A(KLENE), A(KNLNK), A(KNATR),
     &     A(KLINK), A(KATRIB),
     &     A(KIDNS), A(KNNNS), A(KNDNPS), A(KIXNNS), A(KIXDNS),
     $     A(KLTNNS), A(KFACNS), C(KNMNS), 
     &     A(KIDSS), A(KNESS), A(KNDSS), A(KIXESS), A(KIXNSS),
     &     A(KLTESS), A(KLTSSS), A(KFACSS), C(KNMSS), 
     &     A(KIEVOK), A(KNSVOK), A(KSSVOK), A(KTIMES), 
     &     A(KVARGL), A(KVARNP), A(KVAREL), A(KVARNS), A(KVARSS), 
     &     A(KLISNP), A(KNLISE), A(KLISEL), A(KLISNS), A(KLISSS),
     &     A(KLISGV), A(KLISNV), A(KLISEV), A(KLISMV), A(KLISSV))
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      
      GOTO 110
      
 100  CONTINUE
      CALL MEMERR
      GOTO 110
      
 110  CONTINUE
      
      call exclos(ndb, ierr)

 120  CONTINUE
      IF (ANYPRT) CLOSE (NPRT, IOSTAT=IDUM)
      
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))
      
      END
