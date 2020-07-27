C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      PROGRAM EXPLORE
C=======================================================================

C                         *** EXPLORE 2.00 ***

C   --*** EXPLORE *** (EXPLORE) GENESIS/EXODUS database examination program
C   --
C   --EXPLORE is a post-processing program to examine the output of a
C   --finite element analysis, which is in the GENESIS or EXODUS database
C   --format.  EXPLORE allows the user to examine any values in the database.
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

C   --Developed at Sandia National Laboratories.
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C   --
C   --Documentation:
C   --   "User's Manual for EXPLORE"

      include 'exodusII.inc'
      include 'exp_progqa.blk'
      include 'exp_outfil.blk'
      include 'exp_dbase.blk'
      include 'exp_dbtitl.blk'
      include 'exp_dbnums.blk'
      include 'argparse.inc'

C      --A - the dynamic numeric memory base array
      DIMENSION A(1)
      INTEGER IA(1)
      EQUIVALENCE (A(1), IA(1))
      CHARACTER*1 C(1)

      CHARACTER*2048 DBNAME
      CHARACTER*2048 SCRATCH
C      --DBNAME - the database name, needed because the database may be closed
      character*256  option, value

      LOGICAL EXODUS
C      --EXODUS - true iff EXODUS file versus GENESIS file
      LOGICAL MAPND, MAPEL

      LOGICAL ISEOF
      LOGICAL CHECK
      CHARACTER*1 cdum

      include 'exp_qainfo.blk'
      CALL STRTUP (QAINFO)

C   --Set up the print file

      NCRT = -1
      NOUT = NCRT
      NPRT = 20
      ANYPRT = .FALSE.

C   --Print banner to CRT and print file

      CALL BANNER (0, QAINFO,
     &   'A GENESIS/EXODUS DATABASE EXPLORATION PROGRAM',
     &   ' ', ' ')

      call cpyrgt (0, "2008")

C   --Open the database

      NDB = 11

      CMPSIZ = 0
      IOWS   = 0
      DBNAME  = ' '

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .eq. 0) then
        CALL PRTERR ('FATAL', 'Filename not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "explore [-[no]map node|element|all] filename"')
        GOTO 120
      end if

      CALL get_argument(narg,DBNAME, LNAM)
      NDB = exopen(dbname(:lnam), EXREAD, CMPSIZ, IOWS, vers, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//DBNAME(:LNAM)//'" does not exist.'
         CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
         GOTO 120
      END IF

      CALL EXOPTS(EXVRBS,IERR)

C ... By default, ultimately map both nodes and elements
C     HOWEVER, in the transition time do not map either unless requested...
      mapel = .false.
      mapnd = .false.
      check = .false.

      if (narg .gt. 1) then
        do i=1, narg-1, 2
          CALL get_argument(i+0,option, lo)
          CALL get_argument(i+1,value,  lv)
          if (option(:lo) .eq. '-nomap' .or.
     *      option(:lo) .eq. '--nomap') then
            if (value(1:1) .eq. 'n' .or. value(1:1) .eq. 'N')
     *        mapnd = .false.
            if (value(1:1) .eq. 'e' .or. value(1:1) .eq. 'E')
     *        mapel = .false.
            if (value(1:1) .eq. 'a' .or. value(1:1) .eq. 'A') then
              mapnd = .false.
              mapel = .false.
            end if
          else if (option(:lo) .eq. '-map' .or.
     *      option(:lo) .eq. '--map') then
            if (value(1:1) .eq. 'n' .or. value(1:1) .eq. 'N')
     *        mapnd = .true.
            if (value(1:1) .eq. 'e' .or. value(1:1) .eq. 'E')
     *        mapel = .true.
            if (value(1:1) .eq. 'a' .or. value(1:1) .eq. 'A') then
              mapnd = .true.
              mapel = .true.
            end if
          else if (option(:lo) .eq. '-check' .or.
     *      option(:lo) .eq. '--check') then
            check = .TRUE.
          end if
        end do
      end if

        write (*,9999)
 9999   FORMAT(/,
     *    1x,'NOTE: This version has the option to use global',
     *    ' ids for both node and element ids.',/,
     *    1x,'      To see the mapping from local to global, use',
     *    ' the commands:',/,
     *    1x,'          "LIST MAP" (element map), or ',
     *    '"LIST NODEMAP" (node map)',/,
     *    1x,'      To disable the maps and use local ids, restart',
     *    ' explore with "-nomap node|element|all"',//,
     *    1x,'      To enable the maps and use global ids, restart',
     *    ' explore with "-map node|element|all"',//,
     *    1x,'      Notify gdsjaar@sandia.gov if bugs found')

        if (mapel .and. mapnd) then
          WRITE (*, 10010) 'Nodes and Elements using Global Ids'
        else if (mapel) then
          WRITE (*, 10010) 'Elements use Global Ids, Node Ids are Local'
        else if (mapnd) then
          WRITE (*, 10010) 'Element use Local Ids, Node Ids are Global'
        else
          WRITE (*, 10010) 'Nodes and Elements using Local Ids'
        end if
10010   FORMAT (/, 1X, 5A)

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
         call exinq(ndb, EXSSNL, lessnl, rdum, cdum, ierr)
         IF (IERR .NE. 0) GOTO 110
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
      if (mapel) then
        kdbmapel = kmapel
      else
        call mdrsrv ('DBMAPEL', kdbmapel, numel)
      endif
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      if (mapel) then
        CALL RDMAP (NDB, NUMEL, IA(KMAPEL), ISEOF)
      else
        CALL RDMAP (NDB, NUMEL, IA(KdbMAPEL), ISEOF)
        call iniseq(numel, ia(kmapel))
      end if

C ... Read node map
      CALL MDRSRV ('MAPNO', KMAPNO, NUMNP)
      if (mapnd) then
        kdbmapno = kmapno
      else
        call mdrsrv ('DBMAPNO', kdbmapno, numnp)
      endif
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
      if (mapnd) then
        CALL RDNMAP (NDB, NUMNP, IA(KMAPNO), ISEOF)
      else
        CALL RDNMAP (NDB, NUMNP, IA(KDBMAPNO), ISEOF)
        call iniseq(numnp, ia(kmapno))
      end if

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

C ... If NELBLK .gt. 0 and NUMEL .eq. 0, then
C     possible reading an experimental exodus file
C     containing only results data...
C     Calculate NUMEL based on NUMELB entries.
      if (NELBLK .gt. 0 .and. NUMEL .eq. 0) then
        do i=1, NELBLK
          numel = numel + IA(KNELB+I-1)
        end do
      end if

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
10000       FORMAT (1X, 'Number of time steps on the database =', I12)
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
     $     C(KVNAMO+NAMLEN*(IXSS-1)), A(KCORD),
     *     IA(KMAPEL), IA(KDBMAPEL), IA(KMAPNO), IA(KDBMAPNO),
     *     mapnd, mapel, check,
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

      CALL MDDEL ('CORD')
      CALL MCDEL ('NAMECO')
      CALL MDDEL ('MAPEL')
      CALL MDDEL ('MAPNO')
      CALL MDDEL ('IDELB')
      CALL MDDEL ('NUMELB')
      CALL MDDEL ('NUMLNK')
      CALL MDDEL ('NUMATR')
      CALL MDDEL ('LENE')
      CALL MCDEL ('EBTYPE')
      CALL MCDEL ('EBNAME')
      CALL MDDEL ('IDNPS')
      CALL MDDEL ('NNNPS')
      CALL MDDEL ('NDNPS')
      CALL MDDEL ('IXNNPS')
      CALL MDDEL ('IXDNPS')
      CALL MDDEL ('LTNNPS')
      CALL MDDEL ('FACNPS')
      CALL MCDEL ('NSNAME')
      CALL MDDEL ('LISNP')
      CALL MDDEL ('NLISEL')
      CALL MDDEL ('LISEL')
      CALL MDDEL ('LISBEL')
      CALL MDDEL ('LISNPS')
      CALL MDDEL ('LISESS')
      CALL MDDEL ('LINK')
      CALL MDDEL ('ATRIB')
      call mcdel ('ATRNM')
      if (exodus) then
        call mddel('LISGV')
        call mddel('LISNV')
        call mddel('LISEV')
        call mddel('LISMV')
        call mddel('LISSV')
        call mddel('VARGL')
        call mddel('VARNP')
        call mddel('VAREL')
        call mddel('VARNS')
        call mddel('VARSS')
        call mddel('TIMES')
        call mcdel('VNAMEI')
        call mcdel('VNAMEO')
        CALL MDDEL('ISEVOK')
        CALL MDDEL('ISNSVOK')
        CALL MDDEL('ISSSVOK')
      endif
      if (.not. mapnd) then
        call MDDEL('DBMAPNO')
      end if
      if (.not. mapel) then
        call MDDEL('DBMAPEL')
      end if
      call mcdel('QAREC')
      call mcdel('INFREC')

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

      subroutine iniseq(icnt, map)
      integer map(*)
      do i=1, icnt
        map(i) = i
      end do
      return
      end
