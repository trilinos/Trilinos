C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c=======================================================================
      PROGRAM BLOTII
C=======================================================================

C                         *** BLOT ***

C   --*** BLOT *** (BLOT) Post-processing plot program
C   --
C   --BLOT is a graphics program for post-processing of finite element
C   --analyses output in the EXODUS II database format.  BLOT combines the
C   --plotting capabilities of DETOUR, TPLOT, and SPLOT.  It is command
C   --driven with free-format input.  BLOT can drive any graphics device
C   --supported by the Sandia Virtual Device Interface.
C   --
C   --DETOUR produces mesh plots with various representations of the
C   --analysis output variables. The major capabilities of DETOUR are
C   --deformed mesh plots, line contours, filled (painted) contours,
C   --vector plots of two variables (e.g., velocity vectors), and symbol
C   --plots of scalar variables (e.g., discrete cracks). DETOURs features
C   --include element selection by element block, element birth and death,
C   --multiple views for combining several displays on each plot, and
C   --symmetry mirroring.
C   --
C   --TPLOT generates time-versus-variable plots or variable-versus-variable
C   --plots where the variables are user-selected global, nodal or element
C   --EXODUS database variables.  TPLOTs features include multiple curve
C   --plots and neutral file output.
C   --
C   --SPLOT generates distance-versus-variable plots at selected time steps
C   --where the distance is the accumulated distance between pairs of nodes
C   --or element centers and the variable is a nodal or element EXODUS
C   --database variable.  SPLOTs features include multiple curve plots
C   --and neutral file output.
C   --
C   --Plot times from the database are processed in one of the following modes:
C   --   Uniform time interval - plot at every DELT time interval from
C   --      TMIN to TMAX.
C   --   All available times - plot at all available times from TMIN
C   --      to TMAX.
C   --   User selected times - plot at user selected times.
C   --The closest available time step is chosen.
C   --
C   --Expected input:
C   --   o The commands on the standard input device.
C   --   o The input EXODUS database on unit 11.
C   --
C   --Output:
C   --   o A listing of the input database information and any errors
C   --     found on the standard output device.
C   --   o The plots on the specified graphics device.
C   --   o A GRAFAID neutral file on unit 20.
C   --   o A EXPLORE listing file on unit 21.

C   --Developed at Sandia National Laboratories.
C   --
C   --Current author and code sponsor: John Glick
C   --
C   --Revision History:
C   --   05/88  Added PATHLINE (Amy Gilkey)
C   --   03/88  Added Master/Slave logic (Amy Gilkey)
C   --   10/87  Added EXPLORE (Amy Gilkey)
C   --   10/87  Converted from SEACO to EXODUS database (Amy Gilkey)
C   --   07/87  Combined DETOUR, TPLOT, and SPLOT (Amy Gilkey)
C   --DETOUR:
C   --   05/86  Started to add 3D (Amy Gilkey)
C   --   03/85  New sponsor (Amy Gilkey)
C   --   11/82  Created (Dennis Flanagan)
C   --TPLOT:
C   --   xx/xx  Neutral file code (Greg Sjaardema)
C   --   01/86  New sponsor (Amy Gilkey)
C   --   xx/xx  Modified for the VAX11/780 (Johnny Biffle)
C   --   08/80  Created (Zelma Beisinger)
C   --SPLOT:
C   --   xx/xx  Neutral file code (Greg Sjaardema)
C   --   01/86  New sponsor (Amy Gilkey)
C   --   xx/xx  Created (Mary Sagartz)
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   SVDI graphics package
C   --   PLT graphics package
C   --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C   --
C   --Runs on VAX VMS, sun, Unicos, Ultrix
C   --
C   --Documentation: SAND86-0914, printed July 1987                       f
C   --   "DETOUR - A Deformed Mesh / Contour Plot Program"
C   --Documentation: SAND86-0883, printed August 1986
C   --   "TPLOT - A Time History or X-Y Plot Program for the Output
C   --   of a Finite Element Analysis"
C   --Documentation: SAND86-0882, printed August 1986
C   --   "SPLOT - A Distances-versus-Variable Plot Program for the Output
C   --   of a Finite Element Analysis"

C   --NOTE: All parameters and common areas are defined in BLKDAT

      include 'exodusII.inc'

C      --These parameters define the indices of 2D and 3D limit arrays
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

C      --These parameters define the mesh display (see MSHLIN of /MSHOPT/)
      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)

      PARAMETER (MXQARC=MXSTLN)
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbname.blk'
      include 'dbtitl.blk'
      include 'dbnums.blk'
      include 'dbnumgq.blk'
      include 'd3nums.blk'
      include 'dbnams.blk'
      include 'times.blk'
      include 'layout.blk'
      include 'deform.blk'
      include 'mshlim.blk'

      include 'selne.blk'
      include 'neutr.blk'
      include 'csv.blk'
      include 'outfil.blk'
      include 'sizes.blk'

      include 'light.blk'
      include 'icrnbw.blk'
      include 'legopt.blk'

      include 'argparse.inc'

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

C ... Used for hardcopy device drivers...
      common /blotans/ BLTANS
      character*2 BLTANS

      DIMENSION A(1),IA(1)
      EQUIVALENCE (A(1),IA(1))
      CHARACTER*1 C(1)
C      --A - the dynamic memory base array

      character*2048 scratch, value
      character*256  option

      LOGICAL MESHOK, DTOK, LNOK, SPOK, TPOK
      LOGICAL MAPND, MAPEL
      CHARACTER*(MXSTLN) CURPRO
      INTEGER NEWELB
      CHARACTER NEWPRO

      INTEGER IDUM
      REAL RDUM
      CHARACTER*8 CDUM
C     external sample_handler
      EXTERNAL TPREAD
      EXTERNAL BLKDAT

C     --The compute word size and I/O word size
      INTEGER CMPSIZ, IOWS

C ... Initialize Rainbow and Light (ICRNBW and LIGHT common blocks)
C     Block Data doesn't seem to work reliably on all systems

C ... Lights are stored as x, y, z, brightness, Vector is normalized in
C     shade.f.  NLIT is the number of lights.
      RMULT = 1.0
      GMULT = 1.0
      BMULT = 1.0
      LITE(1,1) = -1.0
      LITE(2,1) =  1.0
      LITE(3,1) =  1.0
      LITE(4,1) =  1.0
      LITE(5,1) = -0.57735
      LITE(6,1) =  0.57735
      LITE(7,1) =  0.57735
      LITE(8,1) =  1.0
      NLIT = 1
      AMBIENT = 0.2

      cdebug = '        '
      idebug = 0
      CMPSIZ = 0
      IOWS   = 0
      NEUTRL = 0

      NDB = 11
      NEU = 20
      NEUOPN = .FALSE.
      NEUGRF = 22
      GRFOPN = .FALSE.
      NCSV = 23
      CSVOPN = .FALSE.
      call version(qainfo)

C     Initialize variables
      VERS = 0.0
      IDUM = 0

      CALL STRTUP (QAINFO)

      draw(1) = '                                '
      draw(2) = '                                '
      draw(3) = '                                '
      draw(4) = '                                '
      DRAW(1) = QAINFO(1)
      DRAW(2) = '  ' // QAINFO(3)(3:8)
      DRAW(3)(1:8) = QAINFO(5)(1:8)
      DRAW(4)(1:8) = QAINFO(6)(1:8)

      WRITE (*, 10020)
      CALL BANNER (0, QAINFO,
     &  'A DEFORMED MESH / CONTOUR PLOT PROGRAM',
     &  'WITH X-Y PLOTTING CAPABILITIES',
     &  'FOR POST-PROCESSING OF FINITE ELEMENT ANALYSES'
     &  )
      CALL CPYRGT(0, '2009')

C   --Initialize dynamic memory

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDFILL (-999)
      CALL MCFILL ('Z')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

C   --Open database file

      NARG = argument_count()
      if (narg .eq. 0) then
        CALL PRTERR ('FATAL', 'Filename not specified.')
        CALL PRTERR ('CMDSPEC',
     *    'Syntax is: "blot.dev [-basename basename] [-ps_option num]'//
     *    ' [-nomap node|element|all] filename"')
        GOTO 190
      end if

      CALL get_argument(narg,dbname, lfil)
      dbname(lfil+1:) = ' '
      ndb = exopen(dbname(:lfil), EXREAD, CMPSIZ, IOWS, vers, ierr)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//dbname(:lfil)//'" does not exist.'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        CALL PRTERR ('CMDSPEC',
     *    'Syntax is: "blot.dev [-basename basename] [-ps_option num]'//
     *    ' [-nomap node|element|all] [-show_filename] filename"')
        GOTO 170
      END IF
      EXODUS = .FALSE.

      CALL INISTR (3, ' ', CAPTN(1,1))
      CALL INISTR (3, ' ', CAPTN(1,2))
      call exinq(ndb, EXDBMXUSNM, namlen, rdum, cdum, ierr)
      call exmxnm(ndb, namlen, ierr)

C ... Get basename of the database file to use for csv, neu and other output files
C     Assume first that the basename is from the database name.  It may be replaced
C     later if the user added a -hardcopy or -basename argument.
      last = indexr(dbname(:lfil), '.')
      if (last .gt. 2) then
        basenam = dbname(:last-1)
      else
        basenam = dbname
      end if

C ... Now see if a command-line option overrides this.
C     Options should all be of the form "-option arg"
      bltans = '7'

C ... By default, map both nodes and elements
      mapel = .true.
      mapnd = .true.

      if (narg .gt. 1) then
        i = 1
        do
          CALL get_argument(i,option, lo)
          i = i + 1
          if (option(:lo) .eq. '--show_filename' .or.
     *      option(:lo) .eq. '-show_filename') then
            captn(3,1) = dbname(:lfil)
            captn(3,2) = dbname(:lfil)

          else if (option(:lo) .eq. '-hardcopy' .or.
     *      option(:lo) .eq. '--hardcopy' .or.
     *      option(:lo) .eq. '-basename' .or.
     *      option(:lo) .eq. '--basename') then
            CALL get_argument(i,value,  lv)
            i = i + 1
            basenam = value(:lv)

          else if (option(:lo) .eq. '-ps_option' .or.
     *        option(:lo) .eq. '--ps_option') then
            CALL get_argument(i,value,  lv)
            i = i + 1
            if (lv .le. 2) then
              bltans = value(:lv)
            end if

          else if (option(:lo) .eq. '-nomap' .or.
     *        option(:lo) .eq. '--nomap') then
            CALL get_argument(i,value,  lv)
            i = i + 1
            if (value(1:1) .eq. 'n' .or. value(1:1) .eq. 'N')
     *        mapnd = .false.
            if (value(1:1) .eq. 'e' .or. value(1:1) .eq. 'E')
     *        mapel = .false.
            if (value(1:1) .eq. 'a' .or. value(1:1) .eq. 'A') then
              mapnd = .false.
              mapel = .false.
            end if
          end if
          if (i .gt. narg) exit
        end do
      end if

C   --Set error reporting level
      CALL EXOPTS(EXABRT,IERR)

C   --Initialize graphics
      CALL GRINIT (DBORD0, CHLSIZ)

C   --Read and print database header

      call exgini(ndb, title,
     *  ndim, numnp, numel, nelblk,
     *  numnps, numess, ierr)
      if (numnps .gt. 0) then
        call exinq(ndb, EXNSNL, lnpsnl, rdum, cdum, ierr)
        call exinq(ndb, EXNSDF, lnpsdf, rdum, cdum, ierr)
      else
        lnpsnl = 0
        lnpsdf = 0
      end if
      if (numess .gt. 0) then
        call exinq(ndb, EXSSNL, lessnl, rdum, cdum, ierr)
        call exinq(ndb, EXSSEL, lessel, rdum, cdum, ierr)
        call exinq(ndb, EXSSDF, lessdf, rdum, cdum, ierr)
      else
        lessnl = 0
        lessel = 0
        lessdf = 0
      end if

      NUMNPF = NUMNP

      CALL PRINIT ('NTIS', -1, NDB, DBNAME, TITLE,
     &  NDIM, NUMNP, NUMEL, NELBLK,
     &  NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &  NVARGL, NVARNP, NVAREL, NVARNS, NVARSS)

C   --Check whether mesh may be plotted

      CALL MSCHK (.FALSE., MESHOK)
      IF (.NOT. MESHOK)
     &  CALL PRTERR ('WARNING', 'No mesh is defined')

      IS3DIM = (NDIM .GE. 3)

C   --Read coordinates (full set)

C   --SCALER uses MDFIND to find XN, YN, ZN
      CALL MDGET (NDIM * NUMNP)
      CALL MDRSRV ('XN', KXN, NUMNP)
      CALL MDRSRV ('YN', KYN, NUMNP)
      IF (IS3DIM) THEN
        CALL MDRSRV ('ZN', KZN, NUMNP)
      ELSE
        KZN = 1
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      call exgcor(ndb, a(kxn), a(kyn), a(kzn), ierr)

C     Element Block IDs
      CALL MDRSRV ('IDELB', KIDELB, NELBLK)
C     Element Block IDs of HEXSHELL element blocks
      CALL MDRSRV ('HEXID', KHEXID, NELBLK)
C     Element Block types
      CALL MCRSRV ('NAMELB', KNMLB, NELBLK * MXSTLN)

C     Element Block names
      CALL MCRSRV ('EBNAME', KNMEB, NELBLK * NAMLEN)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      CALL EXGEBI (NDB, IA(KIDELB), IERR)
      if (nelblk .gt. 0) then
        CALL MDLONG ('IDELB', KIDELB, NELBLK)
        CALL MCLONG ('NAMELB', KNMLB, NELBLK * MXSTLN)
        CALL MCLONG ('EBNAME', KNMEB, NELBLK * NAMLEN)
      end if
C     Number of elements in element block
      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
C     Number of nodes per elements in element block
      CALL MDRSRV ('NUMLNK', KNLNKE, NELBLK)
C     Number of attributes in element block
      CALL MDRSRV ('NUMATR', KNATR, NELBLK)
C     Connectivity array
      CALL MDRSRV ('LINK', KLINKE, 0)
C     Connectivity pointer array
      CALL MDRSRV ('LPTR', KLPTR, NELBLK)
C     Attribute array
      CALL MDRSRV ('ATRIB', KATRIB, 0)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

C     Initialize element block integer arrays
      CALL INIINT(NELBLK, 0, IA(KIDELB))
      CALL INIINT(NELBLK, 0, IA(KNELB))
      CALL INIINT(NELBLK, 0, IA(KNLNKE))
      CALL INIINT(NELBLK, 0, IA(KNATR))

      CALL INISTR (NELBLK, ' ', C(KNMLB))

      CALL EXGEBI (NDB, IA(KIDELB), IERR)
C   --Read element block connectivity
      CALL DBIELB (NDB, '*', 1, NELBLK, IA(KIDELB), IA(KNELB),
     &  IA(KNLNKE), IA(KNATR), A, IA, KLINKE, KATRIB, C(KNMLB),
     &  C(KNMEB), IA(KLPTR), NAMLEN, *170)

C     Count the number of element blocks that contain HEXSHELLs. Store
C     the element block ID of HEXSHELL element blocks.
      ISHEX = 0
      CALL RDTYPE (NELBLK, C(KNMLB), IA(KIDELB), IA(KNELB),
     $  ISHEX, IA(KHEXID), NSHL)

      if (ishex .gt. 0) then
        nebsiz = nelblk + ishex
        CALL MDLONG ('IDELB', KIDELB, NEBSIZ)
        CALL MCLONG ('NAMELB', KNMLB, NEBSIZ * MXSTLN)
        CALL MCLONG ('EBNAME', KNMEB, NEBSIZ * NAMLEN)
C     Number of elements in element block
        CALL MDLONG ('NUMELB', KNELB, NEBSIZ)
C     Number of nodes per elements in element block
        CALL MDLONG ('NUMLNK', KNLNKE, NEBSIZ)
C     Number of attributes in element block
        CALL MDLONG ('NUMATR', KNATR, NEBSIZ)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

C     Initialize element block integer arrays
        CALL INIINT(ISHEX, 0, IA(KIDELB+NELBLK))
        CALL INIINT(ISHEX, 0, IA(KNELB+NELBLK))
        CALL INIINT(ISHEX, 0, IA(KNLNKE+NELBLK))
        CALL INIINT(ISHEX, 0, IA(KNATR+NELBLK))

C     Check for HEXSHELL - split HEXSHELL element block into a 'HEX'
C     element block and a 'SHELL' element block.
        CALL MDRSRV ('IDSCR', KIDSCR, NELBLK)
        do 350 j = 1, nelblk
          IA(KIDSCR+j-1) = IA(KIDELB+j-1)
 350    continue

        CALL PROCHS(A, IA, NELBLK, IA(KIDELB), IA(KIDSCR), IA(KNELB),
     &    IA(KNLNKE), IA(KNATR), KLINKE, KATRIB, C(KNMLB),
     &    IA(KLPTR), ISHEX, IA(KHEXID), *170)

C     Reset number of element blocks
        NELBLK = NELBLK + ISHEX
C     Reset number of elements
        INEL = 0
        do 20 I = 1, nelblk
          INEL = INEL + IA(KNELB+I-1)
 20     CONTINUE
        NUMEL = INEL

      ENDIF

C   --Scan element number map (global id)
      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160
      if (mapel .and. numel .gt. 0) then
        call exgenm (ndb, ia(kmapel), ierr)
      else
        call iniseq(numel, ia(kmapel))
      end if

C   --Read node number map (global id)
      CALL MDRSRV ('MAPND', KMAPND, NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160
      if (mapnd) then
        call exgnnm (ndb, ia(kmapnd), ierr)
      else
        call iniseq(numnp, ia(kmapnd))
      end if

C   --Change number of elements per element block to block index

C   --SCALER and MSMEMY and MSGEOM use MDFIND to find LENE
      CALL MDRSRV ('LENE', KLENE, 1+NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      CALL NUM2IX (NELBLK, A(KNELB), A(KLENE))

C   --Special fix for 8-node 2D elements - order nodes consecutively
      IF (.NOT. IS3DIM) THEN
        CALL ORD8NP (NELBLK, A(KLENE), A(KNLNKE), A(KLINKE))
      END IF

C   --Make up the element to element block index

C   --MSGEOM and MSSTEP use MDFIND to find IE2ELB
      CALL MDRSRV ('IE2ELB', KE2ELB, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      CALL MAKE2B (NELBLK, A(KLENE), A(KE2ELB))

C   --Read node sets and side sets

      CALL MDRSRV ('IDNPS',  KIDNS, NUMNPS)
      CALL MDRSRV ('NNNPS', KNNNS, NUMNPS)
      CALL MDRSRV ('NDNPS', KNDNPS, NUMNPS)
C   --DBLIST and SSMAIN use MDFIND to find IXNNPS, LTNNPS, and FACNPS
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
C   --In exgcns call, must have space even if no dist factors defined
      CALL MDRSRV ('IXDNPS', KIXDNS, NUMNPS)
      CALL MCRSRV ('NSNAME', KNMNS,  NUMNPS * NAMLEN)
      CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
      CALL MDRSRV ('FACNPS', KFACNS, LNPSDF)
      CALL MDRSRV ('IDESS', KIDSS, NUMESS)
      CALL MDRSRV ('NEESS', KNESS, NUMESS)
      CALL MDRSRV ('NNESS', KNNSS, NUMESS)
      CALL MDRSRV ('NDESS', KNDSS, NUMESS)
C   --DBLIST and SSMAIN use MDFIND to find IXEESS, IXNESS, LTEESS, LTNESS,
C   --and FACESS
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
      CALL MDRSRV ('IXNESS', KIXNSS, NUMESS)
      CALL MDRSRV ('IXDESS', KIXDSS, NUMESS)
      CALL MDRSRV ('LTNNSS', KLTNNN, LESSEL)
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
      CALL MDRSRV ('LTNESS', KLTNSS, LESSNL)
      CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
      CALL EXINQ (NDB, EXSSDF, LESSDL, rdum, cdum, ierr)
      CALL MDRSRV ('FACESS', KFACSS, LESSDL)
      CALL MCRSRV ('SSNAME', KNMSS,  NUMESS * NAMLEN)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      if (numnps .gt. 0) then
        call exgcns(ndb, a(kidns), a(knnns), a(kndnps), a(kixnns),
     *    a(kixdns), a(kltnns), a(kfacns), ierr)
        call getnams(ndb, EXNSET, numnps, c(knmns), namlen)
      end if
      if (numess .gt. 0) then
        call exgcssf(ndb, a(kidss), a(kness), a(kndss), a(kixess),
     *    a(kixdss), A(KLTeSS), a(kltsss), ierr)
        call getnams(ndb, EXSSET, numess, c(knmss), namlen)
        if (ierr .ne. 0) go to 170
      end if

C   --Read QA records

      call exinq(ndb, EXQA,   nqarec, rdum, cdum, ierr)
      call exinq(ndb, EXINFO, ninfo,  rdum, cdum, ierr)
      call mcrsrv('QAREC', kqarec, nqarec * 4 * MXSTLN)
      call mcrsrv('INFREC', kinfo, ninfo * MXLNLN)
      CALL MCSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160
      exodus = .true.
      if (nqarec .gt. 0) then
C ... Wrapper to get strings the right length
        call exgqaw(ndb, c(kqarec), ierr)
      end if
      if (ninfo .gt. 0) then
C ... Wrapper to get info record the right length
        call exginw(ndb, c(kinfo), ierr)
      end if

      CALL INISTR (4, ' ', CREATE)
      CALL INISTR (4, ' ', MODIFY)
C ... NOTE: cpyst8 must be called since cpystr knows that c() is only
C           a single character array.
      IF (NQAREC .GT. 0) CALL CPYST8 (4, C(KQAREC), CREATE)
      IF (NQAREC .GT. 1) THEN
        CALL CPYST8 (4, C(KQAREC+4*MXSTLN*(nqarec-1)), MODIFY)
      END IF

C   --Initialize (for GENESIS) and read database names

      NAMECO(1) = 'X'
      NAMECO(2) = 'Y'
      IF (IS3DIM) NAMECO(3) = 'Z'

      call exinq(ndb, EXTIMS, NTIMST, RDUM, CDUM, IERR)
      exodus = (ntimst .gt. 0)

      IF (EXODUS) THEN
C      --MSGEOM uses MDFIND to find ISEVOK (reserved in DBINAM)
        CALL DBINAM (NDB, 'CBVT', NDIM, NELBLK, NNDIM, NNELB,
     &    NVARHI, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &    NAMECO, KNAMHV,
     &    KNAMGV, KNAMNV, KNAMEV, KNAMNS, KNAMSS,
     &    A, IA, KIEVOK, C, KNAMES,
     &    EXODUS, IA(KIDELB), ISHEX, KHEXID, NAMLEN, *110)

        GOTO 120

C      --Handle error reading variable names
 110    CONTINUE
        NVARHI = 0
        NVARGL = 0
        NVARNP = 0
        NVAREL = 0
        NVARNS = 0
        NVARSS = 0
        EXODUS = .FALSE.
 120    CONTINUE

      ELSE
        NNDIM = -999
        NNELB = -999
      END IF

      NVARHI = MAX (0, NVARHI)
      NVARGL = MAX (0, NVARGL)
      NVARNP = MAX (0, NVARNP)
      NVAREL = MAX (0, NVAREL)
      NVARNS = MAX (0, NVARNS)
      NVARSS = MAX (0, NVARSS)

      IF (EXODUS) THEN
        CALL PRNAME (-1, namlen,
     *    NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &    C(KNAMES+NAMLEN*(KNAMGV-1)),
     &    C(KNAMES+NAMLEN*(KNAMNV-1)), C(KNAMES+NAMLEN*(KNAMEV-1)),
     &    C(KNAMES+NAMLEN*(KNAMNS-1)), C(KNAMES+NAMLEN*(KNAMSS-1)))
      END IF

C   --Get database time (allocates dynamic memory)
C   --SCALER uses MDFIND to find TIMES and WHOTIM

      if (EXODUS) then
        CALL DBITIM (NDB, '*', EXODUS,
     &    NVARNP, NELBLK, NVAREL, A(KIEVOK),
     &    NSTEPS, NSTEPW, A, KTIMES, KWHOLE, *170)
        CALL PRTIMS ('NM', -1, .TRUE., .TRUE.,
     &    NSTEPS, A(KTIMES), A(KWHOLE))
      else
C ... The 'TIMES' and 'WHOTIM' arrays are accessed even if there are
C     no timesteps on the model (See plcomd.f).
        call mdrsrv('TIMES', KTIMES, 1)
        call mdrsrv('WHOTIM', KWHOLE, 1)
      end if

      IF (.NOT. EXODUS) THEN
        WRITE (*, 10010) 'Database is in the GENESIS format'
      END IF

C   --Calculate the element centers

C      --SCALER uses MDFIND to find XE, YE, ZE
      CALL MDGET (NDIM * NUMEL)
      CALL MDRSRV ('XE', KXE, NUMEL)
      CALL MDRSRV ('YE', KYE, NUMEL)
      IF (IS3DIM) THEN
        CALL MDRSRV ('ZE', KZE, NUMEL)
      ELSE
        KZE = 1
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      CALL ELECOR (NDIM, NELBLK, A(KLENE), A(KNLNKE), A(KLINKE),
     &  A(KXN), A(KYN), A(KZN), A(KXE), A(KYE), A(KZE))

C   --Break 3D elements into faces, and sort faces by element block

      IF (MESHOK) THEN
        IF (IS3DIM) THEN
          WRITE (*, 10010)
     &      'Please wait while the geometry is processed'
        END IF

        IF (.NOT. IS3DIM) THEN
          CALL MDRSRV ('LENF', KLENF, 1 + NELBLK + 2)
          CALL INIINT (1+NELBLK+2, 0, IA(KLENF))
        ELSE
          CALL MDRSRV ('LENF', KLENF, 1 + NELBLK + 4)
          CALL INIINT (1+NELBLK+4, 0, IA(KLENF))
        END IF
        CALL MDRSRV ('NLNKF', KNLNKF, NELBLK)
        CALL MDRSRV ('LINKF', KLINKF, 0)
        CALL MDRSRV ('IF2EL', KIF2EL, 0)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

        CALL MSSURF (A, IA, A(KLENE), A(KNLNKE), A(KLINKE),
     &    A(KLENF), A(KNLNKF), KLINKF, KIF2EL, C(KNMLB))
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KNLNKF = 1
        KLENF = 1
        KLINKF = 1
        KIF2EL = 1
      END IF

C   --Group lines connecting nodes for efficiency

      IF (MESHOK) THEN
        CALL MDRSRV ('LENL', KLENL, 3 + NELBLK)
        CALL INIINT (3+NELBLK, 0, IA(KLENL))
        CALL MDRSRV ('LINSET', KLNSET, 0)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

        CALL MSLINS (A,
     &    A(KLENF), A(KNLNKF), A(KLINKF), A(KIF2EL),
     &    A(KLENL), KLNSET)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KLENL = 1
        KLNSET = 1
      END IF

C   --Calculate (and initialize) the deformed mesh limits

      RDMESH(1) = 0.0
      RDMESH(2) = 0.0
      RDMESH(3) = 0.0
      RDMESH(4) = 0.0

      IF (MESHOK) THEN

C      --Compute NPSURF nodes that determine the mesh limits

        NUMNPF = NUMNP

C      --QNPICK uses NPFIND to find NPSURF
        CALL MDRSRV ('NPSURF', KNPSUR, NUMNPF)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

        IF (.NOT. IS3DIM) THEN
          CALL MAKSU2 (A(KLENL), A(KLNSET), MSHBOR,
     &      .FALSE., IDUM, A(KNPSUR))
        ELSE
          CALL MAKSUR (A(KLENF), A(KNLNKF), A(KLINKF),
     &      .FALSE., IDUM, A(KNPSUR))
        END IF

        CALL MDLONG ('NPSURF', KNPSUR, NNPSUR)

C         --Find the displacement variables

        CALL FNDDIS (NAMECO, C(KNAMES+NAMLEN*(KNAMNV-1)),
     &    DEFOK, IXDEF, IYDEF, IZDEF, DEFFAC, NAMLEN)

C         --Find the undeformed mesh limits

        CALL MINMXS (NNPSUR, A(KNPSUR), A(KXN),
     &    UNMESH(KLFT), UNMESH(KRGT))
        CALL MINMXS (NNPSUR, A(KNPSUR), A(KYN),
     &    UNMESH(KBOT), UNMESH(KTOP))
        IF (IS3DIM) CALL MINMXS (NNPSUR, A(KNPSUR), A(KZN),
     &    UNMESH(KNEA), UNMESH(KFAR))

        CALL CPYREA (2*NDIM, UNMESH, ALMESH)

C         --Calculate the deformed mesh limits

        CALL DEFLIM (A, A(KWHOLE),
     &    A(KXN), A(KYN), A(KZN), A(KNPSUR))
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KNPSUR = 1
      END IF

C   --Reserve memory for all programs

      CALL MDGET (1+NUMNP + 1+NELBLK + 1+NUMEL
     &  + 1+NUMNPS + 1+NUMESS
     &  + 1+NVARHI + 1+NVARGL + 1+NVARNP + 1+NVAREL
     &  + 2+NVARHI+NVARGL)
      CALL MDRSRV ('LISNP', KLINP, 1+NUMNP)
      CALL MDRSRV ('NLISEL', KNLIEL, 1+NELBLK)
      CALL MDRSRV ('LISEL', KLIEL, 1+NUMEL)
      CALL MDRSRV ('LISNPS', KLINPS, 1+NUMNPS)
      CALL MDRSRV ('LISESS', KLIESS, 1+NUMESS)
      IF (EXODUS) THEN
        CALL MDRSRV ('LISHV', KLIHV, 1+NVARHI)
        CALL MDRSRV ('LISGV', KLIGV, 1+NVARGL)
        CALL MDRSRV ('LISNV', KLINV, 1+NVARNP)
        CALL MDRSRV ('LISEV', KLIEV, 1+NVAREL)
        CALL MDRSRV ('LIDSP', KLIDP, 2+NVARHI+NVARGL)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KLIHV = 1
        KLIGV = 1
        KLINV = 1
        KLIEV = 1
        KLIDP = 1
      END IF

      NPTIMS = 0
      NPTIMW = 0
      CALL MDRSRV ('IPTIMS', KPTIMS, MAX (NSTEPS, 1))
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

C   --Reserve memory for mesh plots

      CALL MDRSRV ('BLKCOL', KBKCOL, 1+NELBLK)
      CALL MDRSRV ('IELBST', KELBST, NELBLK)
      CALL MDRSRV ('ISSNPS', KSSNPS, NUMNPS*4)
      CALL MDRSRV ('ISSESS', KSSESS, NUMESS*4)
      CALL MDRSRV ('SHDCOL', KSHDCL, NELBLK*7)
      CALL MDRSRV ('ISHDCL', KISHCL, NELBLK*3)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160
      CALL INIREA (NELBLK*7, 0.0,  A(KSHDCL))
      CALL INIINT (NELBLK*3, 0,   IA(KISHCL))

C   --Reserve memory for DETOUR, if able to run

      CALL DTCHK (.FALSE., DTOK)

      IF (DTOK) THEN
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      END IF

C   --Reserve memory for PATHLINE, if able to run

      CALL LNCHK (.FALSE., LNOK)

      IF (LNOK) THEN
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      END IF

C   --Reserve memory for TPLOT, if able to run

      CALL TPCHK (.FALSE., TPOK)

      IF (TPOK) THEN
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      END IF

C   --Reserve memory for SPLOT, if able to run

      CALL SPCHK (.FALSE., SPOK)

      IF (SPOK) THEN
        MAXNE = MAX (NUMNP, NUMEL)
        CALL MDRSRV ('NENUM', KNENUM, MAXNE)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KNENUM = 1
      END IF

      IF (LNOK .OR. SPOK) THEN
        CALL MDRSRV ('IPATH', KIPATH, 0)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160
      ELSE
        KIPATH = 1
      END IF

C     Initialize array containing list of display variables
      IF (EXODUS) CALL DISPV (.TRUE., ' ', IDUM, IDUM,
     &  ' ', C(KNAMES), A(KLIDP), NAMLEN)

C     Initialize BLKCOL array.
      if (meshok) then
        CALL BCOLOR (.TRUE., ' ', IDUM, IDUM, IDUM,' ', A(KBKCOL))
      end if

C     Initialize line thicknesses for mesh plots
      CALL LINTHK (CDUM, IDUM, IDUM, IDUM, RDUM, CDUM, .TRUE.)

        write (*,9999)
 9999   FORMAT(/,
     *    10x,'NOTE: This version of blot uses global ids for both',
     *    ' node and element ids by default.',/,
     *    10x,'      To see the mapping from local to global, use',
     *    ' the commands:',/,
     *    10x,'          "LIST MAP" (element map), or ',
     *    '"LIST NODEMAP" (node map)',/,
     *    10x,'      To disable the maps and use local ids, restart',
     *    ' blot with "-nomap node|element|all"',//,
     *    10x,'      Notify gdsjaar@sandia.gov if bugs found')

        if (mapel .and. mapnd) then
          WRITE (*, 10010) 'Nodes and Elements using Global Ids'
        else if (mapel) then
          WRITE (*, 10010) 'Elements use Global Ids, Node Ids are Local'
        else if (mapnd) then
          WRITE (*, 10010) 'Element use Local Ids, Node Ids are Global'
        else
          WRITE (*, 10010) 'Nodes and Elements using Local Ids'
        end if

 130  CONTINUE
      IF (.TRUE.) THEN
        CALL MDLONG ('IPTIMS', KPTIMS, MAX (NSTEPS, 1))
        IF (SPOK) THEN
          CALL MDLONG ('NENUM', KNENUM, MAXNE)
        END IF
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

        CALL COMAND (A, CURPRO, C(KQAREC), C(KINFO),
     &    NAMECO, C(KNMLB), C(KNAMES), A(KTIMES), A(KWHOLE),
     *    A(KPTIMS), A(KMAPEL), A(KMAPND),
     &    A(KIDELB), NEWELB, A(KELBST), A(KE2ELB),
     &    A(KLENE), A(KNLNKE), A(KLINKE),
     &    A(KXN), A(KYN), A(KZN), A(KXE), A(KYE), A(KZE),
     &    A(KIEVOK),
     &    A(KSSNPS), A(KIDNS), A(KNNNS),
     &    A(KSSESS), A(KIDSS), A(KNESS), A(KNNSS),
     &    NCSTEP, A(KLINP), A(KNLIEL), A(KLIEL),
     &    A(KLINPS), A(KLIESS),
     &    A(KLIHV), A(KLIGV), A(KLINV), A(KLIEV), A(KLIDP),
     &    A(KBKCOL), A(KNENUM), NEUTRL, NEWPRO, A(KSHDCL),
     *    A(KISHCL), C(KNMEB), C(KNMNS), C(KNMSS), NAMLEN, *150)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

C         --Calculate the undeformed and deformed mesh limits

        IF ((CURPRO .EQ. 'DETOUR') .AND. (DEFFAC .LT. 0.0)) THEN
          CALL DEFLIM (A, A(KWHOLE),
     &      A(KXN), A(KYN), A(KZN), A(KNPSUR))
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 160
        END IF

        CALL MDLONG ('IPTIMS', KPTIMS, NPTIMS)
        IF (SPOK) THEN
          CALL MDLONG ('NENUM', KNENUM, NNENUM)
        END IF
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 160

        IF (NEWPRO .EQ. 'P') THEN
          CALL PLOTER (A, CURPRO, NEUTRL,
     &      C(KNAMES), NPTIMS, A(KPTIMS), A(KTIMES), A(KWHOLE),
     &      A(KNLNKE), A(KIEVOK), A(KNENUM),
     &      A(KXN), A(KYN), A(KZN), A(KXE), A(KYE), A(KZE),
     &      NAMECO, A(KLENF), A(KNLNKF), KLINKF,
     &      A(KLENL), KLNSET,
     &      A(KE2ELB), NEWELB, A(KELBST), KNPSUR,
     &      A(KSSNPS), A(KIDNS), A(KSSESS), A(KIDSS),
     &      A(KLIDP), A(KBKCOL), A(KIDELB),C(KNMLB), NAMLEN,
     *      A(KMAPEL), A(KMAPND))
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 160
        END IF

        GOTO 130
      END IF

 150  CONTINUE

C   --Finish graphics
      CALL GREXIT

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 160

      GOTO 180

 160  CONTINUE
      WRITE (*, *) 'Total memory used = ', MEM, ' words.'
      CALL MEMERR
      GOTO 180

 170  CONTINUE
      GOTO 180

 180  CONTINUE

C   --Close files

      CALL EXCLOS(NDB, IERR)
      IF (NEUOPN) THEN
        WRITE (NEU, 10030) 'xmax', wxmax
        WRITE (NEU, 10030) 'xmin', wxmin
        WRITE (NEU, 10030) 'ymax', wymax
        WRITE (NEU, 10030) 'ymin', wymin
        CLOSE (NEU, IOSTAT=IDUM)
      END IF
      IF (ANYPRT) CLOSE (NPRT, IOSTAT=IDUM)
      IF (GRFOPN) CLOSE(NEUGRF)
      IF (CSVOPN) CLOSE(NCSV)

 190  continue
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

10010 FORMAT (/, 1X, 5A)
10020 FORMAT (/
     &  21X,'BBBBBBBB    LL            OOOOOO    TTTTTTTTTT', /
     &  20X,'BBBBBBBBB   LL           OOOOOOOO   TTTTTTTTTT', /
     &  19X,'BB      BB  LL          OO      OO      TT    ', /
     &  18X,'BB      BB  LL          OO      OO      TT    ', /
     &  17X,'BBBBBBBBB   LL          OO      OO      TT    ', /
     &  16X,'BBBBBBBBB   LL          OO      OO      TT    ', /
     &  15X,'BB      BB  LL          OO      OO      TT    ', /
     &  14X,'BB      BB  LL          OO      OO      TT    ', /
     &  13X,'BBBBBBBBB   LLLLLLLLLL   OOOOOOOO       TT    ', /
     &  12X,'BBBBBBBB    LLLLLLLLLL    OOOOOO        TT    II-2')
10030 FORMAT ('@ world ',A4,1x,1pe15.7E3)

      END

      subroutine inimap(num, iar)
      dimension iar(*)
      do 10 i=1, num
        iar(i) = i
 10   continue
      return
      end

C=======================================================================
      SUBROUTINE CPYST8 (LEN, IFROM, ITO)
C=======================================================================
C   --Parameters:
C   --   LEN   - IN  - the number of strings in the list
C   --   IFROM - IN  - the input list
C   --   ITO   - OUT - the copied list
C   --   ISTLN - IN  - string length

      INTEGER LEN
      include 'params.blk'
      CHARACTER*(MXSTLN) IFROM(*), ITO(*)

      DO 100 I = 1, LEN
         ITO(I) = IFROM(I)
  100 CONTINUE

      RETURN
      END

C ... Wrapper to get strings the right length
      subroutine exgqaw(ndb, qarec, ierr)
      include 'params.blk'
      character*(mxstln) qarec(4, *)
      call exgqa(ndb, qarec, ierr)
      return
      end
      subroutine exginw(ndb, info, ierr)
      include 'params.blk'
      character*(mxlnln) info(*)
      call exginf(ndb, info, ierr)
      return
      end

C=======================================================================
      SUBROUTINE RDTYPE (NELBLK,NAMELB,IDELB,NUMELB,ISHEX,HEXID,NSHL)
C=======================================================================

      INTEGER NELBLK, ISHEX
      INTEGER IDELB(*), HEXID(*), numelb(*)
      CHARACTER*32 NAMELB(*)

C     Initialize HEX return variables
      ISHEX = 0
      NSHL = 0
C     Mark HEXID array with corresponding element blocks that have HEXSHELLs
      DO 310 I = 1, NELBLK
C        Check if element block contains HEXSHELLs
         IF (NAMELB(I)(1:8) .EQ. 'HEXSHELL') THEN
C           ISHEX - counter: how many HEXSHELL element blocks
C           index - HEXID(ISHEX) - stores HEXSHELL element block id
            ISHEX = ISHEX + 1
            HEXID(ISHEX) = IDELB(I)
         ELSE IF (NAMELB(I)(1:5) .EQ. 'SHELL') THEN
            NSHL = NSHL + NUMELB(I)
         ENDIF
 310  CONTINUE

      RETURN
      END

C=======================================================================
      SUBROUTINE GETNAMS (NDB, ITYPE, NUM, NAMES, NAMLEN)
C=======================================================================

      CHARACTER*(NAMLEN) NAMES(*)
      CALL EXGNAMS(NDB, ITYPE, NUM, names, ierr)
      RETURN
      END

      subroutine iniseq(icnt, map)
      integer map(*)
      do i=1, icnt
        map(i) = i
      end do
      return
      end
