C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.

C See packages/seacas/LICENSE for details

      PROGRAM NUMBER

C ... Program to calculate the centroid location and the
C     mass moment of inertia for a axisymmetric and plane 2-D mesh
C     and a 3-D mesh defined in the genesis format.

C     LINK WITH SUPES LIBRARY
      include 'exodusII.inc'

      include 'nu_progqa.blk'
      include 'nu_numg.blk'
      include 'nu_varcnt.blk'
      include 'nu_mass.blk'
      include 'nu_cvty.blk'
      include 'nu_logs.blk'
      include 'nu_ptim.blk'
      include 'nu_nset.blk'
      include 'nu_io.blk'
      include 'nu_ndisp.blk'
      include 'argparse.inc'

      CHARACTER*2048 DBNAME, SCRATCH

      CHARACTER*(MXLNLN)  TITLE
      DIMENSION A(1), IA(1)
      EQUIVALENCE (A(1),IA(1))
      CHARACTER*1 C(1)

      PARAMETER (MXNAM = 256)
      CHARACTER*(MXSTLN) NAMECO(6), NAMES(MXNAM)

      integer cmpsiz, iows

      include 'nu_qainfo.blk'

C        UNIT 6 = STANDARD OUTPUT
C             7 = ASCII OUTPUT
C             9 = BINARY MESH INPUT (GENESIS)

      ITERM = 6
      IHARD = 7
      NDB   = 0

      CALL STRTUP (QAINFO)
      CALL BANNER (ITERM, QAINFO,
     &  'A GENESIS/EXODUS DATABASE INFORMATION PROGRAM',
     &  ' ', ' ')
      CALL BANNER (IHARD, QAINFO,
     &  'A GENESIS/EXODUS DATABASE INFORMATION PROGRAM',
     &  ' ', ' ')
      call cpyrgt (ITERM, '1988')
      call cpyrgt (IHARD, '1988')

C ... GET FILENAMES:

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()

      if (narg .lt. 1) then
        CALL PRTERR ('FATAL', 'Filename(s) not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "numbers filename [output]"')
        GOTO 60
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "numbers filename [output]"')
        GOTO 60
      end if

      if (narg .eq. 2) then
        CALL get_argument(2,dbname, lfil)
      else
        dbname = "numbers.o"
        lfil=lenstr(dbname)
      end if

      open(unit=ihard, file=dbname(:lfil), iostat=ierr)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Could not create "'//dbname(:LFIL)//'"'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 60
      END IF

      CMPSIZ = 0
      IOWS   = 0
      DBNAME  = ' '
      CALL get_argument(1,DBNAME, LNAM)
      NDB = exopen(dbname(:lnam), EXREAD, CMPSIZ, IOWS, vers, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//dbname(:lnam)//'" does not exist.'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 60
      END IF

      CALL MDINIT (A)
      CALL MCINIT (C)

      call exgini(ndb, title, ndim, numnp, numel, nelblk,
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

      CALL DBPINI ('TIS', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     *  NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF,
     *  0, 0, 0, DBNAME(:LNAM))

      AXI    = .TRUE.
      EXODUS = .FALSE.
      NNODES = 2**NDIM

      CALL MDRSRV ('CRD', IR, NUMNP*NDIM)
      IRX = IR
      IRY = IR + NUMNP
      IRZ = IR + 2 * NUMNP

      CALL MDRSRV ('MAT',  IM, 6*NELBLK)
      CALL MDRSRV ('LINK', IX, 0)

      CALL MDRSRV ('DENS',ID, NELBLK)
      CALL MDRSRV ('WAVE',IW, NELBLK)
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) GO TO 55

      CALL INIREA ( NELBLK, 0.0, A(ID))
      CALL INIREA ( NELBLK, 0.0, A(IW))

      call exgcor (ndb, a(irx), a(iry), a(irz), ierr)
      if (ierr .ne. 0) go to 60
      call exgcon (ndb, nameco, ierr)
      if (ierr .ne. 0) go to 60
      do 10 i=1, ndim
        call exupcs(nameco(i))
 10   continue

C ... Scratch space for block info

      CALL MDRSRV ('IDELB',  IDELB,  NELBLK)
      CALL MDRSRV ('NUMELB', NUMELB, NELBLK)
      CALL MDRSRV ('NUMLNK', NUMLNK, NELBLK)
      CALL MDRSRV ('NUMATR', NUMATR, NELBLK)
      CALL MCRSRV ('NAMELB', KNAMEL, MXSTLN*NELBLK)
      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) GO TO 55

      CALL DBIELB (NDB, 'HIC', 1, NELBLK, IA(IDELB), IA(NUMELB),
     *  IA(NUMLNK), IA(NUMATR), A, IA, IX, IDUM, C(KNAMEL), *60)
      CALL TRBLK (IA(IDELB), IA(NUMELB), IA(NUMLNK), IA(IM),
     &  NELBLK, NNODES)
      CALL SORBLK (IA(IDELB), IA(NUMLNK), IA(IM), NELBLK)
      CALL MDDEL ('IDELB' )
      CALL MDDEL ('NUMLNK')
      CALL MDDEL ('NUMELB')
      CALL MDDEL ('NUMATR')

C ... BOUNDARY CONDITION FLAGS

C -- Node Sets:
C    INS1 = IDNPS  (NUMNPS) NODAL POINT SET IDS
C    INS2 = NNNPS  (NUMNPS) NODAL POINT SET COUNTS
C    INS3 = IPTNPS (NUMNPS) NODAL POINT SET POINTER
C    INS4 = LSTNPS (LNPSNL) NODAL POINT SET NODE LIST
C    INS5 = FACNPS (LNPSNL) NODAL POINT DISTRIBUTION FACTORS

C -- Element Side Sets:
C    IBC1 = IDESS  (NUMESS) ELEMENT SIDE SET IDS
C    IBC2 = NEESS  (NUMESS) ELEMENT SIDE SET ELEMENT COUNTS
C    IBC3 = NNESS  (NUMESS) ELEMENT SIDE SET NODE    COUNTS
C    IBC4 = IPEESS (NUMESS) ELEMENT SIDE SET ELEMENT POINTERS
C    IBC5 = IPNESS (NUMESS) ELEMENT SIDE SET NODE    POINTERS
C    IBC6 = LTEESS (LESSEL) ELEMENT SIDE SET ELEMENT LIST
C    IBC7 = LTNESS (LESSNL) ELEMENT SIDE SET NODE    LIST
C    IBC8 = FACESS (LESSNL) ELEMENT SIDE SET DISTRIBUTION FACTORS

      CALL MDRSRV ('IDNPS',  INS1, NUMNPS)
      CALL MDRSRV ('NNNPS',  INS2, NUMNPS)
      CALL MDRSRV ('NDNPS',  INS6, NUMNPS)
      CALL MDRSRV ('IPTNPS', INS3, NUMNPS)
      CALL MDRSRV ('IPTNDS', INS7, NUMNPS)
      CALL MDRSRV ('LSTNPS', INS4, LNPSNL)
      CALL MDRSRV ('FACNPS', INS5, LNPSNL)

      CALL MDRSRV ('IDESS',  IBC1, NUMESS)
      CALL MDRSRV ('NEESS',  IBC2, NUMESS)
      CALL MDRSRV ('NNESS',  IBC3, NUMESS)
      CALL MDRSRV ('NDESS',  IBC11,NUMESS)
      CALL MDRSRV ('IPEESS', IBC4, NUMESS)
      CALL MDRSRV ('IPNESS', IBC5, NUMESS)
      CALL MDRSRV ('IXDESS', IBC9, NUMESS) !
      CALL MDRSRV ('LTNNSS', KLTNNN, LESSEL)
      CALL MDRSRV ('LTEESS', IBC6, LESSEL)
      CALL MDRSRV ('LTNESS', IBC7, LESSNL)
      CALL MDRSRV ('LTSESS', IBC10,LESSEL)
      CALL MDRSRV ('FACESS', IBC8, LESSNL)

      CALL MDSTAT (NERRS, NUSED)
      IF (NERRS .GT. 0) GO TO 55

      if (numnps .gt. 0) then
         call exgcns(ndb, a(ins1), a(ins2), a(ins6), a(ins3), a(ins7),
     &        a(ins4), a(ins5), ierr)
         if (ierr .ne. 0) go to 60
      end if

      if (numess .gt. 0) then
        call exgcss(ndb, a(ibc1), a(ibc2), a(ibc11), a(ibc4),
     *    a(ibc9), A(ibc6), a(ibc10), a(ibc8), ierr)
c     ... Now convert sides to nodes....
C     ... This code stolen from ex2ex1v2, Vic Yarberry
C     offset into element list for current side set
        isoff = 0
C     node count for current side set
        nodcnt = 0
        do 104 i=0,numess-1
C     update index array
          ia(ibc5+i)=nodcnt+1
C     get num of sides & df
          call exgsp(ndb,ia(ibc1+i),nsess,ndess,nerr)

C     get side set nodes
          if (nsess .gt. 0) then
            call exgssn(ndb,ia(ibc1+i),ia(kltnnn+isoff),
     &        ia(ibc7+nodcnt),nerr)
            if (nerr .gt. 0) goto 60
          end if
          nness = 0
C     sum node counts to calculate next index
          do 102 ii=0,nsess-1
            nness=nness+ia(kltnnn+isoff+ii)
 102      continue
          ia(ibc3+i)=nness
          nodcnt=nodcnt+nness
          isoff=isoff+nsess
 104    continue

        if (ierr .ne. 0) go to 60
      end if

C ... TRY TO READ QA RECORDS.  IF EOF THEN NOT EXODUS FORMAT

      call exinq(ndb, EXQA,   nqarec, rdum, cdum, ierr)
      call exinq(ndb, EXINFO, ninfo,  rdum, cdum, ierr)
      call mcrsrv('QAREC', kqarec, nqarec * 4 * MXSTLN)
      call mcrsrv('INFREC', kinfo, ninfo * MXLNLN)
      CALL MCSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 60
      if (nqarec .gt. 0) then
C     ... Wrapper to get strings the right length
        call exgqaw(ndb, c(kqarec), ierr)
      end if
      if (ninfo .gt. 0) then
C     ... Wrapper to get info record the right length
        call exginw(ndb, c(kinfo), ierr)
      end if

      CALL INISTR (MXNAM, ' ', NAMES)

      call exgvp(ndb, 'G', nvargl, ierr)
      call exgvp(ndb, 'N', nvarnp, ierr)
      call exgvp(ndb, 'E', nvarel, ierr)
      EXODUS = (nvargl + nvarnp + nvarel) .gt. 0

      ixgv = 1
      ixnv = ixgv + nvargl
      ixev = ixnv + nvarnp

      if (nvargl .gt. 0) then
        call exgvan(ndb, 'G', nvargl, names(ixgv), ierr)
      end if
      if (nvarnp .gt. 0) then
        call exgvan(ndb, 'N', nvarnp, names(ixnv), ierr)
      end if
      if (nvarel .gt. 0) then
        call exgvan(ndb, 'E', nvarel, names(ixev), ierr)
      end if
      do 40 i=1, ixev+nvarel
        call exupcs(names(i))
 40   continue

      IF (EXODUS) THEN
C ... Read truth table
        call mdrsrv ('ISEVOK', IISEV, nvarel*nelblk)
        CALL MDSTAT (NERRS, NUSED)
        IF (NERRS .GT. 0) GO TO 55
        call exgvtt (ndb, nelblk, nvarel, ia(iisev), ierr)
        CALL DBPNAM ('*', NVARGL, NVARNP, NVAREL,
     *    NAMES(IXGV), NAMES(IXNV), NAMES(IXEV))

C     --Read the database time steps
c     determine how many time steps are stored
        call exinq (ndb, EXTIMS, NSTEP, rdum, cdum, ierr)
        call mdrsrv('TIMES', ITIME, max(NSTEP,1))
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 55

c     read time values at all time steps
        if (nstep .gt. 0) then
          call exgatm (ndb, a(itime), ierr)
        end if

        CALL DBPTIM ('NM', NSTEP, A(ITIME))
        CALL FNDDIS (NAMECO, NAMES(IXNV), ISDIS, NDIM, NVARNP,
     *    NDISP(1), NDISP(2), NDISP(3))

C ...    DONE READING, START CALCULATING

        CALL MDRSRV ('ITMSEL', ITSEL, max(1,NSTEP))
        CALL MDRSRV ('DISPS',  IDSP, NDIM * NUMNP)
      ELSE
        ITSEL = 1
        IDSP = 1
      END IF

C ... CALCULATE ELEMENT CENTROIDS FOR LATER USE

      CALL MDRSRV ('ELCEN',  IECEN, NDIM*NUMEL)
      CALL MDSTAT (MNERRS, MNUSED)
      IF (NERRS .GT. 0) GO TO 55

      CALL ELCENT ( A(IECEN), A(IX), A(IR), NDIM, NUMEL, NNODES, NUMNP)
      CALL HEADER (NDIM, TITLE, NUMEL, NUMNP, AXI, DBNAME(:LNAM))
      CALL COMMAND (A, IA, TITLE, A(ITIME), A(ITSEL), A(IM),
     *  A(IDSP), A(IR), A(IX), A(ID), A(IW), A(IISEV),
     *  NAMES(IXGV), NAMES(IXNV), NAMES(IXEV),
     *  NQAREC, C(KQAREC), NINFO, C(KINFO), DBNAME(:LNAM))
      GO TO 60

 55   CONTINUE
      CALL MEMERR
      GO TO 60

 60   CONTINUE
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))
      if (ndb .gt. 0) then
        call exclos(ndb, ierr)
      end if
      STOP
      END

      subroutine exgqaw(ndb, qarec, ierr)
      include 'exodusII.inc'
      character*(mxstln) qarec(4, *)
      call exgqa(ndb, qarec, ierr)
      return
      end

      subroutine exginw(ndb, info, ierr)
      include 'exodusII.inc'
      character*(mxlnln) info(*)
      call exginf(ndb, info, ierr)
      return
      end
