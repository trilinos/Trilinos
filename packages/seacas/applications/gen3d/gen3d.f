C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      PROGRAM GEN3D
C=======================================================================

C   --*** GEN3D *** (GEN3D) GENESIS 2D to 3D Program
C   --
C   --GEN3D inputs a 2D GENESIS database and outputs a 3D GENESIS database.
C   --The input mesh is either translated along the Z coordinate or rotated
C   --around a Y axis.  The user specifies the number of translations or
C   --rotations.  For rotations, the total number of degrees the 2D mesh
C   --is to be rotated is also user-specified.  The mesh may be rotated
C   --around an edge of the input mesh (creating part of a cylinder-shaped
C   --mesh) or around a line outside the input mesh (creating part of a
C   --cylinder-shaped mesh with a hole in a middle).
C   --
C   --Expected input:
C   --   o The commands on the standard input device.
C   --   o The 2D GENESIS database on unit 9
C   --     (must have 4 nodes per element,
C   --
C   --Output:
C   --   o A listing of the input database information and any errors
C   --     found on the standard output device.
C   --   o The 3D GENESIS database on unit 10.

C   --Developed at Sandia National Laboratories.
C   --
C   --Current author and code sponsor: Gregory D. Sjaardema
C   --
C   --Revision History:
C   --   04/86 Created (Amy Gilkey)
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C   --
C   --Runs on Unix systems

C   --Documentation:
C   --   "User's Manual for GEN3D"

      include 'exodusII.inc'

      INCLUDE 'g3_progqa.blk'
      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbtitl.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'
      INCLUDE 'g3_params.blk'
      INCLUDE 'g3_xyzoff.blk'
      INCLUDE 'g3_xyzrot.blk'
      INCLUDE 'g3_xyzmir.blk'
      INCLUDE 'g3_twist.blk'
      INCLUDE 'argparse.inc'

      CHARACTER*2048 FILIN, FILOUT, SCRATCH

      CHARACTER*(MXSTLN) NAMECO(6)
C      --NAMECO - the coordinate names

C... String containing name of common element topology in model
C    or 'MULTIPLE_TOPOLOGIES' if not common topology.
      character*(MXSTLN) comtop

      CHARACTER CDUM
      logical l64bit

      INTEGER CMPSIZ, IOWS

      DIMENSION A(1)
      INTEGER IA(1)
      EQUIVALENCE (A(1), IA(1))
      CHARACTER*1 C(1)
C      --A - the dynamic numeric memory base array

      INTEGER IDNSET(0:MAXSET,2)
      INTEGER IDESET(0:MAXSET,2)

      INCLUDE 'g3_qainfo.blk'

      CALL STRTUP (QAINFO)

      WRITE (*, 70)
      WRITE (*, 80)
      CALL BANNER (0, QAINFO,
     &   'AN EXODUSII DATABASE 2D TO 3D CONVERSION PROGRAM',
     &   ' ', ' ')
      call cpyrgt (0, '1989')

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C .. Get filename from command line.  If not specified, emit error message
      l64bit = .false.
      NARG = argument_count()
      iarg = 1

      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "gen3d [-64] 2dfilename 3dfilename"')
        GOTO 60
      else if (narg .eq. 3) then
        CALL get_argument(iarg,FILIN, LNAM)
        if (filin(:lnam) .eq. '-64') then
          l64bit = .true.
        else
          SCRATCH = 'Unrecognized command option "'//FILIN(:LNAM)//'"'
          CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
          CALL PRTERR ('FATAL',
     *      'Syntax is: "gen3d [-64] 2dfilename 3dfilename"')
          GOTO 60
        end if
        iarg = 2;
      else if (narg .gt. 3) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *      'Syntax is: "gen3d [-64] 2dfilename 3dfilename"')
        GOTO 60
      end if

C   --Open the input database and read the initial variables

      NDBIN  = 9
      NDBOUT = 10

      CMPSIZ = 0
      IOWS   = 0

      FILIN = ' '
      CALL get_argument(iarg,FILIN, LNAM)
      NDBIN = exopen(filin(:lnam), EXREAD, CMPSIZ, IOWS, vers, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//FILIN(:LNAM)//'" does not exist.'
         CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
         GOTO 60
      END IF

      call exgini(ndbin, title, ndim, numnp, numel, nelblk,
     *     numnps, numess, ierr)
      if (numnps .gt. 0) then
         call exinq(ndbin, EXNSNL, lnpsnl, rdum, cdum, ierr)
         call exinq(ndbin, EXNSDF, lnpsdf, rdum, cdum, ierr)
      else
         lnpsnl = 0
         lnpsdf = 0
      end if
      if (numess .gt. 0) then
         call exinq(ndbin, EXSSNL, lessnl, rdum, cdum, ierr)
         call exinq(ndbin, EXSSEL, lessel, rdum, cdum, ierr)
         call exinq(ndbin, EXSSDF, lessdf, rdum, cdum, ierr)
      else
         lessnl = 0
         lessel = 0
         lessdf = 0
      end if

      CALL DBPINI ('NTIS', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &     NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL,
     &     LESSDF, IDUM, IDUM, IDUM, FILIN)

      IF (NDIM .NE. 2) THEN
         CALL PRTERR ('FATAL', 'Number of dimensions must be 2')
         GOTO 60
      END IF

C   --Reserve memory for the 2D information

      CALL MDRSRV ('XN', KXN, NUMNP)
      CALL MDRSRV ('YN', KYN, NUMNP)
      KZN = 1

      CALL MDRSRV ('IDELB',  KIDELB, NELBLK)
      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
      CALL MDRSRV ('NUMLNK', KNLNK, NELBLK)
      CALL MDRSRV ('NUMATR', KNATR, NELBLK)
      CALL MCRSRV ('NAMELB', KNMLB, MXSTLN*NELBLK)
      CALL MCRSRV ('BLKTYP', KBKTYP, NELBLK)

      CALL MDRSRV ('LINK', KLINK, 4 * NUMEL)
      CALL MDRSRV ('IDNPS',  KIDNS, NUMNPS)
      CALL MDRSRV ('NNNPS',  KNNNS, NUMNPS)
      CALL MDRSRV ('NDNPS',  KNDNPS, NUMNPS)
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
      CALL MDRSRV ('IXDNPS', KIXDNS, NUMNPS)
      CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
      CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)

      CALL MDRSRV ('IDESS',  KIDSS, NUMESS)
      CALL MDRSRV ('NEESS',  KNESS, NUMESS)
      CALL MDRSRV ('NNESS',  KNNSS, NUMESS)
      CALL MDRSRV ('NDESS',  KNDSS,  NUMESS)
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
      CALL MDRSRV ('IXNESS', KIXNSS, NUMESS)
      CALL MDRSRV ('IXDESS', KIXDSS, NUMESS)
      CALL MDRSRV ('LTNNSS', KLTNNN, LESSEL)
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
      CALL MDRSRV ('LTNESS', KLTNSS, LESSNL)
      CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
      CALL MDRSRV ('FACESS', KFACSS, LESSNL)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Read 2D information from the database and close file

      call exgcor (ndbin, a(kxn), a(kyn), a(kzn), ierr)

C ... Don't warn about no map stored in file
      call exopts (0, ierr)
      call exopts (EXVRBS, ierr)

      CALL INISTR (NDIM, ' ', NAMECO)
      call exgcon (ndbin, nameco, ierr)

      CALL RDELB (A, IA(KIDELB), C(KNMLB), IA(KNELB), IA(KNLNK),
     &  IA(KNATR), IA(KLINK), KATRIB, *60)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      call CHKTOP(NELBLK, C(KNMLB), COMTOP)

      if (numnps .gt. 0) then
         call exgcns(ndbin, ia(kidns), ia(knnns), ia(kndnps),
     &        ia(kixnns), ia(kixdns), ia(kltnns), a(kfacns), ierr)
         if (lnpsdf .eq. 0) then
           call inirea(lnpsnl, 1.0, a(kfacns))
         end if
      end if
      if (numess .gt. 0) then
         call exgcss(ndbin, ia(kidss), ia(kness), ia(kndss),
     &        ia(kixess), ia(kixdss), IA(KLTeSS), ia(kltsss),
     &        a(kfacss), ierr)
         if (lessdf .eq. 0) then
           call inirea(lessnl, 1.0, a(kfacss))
         end if

c     ... Now convert sides to nodes.... ia(kltsss),
C     ... This code stolen from ex2ex1v2, Vic Yarberry
C     offset into element list for current side set
         isoff = 0
C     node count for current side set
         nodcnt = 0
         do 104 i=0,numess-1
C     update index array
           ia(kixnss+i)=nodcnt+1
C     get num of sides & df
           call exgsp(ndbin,ia(kidss+i),nsess,ndess,nerr)

C     get side set nodes
           call exgssn(ndbin,ia(kidss+i),ia(kltnnn+isoff),
     &          ia(kltnss+nodcnt),nerr)
           if (nerr .gt. 0) goto 40
           nness = 0
C     sum node counts to calculate next index
           do 102 ii=0,nsess-1
             nness=nness+ia(kltnnn+isoff+ii)
 102       continue
           ia(knnss+i)=nness
           nodcnt=nodcnt+nness
           isoff=isoff+nsess
 104     continue

         if (ierr .ne. 0) go to 40
      end if

      call exinq(ndbin, EXQA,   nqarec, rdum, cdum, ierr)
      call exinq(ndbin, EXINFO, ninfo,  rdum, cdum, ierr)
      call mcrsrv('QAREC', kqarec, (nqarec+1) * 4 * MXSTLN)
      call mcrsrv('INFREC', kinfo, (ninfo+1) * MXLNLN)
      CALL MCSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40
      if (nqarec .gt. 0) then
C     ... Wrapper to get strings the right length
         call exgqaw(ndbin, c(kqarec), ierr)
      end if
      if (ninfo .gt. 0) then
C     ... Wrapper to get info record the right length
         call exginw(ndbin, c(kinfo), ierr)
      end if

      IF ((NQAREC .GT. 0) .OR. (NINFO .GT. 0)) THEN
         CALL DBPQA ('*', NQAREC,  c(kqarec), NINFO, c(kinfo))
      END IF

C   --Read in runtime parameters

      CALL MDRSRV ('IBPARM', KIBPAR, 4 * NELBLK)
C     -- This assumes only one attribute per element block
      CALL MDRSRV ('ELATTR', KELATT, NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL COMAND (IA(KIDNS), IA(KIDSS), IDNSET, IDESET,
     &   C(KBKTYP), IA(KIBPAR), IA(KIDELB), IA(KNELB), IA(KNLNK),
     &   C(KNMLB), A(KELATT), A(KXN), A(KYN), A, *60)

C   --Get the new numbers for the elements and nodes

      CALL MDRSRV ('IXEL', KIXEL, NUMEL)
      CALL MDRSRV ('INCEL', KINCEL, NUMEL)
      CALL MDRSRV ('NREL', KNREL, NUMEL)
      CALL MDRSRV ('IELCOL', KIECOL, NUMEL)
      CALL MDRSRV ('IXNP', KIXNP, NUMNP)
      CALL MDRSRV ('NRNP', KNRNP, NUMNP)
      CALL MDRSRV ('NPCEN', KNPCEN, NUMCOL * 2 * NUMEL)
      CALL INIINT (NUMCOL*2*NUMEL, 0, IA(KNPCEN))
      CALL MDRSRV ('IELROW', KELROW, NUMEL)
      CALL MDRSRV ('IROT', KROT, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL RENUMB (A, C(KBKTYP), IA(KNELB), IA(KLINK), A(KXN), A(KYN),
     &     IA(KIXEL), IA(KINCEL), IA(KNREL), IA(KIECOL), IA(KIXNP),
     &     IA(KNRNP), IA(KNPCEN), IA(KELROW), IA(KROT))

      CALL MDDEL ('IELROW')
      CALL MDLONG ('NPCEN', KNPCEN, NUMCOL * NUMROW)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Get the node sets

      LNPSNO = INTADD (NUMNPS, IA(KNNNS)) * NNREPL
      CALL MDRSRV ('NNNP3', KNNN3, NUMNPS)
      CALL MDRSRV ('IXNNP3', KIXNN3, NUMNPS)
      CALL MDRSRV ('LTNNP3', KLTNN3, LNPSNO)
      CALL MDRSRV ('FACNP3', KFACN3, LNPSNO)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL NEWNPS (IDNSET(0,1), IDNSET(0,2),
     &   IA(KIDNS), IA(KNNNS), IA(KNNN3), IA(KIXNNS), IA(KIXNN3),
     &   IA(KLTNNS), IA(KLTNN3), A(KFACNS), A(KFACN3),
     &   IA(KIXNP), IA(KNRNP))

      CALL MDDEL ('NNNPS')
      CALL MDDEL ('IXNNPS')
      CALL MDDEL ('LTNNPS')
      CALL MDDEL ('FACNPS')
      CALL MDLONG ('LTNNP3', KLTNN3, LNPSNO)
      CALL MDLONG ('FACNP3', KFACN3, LNPSNO)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Get the side sets, and the front and back side sets

      NSSET = IDESET(0,1)+IDESET(0,2)
      IF (NSSET .GT. 0) THEN
         CALL MDRSRV ('ISSFRO', KISFRO, NUMEL)
         CALL MDRSRV ('ISSBCK', KISBCK, NUMEL)
         CALL MDRSRV ('NSSFRO', KNSFRO, 4*NUMEL)
         CALL MDRSRV ('NSSBCK', KNSBCK, 4*NUMEL)
      ELSE
         KISFRO = 1
         KISBCK = 1
         KNSFRO = 1
         KNSBCK = 1
      END IF

      LESSEO = INTADD (NUMESS, IA(KNESS)) * NEREPL
      LESSNO = 4 * LESSEO
      CALL MDRSRV ('NEES3', KNES3, NUMESS)
      CALL MDRSRV ('NNES3', KNNS3, NUMESS)
      CALL MDRSRV ('IXEES3', KIXES3, NUMESS)
      CALL MDRSRV ('IXNES3', KIXNS3, NUMESS)
      CALL MDRSRV ('LTEES3', KLTES3, LESSEO)
      CALL MDRSRV ('LTSES3', KLTSS3, LESSEO)
      CALL MDRSRV ('LTNES3', KLTNS3, LESSNO)
      CALL MDRSRV ('FACES3', KFACS3, LESSNO)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL NEWESS
     &     (IDESET(0,1), IDESET(0,2),
     &     IA(KLINK), IA(KISFRO), IA(KISBCK), NSSUR, NESUR,
     &     IA(KNSFRO), IA(KNSBCK),
     &     IA(KIDSS), IA(KNESS), IA(KNES3), IA(KNNSS), IA(KNNS3),
     &     IA(KIXESS), IA(KIXES3), IA(KIXNSS), IA(KIXNS3),
     &     IA(KLTESS), IA(KLTES3), IA(KLTSSS), IA(KLTSS3),
     &     IA(KLTNSS), IA(KLTNS3), A(KFACSS), A(KFACS3),
     &     IA(KIXEL), IA(KINCEL), IA(KNREL), IA(KIECOL), IA(KIXNP),
     &     IA(KNRNP), IA(KROT))

      CALL MDDEL ('NEESS')
      CALL MDDEL ('NNESS')
      CALL MDDEL ('IXEESS')
      CALL MDDEL ('IXNESS')
      CALL MDDEL ('LTEESS')
      CALL MDDEL ('LTNESS')
      CALL MDDEL ('FACESS')
      CALL MDDEL ('IROT')
      CALL MDLONG ('LTEES3', KLTES3, LESSEO)
      CALL MDLONG ('LTNES3', KLTNS3, LESSNO)
      CALL MDLONG ('FACES3', KFACS3, LESSNO)
      IF (NSSET .GT. 0) THEN
         CALL MDLONG ('NSSFRO', KNSFRO, NSSUR)
         CALL MDLONG ('NSSBCK', KNSBCK, NSSUR)
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Open the output database

      FILOUT = ' '
      CALL get_argument(iarg+1,FILOUT, LFIL)
      CMPSIZ = 0
      IOWS   = iowdsz()
      MODE = EXCLOB
      if (l64bit) then
        MODE = MODE + EX_ALL_INT64_DB + EX_ALL_INT64_API
      end if
      ndbout = excre(filout(:lfil), MODE, CMPSIZ, IOWS, IERR)
      if (ierr .lt. 0) then
         call exopts (EXVRBS, ierr)
         call exerr('grepos', 'Error from excre', ierr)
         go to 50
      endif
      if (l64bit) then
C ... Compress the output
        call exsetopt(ndbout, EX_OPT_COMPRESSION_LEVEL, 1, ierr)
        call exsetopt(ndbout, EX_OPT_COMPRESSION_SHUFFLE, 1, ierr)
      end if

C   --Write the QA records
      CALL DBOQA (NDBOUT, QAINFO, NQAREC, c(kqarec),
     &     NINFO, c(kinfo), ' Gen3D: ', FILIN)

C   --Write the initial variables

      CALL NEWINI (IDNSET(0,1)+IDNSET(0,2), IDESET(0,1)+IDESET(0,2),
     &   NSSUR, NESUR, C(KBKTYP), IA(KIBPAR))
      call expini (ndbout, title, ndim3, numnp3, numel3, nelbl3,
     &     nnps3, ness3, ierr)
      if (ierr .lt. 0) then
         call exerr('gen3d2', 'Error from expini', exlmsg)
         go to 40
      endif

      CALL DBPINI ('NTIS', NDBOUT, TITLE, NDIM3, NUMNP3, NUMEL3, NELBL3,
     &     NNPS3, LNPSN3, LNPSN3, NESS3, LESSE3, LESSN3, LESSN3,
     &     IDUM, IDUM, IDUM, FILOUT)

C   --Write the coordinates

      CALL MDRSRV ('ZCORD', KZCORD, NNREPL)
      CALL MDRSRV ('SINANG', KSINA, NNREPL)
      CALL MDRSRV ('COSANG', KCOSA, NNREPL)

      CALL MDRSRV ('XN3', KXN3, NUMNP3)
      CALL MDRSRV ('YN3', KYN3, NUMNP3)
      CALL MDRSRV ('ZN3', KZN3, NUMNP3)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL NEWXYZ (A(KXN), A(KYN), A(KXN3), A(KYN3), A(KZN3),
     &     IA(KIXNP), IA(KNRNP), IA(KNPCEN), A(KZCORD),
     &     A(KSINA), A(KCOSA), A)
      call expcor (ndbout, a(kxn3), a(kyn3), a(kzn3), ierr)
      if (ierr .lt. 0) then
         call exerr('gen3d2', 'Error from expcor', exlmsg)
         go to 40
      endif

      NAMECO(1) = 'X'
      NAMECO(2) = 'Y'
      NAMECO(3) = 'Z'
      call expcon(ndbout, nameco, ierr)
      if (ierr .lt. 0) then
         call exerr('gen3d2', 'Error from expcon', exlmsg)
         go to 40
      endif

      CALL MDDEL ('ZCORD')
      CALL MDDEL ('SINANG')
      CALL MDDEL ('COSANG')

      CALL MDDEL ('XN')
      CALL MDDEL ('YN')
      CALL MDDEL ('XN3')
      CALL MDDEL ('YN3')
      CALL MDDEL ('ZN3')
      CALL MDRSRV ('NUMELB3', KNELB3, NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL WRELB (A, IA, C(KBKTYP), C(KNMLB), IA(KIBPAR),
     &     IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &     IA(KLINK), A(KATRIB), A(KELATT),
     &     IA(KIXEL), IA(KINCEL), IA(KNREL), IA(KIECOL), IA(KIXNP),
     &     IA(KNRNP), IA(KNELB3))

      CALL MDDEL ('LINK')
      CALL MDDEL ('ATRIB')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Write the node sets

      CALL WRNPS (A, IA, IDNSET(0,1), IDNSET(0,2),
     &     IA(KIDNS), IA(KNNN3), IA(KIXNN3), IA(KLTNN3), A(KFACN3),
     &     IA(KIXNP), IA(KNRNP), *40)

      CALL MDDEL ('IDNPS')
      CALL MDDEL ('NNNP3')
      CALL MDDEL ('IXNNP3')
      CALL MDDEL ('LTNNP3')
      CALL MDDEL ('FACNP3')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Fixup sides sets if mirrored
      IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
        CALL MDRSRV ('IDXELB', KIDXELB, NELBLK+1)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
        CALL MIRSS (IDESET(0,1), IDESET(0,2),
     &    NESUR, IA(KISFRO), IA(KISBCK), IA(KLTES3), IA(KLTSS3),
     *    COMTOP, C(KNMLB), IA(KNELB3), IA(KIDXELB))
        CALL MDDEL ('IDXELB')
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      END IF
      CALL MDDEL('NUMELB3')
C     --Write the side sets

      CALL WRESS (A, IA, IDESET(0,1), IDESET(0,2),
     &     IA(KISFRO), IA(KISBCK), NSSUR, NESUR, IA(KNSFRO), IA(KNSBCK),
     &     IA(KIDSS), IA(KNES3), IA(KNNS3),
     &     IA(KIXES3), IA(KIXNS3), IA(KLTES3), IA(KLTSS3),
     &     IA(KLTNS3), A(KFACS3), *40)

      IF (NSSET .GT. 0) THEN
         CALL MDDEL ('ISSFRO')
         CALL MDDEL ('ISSBCK')
         CALL MDDEL ('NSSFRO')
         CALL MDDEL ('NSSBCK')
      END IF
      CALL MDDEL ('IDESS')
      CALL MDDEL ('NEES3')
      CALL MDDEL ('NNES3')
      CALL MDDEL ('IXEES3')
      CALL MDDEL ('IXNES3')
      CALL MDDEL ('LTEES3')
      CALL MDDEL ('LTNES3')
      CALL MDDEL ('FACES3')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      GOTO 50

   40 CONTINUE
      CALL MEMERR
      GOTO 50

   50 CONTINUE
      call exclos(ndbin, ierr)
      call exclos(ndbout, ierr)

   60 CONTINUE
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

   70 FORMAT (/
     &   16X,'  GGGGGG    EEEEEEEEEE  NN      NN   3333333   ',
     $     ' DDDDDDD  ', /
     &   15X,' GGGGGGGG   EEEEEEEEEE  NN      NN  333333333  ',
     $     ' DDDDDDDD ', /
     &   14X,'GG      GG  EE          NNN     NN  33      33 ',
     $     ' DD     DD', /
     &   13X,'GG          EE          NNNN    NN          33 ',
     $     ' DD     DD')
   80 FORMAT (
     &   12X,'GG          EEEEEEEE    NN NN   NN      33333  ',
     $     ' DD     DD', /
     &   11X,'GG    GGGG  EEEEEEEE    NN  NN  NN      33333  ',
     $     ' DD     DD', /
     &   10X,'GG    GGGG  EE          NN   NN NN          33 ',
     $     ' DD     DD', /
     &    9X,'GG      GG  EE          NN    NNNN  33      33 ',
     $     ' DD     DD', /
     &    8X,' GGGGGGGG   EEEEEEEEEE  NN     NNN  333333333  ',
     $     ' DDDDDDDD ', /
     &    7X,'  GGGGGG    EEEEEEEEEE  NN      NN   3333333   ',
     $     ' DDDDDDD  II')
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

C...Check whether model contains elements of a single topology.
C   This is currently used in the sideset mirroring code
      subroutine chktop(nelblk, blktyp, comtop)

      include 'exodusII.inc'
      integer nelblk
      character*(MXSTLN) blktyp(nelblk)
      character*(MXSTLN) comtop

      comtop = blktyp(1)(:3)
      do 10 i=2, nelblk
         if (blktyp(i)(:3) .ne. comtop(:3)) then
            comtop = 'MULTIPLE_TOPOLOGIES'
            return
         end if
 10   continue
      return
      end

