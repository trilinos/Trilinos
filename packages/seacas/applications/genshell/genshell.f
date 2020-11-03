C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      PROGRAM GENSHL
C=======================================================================
C     --*** GENSHL *** (GENSHL) GENESIS 2D to 3D Shell Program
C     --
C     --GENSHL inputs a 2D GENESIS database and outputs a 3D GENESIS database.
C     --
C     --Expected input:
C     --   o The commands on the standard input device.
C     --   o The 2D GENESIS database on unit 9
C     --     (must have 4 nodes per element,
C     --     cannot have more than 2048 element blocks).
C     --
C     --Output:
C     --   o A listing of the input database information and any errors
C     --     found on the standard output device.
C     --   o The 3D GENESIS Shell database on unit 10.

C     --Developed at Sandia National Laboratories.
C     --
C     --Current author and code sponsor: Gregory D. Sjaardema
C     --
C     --Revision History:
C     --   04/86 GEN3D Created (Amy Gilkey)
C     --   01/91 GENSHL Created (Greg Sjaardema)
C     --
C     --Source is in FORTRAN 77
C     --
C     --External software used:
C     --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C     --   SUPLIB package (EXODUS file manipulation and other routines)
C     --

C     --Documentation:
C     --   --NONE--

      include 'exodusII.inc'
      INCLUDE 'gs_progqa.blk'
      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbtitl.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'
      INCLUDE 'gs_xyzoff.blk'
      INCLUDE 'gs_xyzrot.blk'
      INCLUDE 'gs_xyzmir.blk'
      INCLUDE 'argparse.inc'

      CHARACTER*2048 FILIN, FILOUT, SCRATCH

      CHARACTER*(MXSTLN) NAMECO(6)
      CHARACTER*(MXSTLN) NAMELB(2048)
C     --NAMECO - the coordinate names
C     --NAMELB - the element block names

      CHARACTER CDUM

      DIMENSION A(1)
      INTEGER IA(1)
      EQUIVALENCE (A(1), IA(1))
      CHARACTER*1 C(1)
C     --A - the dynamic numeric memory base array

      INTEGER IDNSET(0:10,2)
      INTEGER IDESET(0:10,2)

      INCLUDE 'gs_qainfo.blk'

      CALL STRTUP (QAINFO)

      CALL BANNER (0, QAINFO,
     &     'A GENESIS DATABASE 2D TO 3D SHELL CONVERSION PROGRAM',
     &     ' ', ' ')
      call cpyrgt (0, '1991')

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "genshell 2dfilename 3dfilename"')
        GOTO 60
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "genshell 2dfilename 3dfilename"')
        GOTO 60
      end if

C     --Open the input database and read the initial variables

      NDBIN  = 9
      NDBOUT = 10

      CMPSIZ = 0
      IOWS   = 0

      FILIN = ' '
      CALL get_argument(1,FILIN, LNAM)
      NDBIN = exopen(filin(:lnam), EXREAD, CMPSIZ, IOWS, vers, IERR)
      IF (IERR .NE. 0) THEN
        SCRATCH = 'Database "'//FILIN(:LNAM)//'" does not exist.'
        CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
        GOTO 60
      END IF

      call exgini(ndbin, title, ndim, numnp, numel, nelblk,
     *     numnps, numess, ierr)
      if (nelblk .gt. 2048) then
        CALL PRTERR ('PROGRAM', 'Too many element blocks. 2048 MAX')
         GOTO 60
      end if
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
     &     LESSDF, IDUM, IDUM, IDUM)

      IF (NDIM .NE. 2) THEN
         CALL PRTERR ('FATAL', 'Number of dimensions must be 2')
         GOTO 60
      END IF

C     --Reserve memory for the 2D information

      CALL MDRSRV ('XN', KXN, NUMNP)
      CALL MDRSRV ('YN', KYN, NUMNP)
      KZN = 1

      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)

      CALL MDRSRV ('IDELB', KIDELB, NELBLK)
      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
      CALL MDRSRV ('NUMLNK', KNLNK, NELBLK)
      CALL MDRSRV ('NUMATR', KNATR, NELBLK)
      CALL MDRSRV ('LINK', KLINK, 0)
      CALL MDRSRV ('ATRIB', KATRIB, 0)

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

C     --Read 2D information from the database and close file

      call exgcor (ndbin, a(kxn), a(kyn), a(kzn), ierr)
      call exgmap (ndbin, ia(kmapel), ierr)

      CALL INISTR (NDIM, ' ', NAMECO)
      call exgcon (ndbin, nameco, ierr)

      CALL INISTR (2048, ' ', NAMELB)
      CALL DBIELB (NDBIN, '*', 1, NELBLK,
     *  IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     *  A, IA, KLINK, KATRIB, NAMELB, *60)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C ... At the current time, only handle all 4-node or 8-node hexes per model.
      NLNK = IA(KNLNK)
      DO 5 IEL = 2, NELBLK
        if (NLNK .NE. IA(KNLNK+IEL-1)) THEN
          CALL PRTERR ('FATAL',
     *      'Mesh must be ALL 4-node or ALL 8-node quads')
          STOP
        END IF
 5    CONTINUE

      if (numnps .gt. 0) then
         call exgcns(ndbin, a(kidns), a(knnns), a(kndnps), a(kixnns),
     *        a(kixdns), a(kltnns), a(kfacns), ierr)
      end if

      if (numess .gt. 0) then
         call exgcss(ndbin, a(kidss), a(kness), a(kndss), a(kixess),
     *        a(kixdss), A(KLTeSS), a(kltsss), a(kfacss), ierr)
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

C     -- This assumes only max of 7 attributes per element block
      CALL MDRSRV ('ATRIBNW', KATRIBN, 0)
      CALL MDRSRV ('ELATTR', KELATT, NELBLK*7)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --Read in runtime parameters

      CALL COMAND (IA(KIDNS), IA(KIDSS), IDNSET, IDESET,
     &     IA(KIDELB), IA(KNELB), IA(KNLNK), NAMELB, A(KELATT),
     &     A(KXN), A(KYN), A, *60)

C     --Get the new numbers for the elements and nodes

      LNPSNO = LNPSNL

C     --Get the node sets

      CONTINUE

C     --Get the side sets, and the front and back side sets

      NSSET = IDESET(0,1)+IDESET(0,2)
      IF (NSSET .GT. 0) THEN
         CALL MDRSRV ('ISSFRO', KISFRO, NUMEL)
         CALL MDRSRV ('ISSBCK', KISBCK, NUMEL)
         CALL MDRSRV ('NSSFRO', KNSFRO, NLNK*NUMEL)
         CALL MDRSRV ('NSSBCK', KNSBCK, NLNK*NUMEL)
      ELSE
         KISFRO = 1
         KISBCK = 1
         KNSFRO = 1
         KNSBCK = 1
      END IF

      LESSEO = INTADD (NUMESS, IA(KNESS))
      LESSNO = INTADD (NUMESS, IA(KNDSS))

      CALL NEWESS
     &     (IDESET(0,1), IDESET(0,2), NSSUR, NLNK,
     &     IA(KLINK), IA(KISFRO), IA(KISBCK), IA(KNSFRO), IA(KNSBCK))

      IF (NSSET .GT. 0) THEN
         CALL MDLONG ('NSSFRO', KNSFRO, NSSUR)
         CALL MDLONG ('NSSBCK', KNSBCK, NSSUR)
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --Open the output database
      FILOUT = ' '
      CALL get_argument(2,FILOUT, LFIL)
      CMPSIZ = 0
      IOWS   = iowdsz()
      ndbout = excre(filout(:lfil), EXCLOB, CMPSIZ, IOWS, IERR)
      if (ierr .lt. 0) then
         call exopts (EXVRBS, ierr)
         call exerr('grepos', 'Error from excre', ierr)
         go to 50
      endif

C     --Write the initial variables

      CALL NEWINI (IDNSET(0,1)+IDNSET(0,2), IDESET(0,1)+IDESET(0,2),
     &     NSSUR, IA(KNATR))
      call expini (ndbout, title, ndim3, numnp3, numel3, nelbl3,
     &     nnps3, ness3, ierr)

      CALL DBPINI ('NTIS', NDBOUT, TITLE, NDIM3, NUMNP3, NUMEL3, NELBL3,
     &  NNPS3, LNPSN3, LNPSN3, NESS3, LESSE3, LESSN3, LESSN3,
     &  IDUM, IDUM, IDUM)

C     --Write the coordinates

      CALL MDRSRV ('XN3', KXN3, NUMNP3)
      CALL MDRSRV ('YN3', KYN3, NUMNP3)
      CALL MDRSRV ('ZN3', KZN3, NUMNP3)

C     -- Since we only read quads, NUMATR should be 0.  For a shell, we
C     have 1 attribute (thickness)

      CALL MDLONG ('ATRIB', KATRIB, NUMEL3)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      CALL NEWXYZ (A(KXN), A(KYN), A(KXN3), A(KYN3), A(KZN3),
     $     A(KATRIB), A )
      call expcor (ndbout, a(kxn3), a(kyn3), a(kzn3), ierr)

      NAMECO(1) = 'X'
      NAMECO(2) = 'Y'
      NAMECO(3) = 'Z'
      call expcon(ndbout, nameco, ierr)

      CALL MDDEL ('XN')
      CALL MDDEL ('YN')
      CALL MDDEL ('XN3')
      CALL MDDEL ('YN3')
      CALL MDDEL ('ZN3')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --Write the element order map - NOTE: New map is same as old map
      call expmap (ndbout, ia(kmapel), ierr)
      CALL MDDEL ('MAPEL')

C     --Fixup connectivity if mirrored
      IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
         CALL DBMIRR (1, NELBLK, IA(KIDELB), IA(KNELB),
     $        IA(KNLNK), IA(KLINK))
      END IF
C     --Write the element block

      ILOFF = 0
      IAOFF = 0
      DO 30 I = 1, NELBLK
C ... Calculate these here since values may change if BEAM...
        II = I - 1
        ILINC = IA(KNELB+II) * IA(KNLNK+II)
        IAINC = IA(KNELB+II) * IA(KNATR+II)
        IF (NAMELB(I)(:4) .EQ. 'QUAD' .OR. NAMELB(I) .EQ. ' ') THEN
          NAMELB(I) = 'SHELL'
        ELSE IF (NAMELB(I)(:4) .EQ. 'BEAM') THEN
          IA(KNATR+I-1) = 7
        ELSE
          CALL PRTERR ('CMDWARN',
     $      'Invalid Element Block Type')
        END IF
        CALL DBOELB (A, NDBOUT, IA(KIDELB+II), IA(KNELB+II),
     $    IA(KNLNK+II), IA(KNATR+II), IA(KLINK+ILOFF),
     *    NAMELB(I), A(KATRIB), A(KELATT+(II*7)))
        ILOFF = ILOFF + ILINC
        IAOFF = IAOFF + IAINC

 30   CONTINUE

      CALL MDDEL ('LINK')
      CALL MDDEL ('ATRIB')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --Write the node sets

      CALL WRNPS (A, IA, IDNSET(0,1), IDNSET(0,2), IA(KIDNS),
     &  IA(KNNNS), IA(KIXNNS), IA(KLTNNS), A(KFACNS), *40)

C     --Fixup sides sets if mirrored
c$$$      IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
c$$$         CALL MIRSS (IDESET(0,1), IDESET(0,2), NLNK,
c$$$     &        NSSUR, IA(KNSFRO), IA(KNSBCK), IA(KLTSSS))
c$$$      END IF
C     --Write the side sets

      CALL WRESS (A, IA, IDESET(0,1), IDESET(0,2),
     &     IA(KISFRO), IA(KISBCK), NSSUR, IA(KNSFRO), IA(KNSBCK),
     &     IA(KIDSS), IA(KNESS), IA(KNDSS), IA(KIXESS), IA(KIXDSS),
     &     IA(KLTESS), IA(KLTSSS), A(KFACSS), *40)

      IF (NSSET .GT. 0) THEN
         CALL MDDEL ('ISSFRO')
         CALL MDDEL ('ISSBCK')
         CALL MDDEL ('NSSFRO')
         CALL MDDEL ('NSSBCK')
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C   --Write the QA records
      CALL DBOQA (NDBOUT, QAINFO, NQAREC, c(kqarec),
     &     NINFO, c(kinfo), ' GenShell: ', FILIN)

      GOTO 50

 40   CONTINUE
      CALL MEMERR
      GOTO 50

 50   CONTINUE
      call exclos(ndbin,  ierr)
      call exclos(ndbout, ierr)

 60   CONTINUE
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

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
