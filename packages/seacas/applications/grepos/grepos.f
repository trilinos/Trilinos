C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      PROGRAM GREPOS
C=======================================================================
C     
C     --*** GREPOS *** (GREPOS) GENESIS Positioning Program
C     --   Written by Greg Sjaardema - revised 03/07/89
C     --   Modified from GEN3D
C     --
C     --Expected input:
C     --   o The commands on the standard input device.
C     --   o The 3D GENESIS database on unit 9
C     --
C     --Output:
C     --   o A listing of the input database information and any errors
C     --     found on the standard output device.
C     --   o The repositioned 3D GENESIS database on unit 10.
C     --
C     --Developed at Sandia National Laboratories.
C     --
C     --Current author and code sponsor: Greg Sjaardema
C     --
C     --Revision History:
C     --   Modified from GEN3D [04/86 Created (Amy Gilkey)]
C     --
C     --Source is in FORTRAN 77
C     --
C     --External software used:
C     --  SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C     --  ExodusII library
C     --  NetCDF library (with modified parameters)      
C     --
C     --Documentation:
C     --   none

      include 'exodusII.inc'

      include 'gp_namlen.blk'
      include 'gp_progqa.blk'
      include 'gp_dbase.blk'
      include 'gp_dbtitl.blk'
      include 'gp_dbnums.blk'
      include 'gp_xyzoff.blk'
      include 'gp_xyzrot.blk'
      include 'gp_xyzmir.blk'
      include 'gp_xyzwrp.blk'
      include 'gp_nsset.blk'
      include 'gp_smooth.blk'
      include 'gp_snap.blk'
      include 'gp_combine.blk'
      include 'gp_deform.blk'
      include 'gp_attrot.blk'
      INCLUDE 'argparse.inc'
      
      CHARACTER*2048 FILIN, FILOUT, SCRATCH, SYNTAX
      CHARACTER*80 SCRSTR

C... String containing name of common element topology in model
C    or 'MULTIPLE_TOPOLOGIES' if not common topology.
      character*(MXSTLN) comtop

      LOGICAL EXODUS, NONQUD, ALLONE
      LOGICAL SMOOTH, SWPSS, USRSUB, CENTRD

      LOGICAL ISATRB, EXECUT

      LOGICAL RENEL, DELEL, DELNP

      INTEGER CMPSIZ, IOWS

      DIMENSION A(1)
      INTEGER IA(1)
      EQUIVALENCE (A(1), IA(1))
      CHARACTER*1 C(1)
      
C     --A - the dynamic numeric memory base array

      INCLUDE 'gp_qainfo.blk'
      CALL STRTUP (QAINFO)

      WRITE (*, 70)
      CALL BANNER (0, QAINFO,
     &     'A GENESIS DATABASE POSITIONING PROGRAM',
     &     ' ', ' ')
      call cpyrgt (0, '1990')

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C .. Get filename from command line.  If not specified, emit error message
      SYNTAX = 'Syntax is: "grepos [-name_length len] file_in file_out"'
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('CMDSPEC', SYNTAX(:LENSTR(SYNTAX)))
        GOTO 60
      else if (narg .gt. 4) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('CMDSPEC', SYNTAX(:LENSTR(SYNTAX)))
        GOTO 60
      end if

C ... Parse name_length option
      name_len = 0
      iarg = 1
      if (narg .eq. 4) then
        CALL get_argument(iarg,FILIN, LNAM)
        if (filin(:lnam) .eq. '-name_length') then
          CALL get_argument(iarg+1,FILIN, LNAM)
          read (filin(:lnam), '(i10)') name_len
          iarg = iarg + 2
        else
          SCRATCH = 'Unrecognized command option "'//FILIN(:LNAM)//'"'
          CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
          CALL PRTERR ('CMDSPEC', SYNTAX(:LENSTR(SYNTAX)))
          GOTO 60
        end if
      end if
C     --Open the input database and read the initial variables

      NDBIN = 9
      NDBOUT = 10
      
      CMPSIZ = 0
      IOWS   = 0

      FILIN  = ' '
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
         call exinq(ndbin, EXSSEL, lessel, rdum, cdum, ierr)
         call exinq(ndbin, EXSSDF, lessdf, rdum, cdum, ierr)
      else
         lessel = 0
         lessdf = 0
      end if

      call exinq(ndbin, EXDBMXUSNM, namlen, rdum, cdum, ierr)
      if (name_len .eq. 0) then
        if (namlen .gt. mxname) then
          namlen = mxname 
          maxnam = mxname
        else
          maxnam = namlen
        end if
      else
C ... Use user-specified length
        namlen = name_len 
        maxnam = name_len
      end if
      
      call exmxnm(ndbin, namlen, ierr)

C     --Reserve memory for the input information
      CALL MDRSRV ('XN', KXN, NUMNP)
      IF (NDIM .GE. 2) THEN
        CALL MDRSRV ('YN', KYN, NUMNP)
      END IF
      IF (NDIM .EQ. 3) THEN
         CALL MDRSRV ('ZN', KZN, NUMNP)
      ELSE
         CALL MDRSRV ('ZN', KZN, 1)
      END IF

      CALL MDRSRV ('MAPNN',  KMAPNN, NUMNP)
      CALL MDRSRV ('MAPEL',  KMAPEL, NUMEL)
      CALL MDRSRV ('ATRIB',  KATRIB, 0)
      CALL MDRSRV ('LINK',   KLINK,  0)

      CALL MDRSRV ('IDELB',  KIDELB, NELBLK)
      CALL MDRSRV ('NUMELB', KNELB,  NELBLK)
      CALL MDRSRV ('NUMLNK', KNLNK,  NELBLK)
      CALL MDRSRV ('NUMATR', KNATR,  NELBLK)
      CALL MCRSRV ('BLKTYP', KBKTYP, MXSTLN*NELBLK)
      CALL MCRSRV ('NAMEEB', KNAMEB, maxnam*NELBLK)

      CALL MDRSRV ('IDNPS',  KIDNS,  NUMNPS)
      CALL MDRSRV ('NNNPS',  KNNNS,  NUMNPS)
      CALL MDRSRV ('NDNPS',  KNDNPS, NUMNPS)
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
      CALL MDRSRV ('IXDNPS', KIXDNS, NUMNPS)
      CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
      CALL MDRSRV ('FACNPS', KFACNS, LNPSDF)
      CALL MCRSRV ('NAMENP', KNAMNP, maxnam*NUMNPS)

      CALL MDRSRV ('IDESS',  KIDSS,  NUMESS)
      CALL MDRSRV ('NEESS',  KNESS,  NUMESS)
      CALL MDRSRV ('NDESS',  KNDSS,  NUMESS)
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
      CALL MDRSRV ('IXDESS', KIXDSS, NUMESS)
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
      CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
      CALL MDRSRV ('LTSSNC', KLTSNC, LESSEL)
      CALL MDRSRV ('FACESS', KFACSS, LESSDF)
      CALL MCRSRV ('NAMESS', KNAMSS, maxnam*NUMESS)
      CALL MCRSRV ('NAMECO', KNAMCO, maxnam*ndim)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --Read information from the database and close file

      call exgcor (ndbin, a(kxn), a(kyn), a(kzn), ierr)

      call exopts (0, ierr)
      call exgnnm (ndbin, ia(kmapnn), ierr)
      call exgenm (ndbin, ia(kmapel), ierr)
      call exopts (EXVRBS, ierr)
      call getcon (ndbin, ndim, C(knamco), ierr)

      CALL DBIELB (NDBIN, '*', 1, NELBLK,
     &     IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &     A(1), IA(1), KLINK, KATRIB, C(KBKTYP), *60)
      
      call CHKTOP(NELBLK, C(KBKTYP), COMTOP)
      call getnam(NDBIN, 1, nelblk, C(KNAMEB))

C ... Attribute names...
C ... get total attribute count...
      numatt = 0
      do i=1, nelblk
        numatt = numatt + ia(knatr+i-1)
      end do
      
      call mcrsrv('NAMATT', KNAMATT, maxnam*NUMATT)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40
      call getattnam(ndbin, nelblk, ia(kidelb), ia(knatr), c(knamatt))
      
      if (numnps .gt. 0) then
         call exgcns(ndbin, ia(kidns), ia(knnns), ia(kndnps),
     &    ia(kixnns), ia(kixdns), ia(kltnns), a(kfacns), ierr)
         call getnam(NDBIN, 2, numnps, C(KNAMNP))
      end if
      if (numess .gt. 0) then
        call getss(ndbin, numess, ia(kidss), ia(kness), ia(kndss),
     &    ia(kixess), ia(kixdss), ia(kltess), ia(kltsss),
     &    ia(kltsnc), kfacss, a, allone, lessdf, *50)
        call getnam(NDBIN, 3, numess, C(KNAMSS))
      end if
      
      EXODUS = .TRUE.
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
      
C     --Read the database time steps
c     determine how many time steps are stored
      call exinq (ndbin, EXTIMS, NSTEPS, rdum, cdum, ierr)
      call mdrsrv('TIMES', KTIMES, NSTEPS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

c     read time values at all time steps
      if (nsteps .gt. 0) then
         call exgatm (ndbin, a(ktimes), ierr)      
      end if

      exodus = (nsteps .gt. 0)

      IF (EXODUS) THEN
        CALL DBINAM (NDBIN, '*', NDIM, NELBLK, NUMNPS, NUMESS,
     &    NNDIM, NNELB, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &    IXGV, IXNV, IXEV, IXNSV, IXSSV,
     &    A, IA, KIEVOK, KNSVOK, KSSVOK, C, KNAMES, EXODUS, *40)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      ELSE
        nvargl = 0
        nvarnp = 0
        nvarel = 0
        nvarns = 0
        nvarss = 0
        ixgv   = 0
        ixnv   = 0
        ixev   = 0
        ixnsv  = 0
        ixssv  = 0
        kievok = 0
        knsvok = 0
        kssvok = 0
        knames = 0
        
      END IF
      
      CALL DBPINI ('NTISV', NDBIN, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &     NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, 
     &     LESSDF, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS, FILIN)

C     --Read in runtime parameters

C     --Reserve memory for offsets

      MBLK = MAX(NELBLK, NUMNPS)
      CALL MDRSRV ('XEXPL', KXEXPL, MBLK)
      IF (NDIM .GE. 2) THEN
        CALL MDRSRV ('YEXPL', KYEXPL, MBLK)
      ELSE
         CALL MDRSRV ('YEXPL', KYEXPL, 1)
      END IF
      IF (NDIM .EQ. 3) THEN
         CALL MDRSRV ('ZEXPL', KZEXPL, MBLK)
      ELSE
         CALL MDRSRV ('ZEXPL', KZEXPL, 1)
      END IF

C     -- Reserve memory for node equiv (equiv node X with Y handled in command)
      
      call mdrsrv ('IXNP', KIXNP, NUMNP)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     .. ISATRB is TRUE if the model contains attributes      
      CALL MDFIND('ATRIB', KATRIB, KALEN)
      ISATRB = (KALEN .NE. 0)

C     .. Allocate memory for attribute scaling
      CALL CNTATR(NELBLK, IA(KNATR), IATCNT)
      CALL MDRSRV('ATRSCL', KATRSC, 2*IATCNT)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40
      
C     .. Set up status arrays for user manipulation of element blocks and sets

      CALL MDRSRV ('IELBST', KIELBS, NELBLK)
      CALL MDRSRV ('INPSST', KINPSS, NUMNPS)
      CALL MDRSRV ('IESSST', KIESSS, NUMESS)
      CALL MDRSRV ('ITIMST', KITIMS, NSTEPS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     .. Set up status arrays for element to nodal variable conversion
      call mdrsrv ('INOD2EL', KINOD2EL, NVARNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

 25   CONTINUE

      CALL COMAND (NDBIN, EXECUT, 
     &     IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &     IA(KIDNS),  IA(KNNNS), IA(KNDNPS),IA(KIXNNS),IA(KIXDNS),
     $     IA(KLTNNS), A(KFACNS),
     &     IA(KIDSS),  IA(KNESS), IA(KNDSS), IA(KIXESS), IA(KIXDSS),
     &     IA(KLTESS), IA(KLTSSS), A(KFACSS),
     &     A(KXN), A(KYN), A(KZN),
     &     A(KXEXPL), A(KYEXPL), A(KZEXPL), MODBLK,
     &     ISATRB, A(KATRSC), IA(KIXNP), IA(KMAPNN),
     &     IA(KIELBS), IA(KINPSS), IA(KIESSS),
     &     NQAREC, C(KQAREC), NINFO, c(kinfo), c(kbktyp),
     *     c(knameb), c(knamnp), c(knamss), c(knamatt),
     &     c(knames+nvargl*maxnam), nvarnp, IA(KINOD2EL),
     &     SWPSS, SMOOTH, USRSUB, CENTRD,
     &     NSTEPS, A(KTIMES), IA(KITIMS), A, IA, *60)

C     ... Snap
      if (numsnp .gt. 0) then
         call mdrsrv('ssblk', issblk, nelblk)
         call mdrsrv('iscrn', iscrnp, numnp)
         call mdrsrv('iscre', iscrep, numel)
         do 30 i=1, numsnp
           if (ismtyp(i) .eq. ISNAP) then
              call snap(ndbin, a, ia, i, a(kxn), a(kyn), a(kzn), ndim,
     *         numnp, numess, ia(kidss), ia(kness), ia(kixess),
     *         nelblk, ia(kidelb), ia(knelb), ia(knlnk))
           else
              call move(ndbin, a, ia, i, a(kxn), a(kyn), a(kzn), ndim,
     *         numnp, numel, numess, ia(kidss), ia(kness), ia(kixess), 
     *         ia(kltess), 
     *         nelblk, ia(kidelb), ia(knelb), ia(knlnk), ia(klink),
     *         ia(issblk), ia(iscrnp), ia(iscrep))
           end if
 30      continue
         call mddel('ssblk')
         call mddel('iscrn')
         call mddel('iscre')
      end if

C     ... Warp
      if (iwarp .ne. 0) then
         if (ndim .ne. 3) then
            CALL PRTERR('PROGRAM',
     *           'Warp cannot be specified for 2D databases')
         else
            call wrpxyz(a(kxn), a(kyn), a(kzn), numnp,
     *           iwarp, nrmwrp, wrpdis)
         end if
      end if
      
C     --Get the new positions for the elements and nodes
      IF (XMIRR * YMIRR * ZMIRR .LT. 0.0) THEN
         CALL DBMIRR (1, NELBLK, IA(KIDELB), IA(KNELB), IA(KNLNK),
     *        IA(KLINK), c(kbktyp), NDIM, NONQUD)

         CALL MDRSRV ('iblock', kiblock, nelblk)
         call chkss(nelblk, ia(knelb), ia(kiblock))

         CALL MIRSS (IA(KIDSS), IA(KNESS), 
     *        IA(KIXESS), IA(KLTESS),
     *        IA(KLTSSS), 
     *        IA(KIBLOCK), c(kbktyp), ALLONE, COMTOP)

         call mddel('iblock')
      END IF

C ... Swaps the orientation of a sideset. Should only be used on shells...
      IF (SWPSS) THEN
         CALL SWPESS (NUMESS, IA(KIDSS), IA(KNESS), IA(KNDSS),
     *        IA(KIXESS), IA(KIXDSS), IA(KLTESS),
     *        IA(KLTSSS), IA(KLTSNC), A(KFACSS), ALLONE, NONQUD,
     *        COMTOP)
      END IF

      IF (USRSUB) THEN
         CALL MDRSRV ('XNEW', kxnew, numnp)
         CALL MDRSRV ('YNEW', kynew, numnp)
         if (ndim .eq. 3) then
            CALL MDRSRV ('ZNEW', kznew, numnp)
         else
            CALL MDRSRV ('ZNEW', kznew, 1)
         end if

         CALL XYZMOD (numnp, ndim, A(KXN), A(KYN), A(KZN),
     *        A(KXNEW), A(KYNEW), A(KZNEW))

         CALL CPYREA(NUMNP, A(KXNEW), A(KXN))
         CALL CPYREA(NUMNP, A(KYNEW), A(KYN))
         if (ndim .eq. 3) CALL CPYREA(NUMNP, A(KZNEW), A(KZN))
         CALL MDDEL ('XNEW')
         CALL MDDEL ('YNEW')
         CALL MDDEL ('ZNEW')
      END IF

      CALL NEWXYZ (A(KXN), A(KYN), A(KZN), NUMNP, NDIM, A)

C     --Offset or scale each block if specified (EXPLODE, SCALE BLOCK)
C     -- MODBLK = 1 if block explode
C     -- MODBLK = 2 if block scale
C     -- MODBLK = 3 if block randomize

      IF (MODBLK .EQ. 1 .OR. MODBLK .EQ. 2 .OR. MODBLK .EQ. 3) THEN
         CALL MDRSRV ('ICONOD', KICOND, NELBLK*NUMNP)
         CALL MDRSRV ('MATMAP', KMTMAP, NELBLK*NELBLK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

         CALL EXPXYZ (A(KXN), A(KYN), A(KZN), IA(KICOND),
     &        A(KXEXPL), A(KYEXPL), A(KZEXPL), IA(KMTMAP),
     &        NELBLK, IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KLINK),
     &        NUMNP, NDIM, MODBLK)
      END IF

C     -- MODBLK = 4 if nodeset randomize
      IF (MODBLK .EQ. 4) THEN
         CALL EXPXYZN (A(KXN), A(KYN), A(KZN), 
     &        A(KXEXPL), A(KYEXPL), A(KZEXPL), 
     &        NUMNPS, IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLTNNS),
     &        NUMNP, NDIM, MODBLK)
      END IF

C     ... Deform
      if (idefst .gt. 0) then
        if (.not. exodus) then
          CALL PRTERR ('FATAL',
     *      'Deform command requires displacements.')
         GOTO 60
        end if
        if (nsteps .lt. idefst) then
          CALL PRTERR ('FATAL',
     *      'Deformation step is larger than database step count.')
         GOTO 60
        end if
        if (nvarnp .lt. ndim) then
          CALL PRTERR ('FATAL',
     *      'There are not enough nodal displacement variables.')
         GOTO 60
        end if
        call mdrsrv('DISPX', kdispx, numnp)
        call mdrsrv('DISPY', kdispy, numnp)
        call mdrsrv('DISPZ', kdispz, numnp)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
        
        call deform(a(kxn), a(kyn), a(kzn), numnp, ndim,
     *    a(kdispx), a(kdispy), a(kdispz), ndbin, idefst)

        call mddel('DISPX')
        call mddel('DISPY')
        call mddel('DISPZ')
      end if

C     ... Scale the attributes
      IF (ISATRB) THEN
         CALL NEWATR (NELBLK, IA(KNATR), A(KATRSC), IA(KNELB),
     &        A(KATRIB))
      END IF

      IF (REVATT) THEN
        CALL ROTATR (NELBLK, NDIM, IA(KIDELB), C(KBKTYP), IA(KNATR),
     *    IA(KNELB), A(KATRIB))
      END IF

      IF (SMOOTH) THEN
         call mdrsrv ('COSIN', ICOSIN, numnp * ndim)
         call mdrsrv ('NODES', inodes, numnp)
         call mdrsrv ('ISBND', iisbnd, numnp)
         
         call dobnd (a(kxn), a(kyn), a(kzn), ia(knelb), ia(knlnk),
     &        ia(kidelb), ia(klink), a(icosin), ia(inodes), ia(iisbnd),
     &        numnp, ndim, nelblk)
         call mddel('COSIN')
         call mdrsrv('XSCR', ixscr, numnp)
         call mdrsrv('YSCR', iyscr, numnp)
         if (ndim .eq. 3) then 
            call mdrsrv('ZSCR', izscr, numnp)
         else
            izscr = 1
         end if
         call smogs (a(kxn), a(kyn), a(kzn), ia(knelb), ia(knlnk),
     &        ia(kidelb), ia(klink), ia(iisbnd), nelblk, numnp,
     &        nit, toler, r0,
     &        a(ixscr), a(iyscr), a(izscr), ia(inodes), ndim)
      END IF

C     ... Incremental execution mode (re-enter COMAND routine)
C     ... NOTE: the code from 'CALL COMAND' to here can be reexecuted.
C     if MDRSRV or MDDEL is called, make sure it can be run
C     multiple times.
      
      if (execut) go to 25
      
C ... Calculate centroid if requested...
      if (centrd) then
        if (exodus) then
          CALL PRTERR ('FATAL',
     *      'Centroid function does not work for exodus databases yet.')
         GOTO 60
        end if
        call mdrsrv ('CENTROID', KCENT, NDIM*NUMEL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

         call centroid(a(kxn), a(kyn), a(kzn),
     *     a(kcent), a(kcent+numel), a(kcent+2*numel), nelblk,
     *     ia(knelb), ia(knlnk), ia(klink), ndim)     
       end if

C     ... See if node equivalencing specified
      DELNP = .FALSE.

      if (equiv) then
        if (eqtoler .ge. 0.0) then
          call mdrsrv ('IX',   KIX,   NUMNP)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 40

          call matxyz(ndim, numnp, a(kxn), a(kyn), a(kzn), ia(kix),
     &      ia(kixnp), nmatch, eqtoler)
          CALL MDDEL('IX')
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 40
        else
          nmatch = 0
          do i=1, numnp
            if (ia(kixnp+i-1) .ne. i) nmatch = nmatch+1
          end do
        end if
        write (*,*) 'Equivalencing ', nmatch, ' nodes'
        DELNP = (NMATCH .NE. 0)
        
        IF (.NOT. DELNP) THEN
          CALL MDDEL ('IXNP')
        END IF
      end if
C     --"Munch" the element blocks

C     ... Save old counts that are needed for writing timesteps
      numel0  = numel
      nelblk0 = nelblk
      numnp0  = numnp
      
C     --location of original numelb, isevok arrays
      kidelb0 = kidelb
      knelb0  = knelb
      kievok0 = kievok
      
      I = INTCNT (0, IA(KIELBS), NELBLK)
      RENEL = (I .LT. NELBLK)
      NUMEL1 = NUMEL

      
      if (renel .or. delnp) then
        CALL MDRSRV ('MSCR', KMSCR, MAX(NUMEL0, NUMNP0))
      end if

      IF (RENEL) THEN
C     ... Reserve space for original NUMELB and ISEVOK arrays and copy
C     old array contents into new (Only needed if EXODUS)       
         IF (EXODUS) THEN
            CALL MDRSRV ('NUMELB0', KNELB0, NELBLK0)
            CALL MDRSRV ('IDELB0',  KIDELB0, NELBLK0)
            CALL MDFIND ('ISEVOK',  IDUM, LIEVOK)
            CALL MDRSRV ('ISEVOK0', KIEVOK0, LIEVOK)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 40
            CALL CPYINT (NELBLK0, IA(KNELB),  IA(KNELB0))
            CALL CPYINT (NELBLK0, IA(KIDELB), IA(KIDELB0))
            CALL CPYINT (LIEVOK,  IA(KIEVOK), IA(KIEVOK0))

C     ... Map from new var to old for mapping variables.
            CALL MDRSRV ('MAPL', KMAPL, NUMEL0)
            CALL MDRSRV ('MAPN', KMAPN, NUMNP0)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 40
            CALL INIMAP(NUMEL0, IA(KMAPL))
            CALL INIMAP(NUMNP0, IA(KMAPN))
         END IF
         
         CALL MDRSRV ('IXEL', KIXEL, NUMEL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

         CALL MDFIND ('LINK', IDUM, LLNK)
         CALL MDRSRV ('LINKO', KLINKO, LLNK)

         CALL MDFIND ('ATRIB', IDUM, LATR)
         CALL MDRSRV ('ATRIBO', KATRO, LATR)

         CALL MDRSRV ('IXELB', KIXELB, NELBLK)
         CALL MDRSRV ('JNELB', KJNELB, NELBLK)
         CALL MDRSRV ('ISCR',  KISCR,  NELBLK)
         CALL MCRSRV ('NAMSCR', KNMSC, MXSTLN*NELBLK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

C     ... NUMEL changed in this block (if elements deleted)
         CALL MUNELB (NELBLK, IA(KIELBS), NUMEL,
     &        IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &        IA(KLINK), A(KATRIB), IA(KLINKO), A(KATRO),
     &        IA(KIXEL), IA(KIXELB), IA(KJNELB), IA(KISCR),
     &        c(kbktyp), c(knmsc), LLINK, LATRIB, c(knamatt),
     *        c(knameb))

C     ... Fix up the truth table if the element block count changes... 
         if (exodus .and. nvarel .gt. 0 .and. nelblk .ne. nelblk0) then
            call muntt(nelblk0, nelblk, nvarel, 
     $           ia(kievok0), ia(kievok), ia(kielbs))
         end if

         CALL MDDEL ('LINKO')
         CALL MDDEL ('ATRIBO')
         CALL MDDEL ('IXELB')
         CALL MDDEL ('JNELB')
         CALL MDDEL ('ISCR')
         CALL MCDEL ('NAMSCR')
         CALL MDLONG('LINK',  KLINK,  LLINK) 
         CALL MDLONG('ATRIB', KATRIB, LATRIB) 
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
      END IF

      CALL MDDEL ('IELBST')

C     --Mark if any elements are deleted

      DELEL = NUMEL .LT. NUMEL1

      IF (DELEL) THEN

C     --Make up an index of nodes in the existing element blocks
         if (.not. delnp) then
            CALL MDRSRV ('IXNP', KIXNP, NUMNP)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 40
         end if

         N = NUMNP
         CALL ZMFIXD (NELBLK, IA(KNELB), IA(KNLNK), IA(KLINK),
     *        N, IA(KIXNP))

         DELNP = (N .LT. NUMNP)

         IF (.NOT. DELNP) THEN
            CALL MDDEL ('IXNP')
         END IF
      END IF

C     --Squeeze the coordinates

      IF (DELNP) THEN
C     ... NUMNP modified in this call
         CALL ZMXYZ (NDIM, NUMNP, IA(KIXNP), A(KXN), A(KYN), A(KZN))

         CALL MDLONG ('XN', KXN, NUMNP)
         CALL MDLONG ('YN', KYN, NUMNP)
         IF (NDIM .EQ. 3) THEN
            CALL MDLONG ('ZN', KZN, NUMNP)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
         CALL REMAP(NUMNP0, IA(KIXNP), IA(KMAPNN), IA(KMSCR))
      END IF

C     --Renumber the element map

      IF (RENEL) THEN
         CALL REMAP (NUMEL0, IA(KIXEL), IA(KMAPEL), IA(KMSCR))
         IF (EXODUS) THEN
            CALL REMAP (NUMEL0, IA(KIXEL), IA(KMAPL),  IA(KMSCR))
         END IF
      END IF

C     --Squeeze the element map

      IF (DELEL) THEN
C     ... Changes first argument
         NSAVE = NUMEL0
         CALL ZMMAP (NSAVE, IA(KMAPEL))
         CALL MDLONG ('MAPEL', KMAPEL, NUMEL)
         NSAVE = NUMEL0
         IF (EXODUS) THEN
            CALL ZMMAP (NSAVE, IA(KMAPL))
            CALL MDLONG ('MAPL',  KMAPL,  NUMEL)
         END IF
      END IF

C     --Renumber the element block nodes

      IF (DELNP) THEN
        NSAVE = NUMNP0
        CALL ZMMAP(NSAVE, IA(KMAPNN))
        CALL MDLONG ('MAPNN', KMAPNN, NUMNP)
        CALL RENELB (NELBLK, -999, IA(KIXNP),
     &    IA(KNELB), IA(KNLNK), IA(KLINK))
      END IF

C     --Renumber the nodal point set nodes

      IF (DELNP) THEN
         CALL RENIX (LNPSNL, -999, IA(KIXNP), IA(KLTNNS), .TRUE.)
         NSAVE = NUMNP0
         IF (EXODUS) THEN
            IF (.not. DELEL) THEN
              CALL MDRSRV ('MSCR', KMSCR, MAX(NUMEL0, NUMNP0))
            ENDIF
            CALL REMAP (NSAVE, IA(KIXNP), IA(KMAPN),  IA(KMSCR))
            NSAVE = NUMNP0
            CALL ZMMAP (NSAVE, IA(KMAPN))
         END IF
      END IF

C     --Renumber the element side set elements

      IF (RENEL) THEN
         CALL RENIX (LESSEL, -999, IA(KIXEL), IA(KLTESS), .TRUE.)
      END IF

      IF (RENEL) THEN
         CALL MDDEL ('IXEL')
      END IF

C     --"Munch" the nodal point sets

      I = INTCNT (0, IA(KINPSS), NUMNPS)

      NUMNPS0 = NUMNPS
      KIDNS0  = KIDNS
      KNNNS0  = KNNNS
      LNPSNL0 = LNPSNL
      KNSVOK0 = KNSVOK
      
      IF ((I .LT. NUMNPS) .OR. DELNP) THEN
         CALL MDRSRV ('LTNNPO', KLTNNO, LNPSNL)
         CALL MDRSRV ('FACNPO', KFACNO, LNPSNL)
         CALL MDRSRV ('IXNNPO', KIXNNO, NUMNPS)
         CALL MDRSRV ('NNNPO',  KNNNO,  NUMNPS)
         CALL MDRSRV ('ISCR',   KISCR,  NUMNPS)
         CALL MDRSRV ('IDNS0',  KIDNS0, NUMNPS0)
         CALL MDRSRV ('NNNPS0', KNNNS0, NUMNPS0)
         CALL MCRSRV ('NAMSCR', KNMSC, MXSTLN*NUMNPS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

         CALL CPYINT (NUMNPS0, IA(KIDNS), IA(KIDNS0))
         CALL CPYINT (NUMNPS0, IA(KNNNS), IA(KNNNS0))

         CALL MUNNPS (NUMNPS, IA(KINPSS), LNPSNL, LNPSDF,
     &        IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLTNNS), A(KFACNS),
     &        IA(KLTNNO), A(KFACNO), IA(KIXNNO), IA(KNNNO), IA(KISCR),
     *        C(KNMSC), C(KNAMNP))

         CALL MDDEL ('LTNNPO')
         CALL MDDEL ('FACNPO')
         CALL MDDEL ('IXNNPO')
         CALL MDDEL ('NNNPO')
         CALL MDDEL ('ISCR')
         CALL MCDEL ('NAMSCR')

C     --Squeeze the nodal point sets

         IF (DELNP) THEN
            CALL ZMNPS (NUMNPS, IA(KINPSS), LNPSNL, LNPSDF,
     *       IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLTNNS), A(KFACNS))
         END IF

C     ... Fix up the truth table if the nodeset count changes... 
         if (exodus .and. nvarns .gt. 0 .and. numnps .ne. numnps0) then
           CALL MDFIND ('ISNSVOK',  IDUM, LNSVOK)
           CALL MDRSRV ('ISNSVOK0', KNSVOK0, LNSVOK)
           CALL MDSTAT (NERR, MEM)
           IF (NERR .GT. 0) GOTO 40
           CALL CPYINT (LNSVOK,  IA(KNSVOK), IA(KNSVOK0))
         
           call muntt(numnps0, numnps, nvarns, 
     $       ia(knsvok0), ia(knsvok), ia(kinpss))

C ... check that the sidesets that are retained contain the same number
C     of faces that the original sidesets contain.  At the current time,
C     can only map sideset variables if the sidesets are the same...
            i1 = 0
            do i=0,numnps0-1
               if (ia(kinpss+i) .eq.0) then
                  if (IA(KNNNS0+i) .ne. IA(KNNNS+i1)) then
                     write (*,900) 'Nodeset', ia(kidns0+i)
                  end if
                  i1 = i1 + 1
               end if
            end do
         end if

         CALL MDLONG ('IDNPS', KIDNS, NUMNPS)
         CALL MDLONG ('NNNPS', KNNNS, NUMNPS)
         CALL MDLONG ('IXNNPS', KIXNNS, NUMNPS)
         CALL MDLONG ('LTNNPS', KLTNNS, LNPSNL)
         CALL MDLONG ('FACNPS', KFACNS, LNPSDF)
      END IF

      CALL MDDEL ('INPSST')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

C     --"Munch" the element side sets

      I = INTCNT (0, IA(KIESSS), NUMESS)

      NUMESS0 = NUMESS
      KIDSS0  = KIDSS
      KNESS0  = KNESS
      LESSEL0 = LESSEL
      KSSVOK0 = KSSVOK

      IF ((I .LT. NUMESS) .OR. DELEL) THEN
         CALL MDRSRV ('LTEESO', KLTESO, LESSEL)
         CALL MDRSRV ('LTSSO',  KLTSSO, LESSEL)
         CALL MDRSRV ('FACS0',  KFACS0, LESSDF)
         CALL MDRSRV ('IXEESO', KIXESO, NUMESS)
         CALL MDRSRV ('IXEDS0', KIXDS0, NUMESS)
         CALL MDRSRV ('NEESO',  KNESO,  NUMESS)
         CALL MDRSRV ('NEDS0',  KNDS0,  NUMESS)
         CALL MDRSRV ('ISCR',   KISCR,  NUMESS)
         CALL MDRSRV ('IDSS0',  KIDSS0, NUMESS0)
         CALL MDRSRV ('NEESS0', KNESS0, NUMESS0)
         CALL MCRSRV ('NAMSCR', KNMSC, MXSTLN*NUMESS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

         CALL CPYINT (NUMESS0, IA(KIDSS), IA(KIDSS0))
         CALL CPYINT (NUMESS0, IA(KNESS), IA(KNESS0))

         CALL MUNESS (NUMESS, IA(KIESSS), LESSEL, LESSDF,
     &     IA(KIDSS), IA(KNESS), IA(KNDSS), IA(KIXESS), IA(KIXDSS),
     &     IA(KLTESS), IA(KLTSSS), A(KFACSS),
     &     IA(KLTESO), IA(KLTSSO), A(KFACS0), IA(KIXESO), IA(KIXDS0),
     &     IA(KNESO), IA(KNDS0), IA(KISCR), C(KNMSC), C(KNAMSS))
         

         CALL MDDEL ('LTEESO')
         CALL MDDEL ('LTSSO')
         CALL MDDEL ('FACS0')
         CALL MDDEL ('IXEESO')
         CALL MDDEL ('IXEDS0')
         CALL MDDEL ('NEESO')
         CALL MDDEL ('NEDS0')
         CALL MDDEL ('ISCR')
         CALL MCDEL ('NAMSCR')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40

C     --Squeeze the element side sets

         IF (DELEL) THEN
            CALL ZMESS (NUMESS, ia(kiesss), LESSEL, LESSDF, 
     &       IA(KIDSS), IA(KNESS), IA(KNDSS), IA(KIXESS),
     *       IA(KIXDSS), IA(KLTESS), IA(KLTSSS), IA(KLTSNC), A(KFACSS))
         END IF

C     ... Fix up the truth table if the sideset count changes... 
         if (exodus .and. nvarss .gt. 0 .and. numess .ne. numess0) then
            CALL MDFIND ('ISSSVOK',  IDUM,    LSSVOK)
            CALL MDRSRV ('ISSSVOK0', KSSVOK0, LSSVOK)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 40
            call cpyint(lssvok, ia(kssvok),ia(kssvok0))

            call muntt(numess0, numess, nvarss, 
     $           ia(kssvok0), ia(kssvok), ia(kiesss))

C ... check that the sidesets that are retained contain the same number
C     of faces that the original sidesets contain.  At the current time,
C     can only map sideset variables if the sidesets are the same...
            i1 = 0
            do i=0,numess0-1
               if (ia(kiesss+i) .eq.0) then
                  if (IA(KNESS0+i) .ne. IA(KNESS+i1)) then
                     write (*,900) 'Sideset', ia(kidss0+i)
                  end if
                  i1 = i1 + 1
               end if
            end do
         end if

         CALL MDLONG ('IDESS', KIDSS, NUMESS)
         CALL MDLONG ('NEESS', KNESS, NUMESS)
         CALL MDLONG ('IXEESS', KIXESS, NUMESS)
         CALL MDLONG ('LTEESS', KLTESS, LESSEL)
         CALL MDLONG ('LTSESS', KLTSSS, LESSEL)
         CALL MDLONG ('FACESS', KFACSS, LESSDF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
      END IF

      CALL MDDEL ('IESSST')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40

      IF (DELNP) THEN
         CALL MDDEL ('IXNP')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
      END IF

      CALL MINMAX (NUMNP, A(KXN), XMIN, XMAX)
      CALL MINMAX (NUMNP, A(KYN), YMIN, YMAX)
      IF (NDIM .EQ. 3) THEN
         CALL MINMAX (NUMNP, A(KZN), ZMIN, ZMAX)
      END IF

      call limits('Output Mesh Limits:', ndim,
     &     xmin, xmax, ymin, ymax, zmin, zmax)

C     --Open the output database
      
      FILOUT = ' '
      CALL get_argument(iarg+1,FILOUT, LFIL)
      CMPSIZ = 0
      IOWS   = iowdsz()
      ndbout = excre(filout(:lfil), EXCLOB, CMPSIZ, IOWS, IERR)
      if (ierr .lt. 0) then
         call exopts (EXVRBS, ierr1)
         call exerr('grepos', 'Error from excre', ierr)
         go to 50
      endif
      call exmxnm(ndbout, maxnam, ierr)

C     --Write the QA records
      CALL DBOQA (NDBOUT, QAINFO, NQAREC, c(kqarec),
     &     NINFO, c(kinfo), ' Grepos:  ', FILIN)

C     --Write the initial variables
      call expini (ndbout, title, ndim, numnp, numel, nelblk, numnps,
     &     numess, ierr)

C     --Write the coordinates
      call expcor (ndbout, a(kxn), a(kyn), a(kzn), ierr)

C     --Write the node/element order maps
      call expenm (ndbout, ia(kmapel), ierr)
      call expnnm (ndbout, ia(kmapnn), ierr)

C     --Write the element block
      CALL DBOELB (NDBOUT, 1, NELBLK,
     &     IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &     IA(KLINK), c(kbktyp), A(KATRIB))
      call putnam(NDBOUT, 1, nelblk, C(KNAMEB))

      call putattnam(ndbout, nelblk, ia(kidelb), ia(knatr), c(knamatt))

C     --Write the node sets
      if (numnps .gt. 0) then
        do 90 i= 0, numnps-1
          if (lnpsdf .eq. 0) then
            ia(kndnps+i) = 0
            ia(kixdns+i) = 0
          else
            ia(kndnps+i) = ia(knnns+i)
            ia(kixdns+i) = ia(kixnns+i)
          end if
 90     continue
         call expcns (ndbout, ia(kidns), ia(knnns), ia(kndnps),
     &        ia(kixnns), ia(kixdns), ia(kltnns), a(kfacns), ierr)
         call putnam(NDBOUT, 2, numnps, C(KNAMNP))
      end if
      
C     --Write the side sets
      if (numess .gt. 0) then
         call putss(ndbout, numess, ia(kidss), ia(kness), ia(kndss),
     *        ia(kixess), ia(kixdss), ia(kltess), ia(kltsss),
     *        a(kfacss), *40)
         call putnam(NDBOUT, 3, numess, C(KNAMSS))
      end if

      inod2el = 0
      do i=1, nvarnp
        inod2el = inod2el + ia(kinod2el+i-1)
      end do
      if (inod2el .gt. 0) then
        call mdrsrv('ELEMTIZE', kelmtz, inod2el*numel)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      end if

      CALL DBPINI ('NTISV', NDBOUT, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &     NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSDF,
     &  NVARGL, NVARNP, NVAREL+inod2el, NVARNS, NVARSS, FILOUT)

C ... See if writing the CENTROID field to the database...
      if (centrd) then
        if (.not. exodus) then
          nvarel = ndim
          CALL MCRSRV('NAMES', KNAMES, maxnam*NVAREL)
          CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
          CALL MCSTAT (NERR, MEM)
          IF (NERR .GT. 0) GOTO 40
          call centnam(c(knames), ia(kievok), nelblk, nvarel, ndim)
          ixev = 1
        end if
      end if
      
C     --Write the database names
      if (inod2el .gt. 0) then
        call dbonam(ndbout, ndim, c(knamco),
     &    nvargl, ixgv, nvarnp, ixnv, 0, ixev,
     *    nvarns, ixnsv, nvarss, ixssv, c(knames))
        call mcrsrv('TNAME', KTNAM, maxnam*(nvarel+inod2el))
        call eltznam(c(ktnam), c(knames), ixnv, nvarnp,
     *    ia(kinod2el), ixev, nvarel)
        call dbonam(ndbout, 0, c(knamco), 
     &    0, ixgv, 0, ixnv, nvarel+inod2el, 1,
     *    0, ixnsv, 0, ixssv, c(ktnam))
        call mcdel('TNAME')
      else
        call dbonam(ndbout, ndim, c(knamco),
     &    nvargl, ixgv, nvarnp, ixnv, nvarel, ixev,
     *    nvarns, ixnsv, nvarss, ixssv, c(knames))
      end if

C     ... Truth Table.
      if (nvarel+inod2el .gt. 0) then
        CALL MDRSRV ('ITMP', ktmp, NELBLK * (NVAREL+INOD2EL))
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
        call dbott(ndbout, 'E', nelblk, nvarel, inod2el,
     *    ia(kievok), ia(ktmp))
        call mddel('itmp')
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      end if
      
      if (nvarns .gt. 0) then
        CALL MDRSRV ('ITMP', ktmp, NUMNPS * NVARNS)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
        call dbott(ndbout, 'M', numnps, nvarns, 0, ia(knsvok), ia(ktmp))
        call mddel('itmp')
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      end if
      
      if (nvarss .gt. 0) then
        CALL MDRSRV ('ITMP', ktmp, NUMESS * NVARSS)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
        call dbott(ndbout, 'S', numess, nvarss, 0, ia(kssvok), ia(ktmp))
        call mddel('itmp')
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 40
      end if
      
      if (centrd .and. .not. exodus) then
        istep = 1
        time = 0.0
        CALL DBOSTE (NDBOUT, ISTEP,
     &    NVARGL, NVARNP, NUMNP, NVAREL, 0, NELBLK,
     &    IA(KNELB), IA(KIEVOK), IA(KIDELB),
     *    NVARNS, NUMNPS, IA(KNNNS), IA(KNSVOK), IA(KIDNS),
     *    NVARSS, NUMESS, IA(KNESS), IA(KSSVOK), IA(KIDSS),
     *    TIME, A(KVARGL), A(KVARNP), A(KCENT), A(KVARNS), A(KVARSS),
     $    A(1))
        
        call mddel('CENTROID')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
      end if
      if (.not. exodus) GOTO 50

C     --Read the database time steps
      MXV = MAX(NVARNP*NUMNP0, NVAREL*NUMEL0,
     *  NVARNS*LNPSNL0, NVARSS*LESSEL0, INOD2EL*NUMNP0)
      CALL MDRSRV ('VARGL', KVARGL, NVARGL)
      CALL MDRSRV ('VARNP', KVARNP, NVARNP * NUMNP0)
      CALL MDRSRV ('VAREL', KVAREL, (NVAREL+INOD2EL) * NUMEL0)
      CALL MDRSRV ('VARNS', KVARNS, NVARNS * LNPSNL0)
      CALL MDRSRV ('VARSS', KVARSS, NVARSS * LESSEL0)
      IF (RENEL .OR. DELNP .OR. DELEL) THEN
         CALL MDRSRV ('VARSCR', KVARSC, MXV)
      ELSE
         CALL MDRSRV ('VARSCR', KVARSC, 0)
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 40
      
      WRITE (*, *)
      WRITE (*, *)
      
C ... NOTE: VARNP and VAREL are treated as doubly-dimensioned arrays
C           dimensioned as (NUMEL, NVAREL)
      do 110 istep = 1, nsteps
        if (ia(kitims+istep-1) .eq. 0) then
           CALL DBISTE (NDBIN, '*', istep,
     &      NVARGL,
     *      NVARNP, NUMNP0,
     *      NVAREL, NELBLK0, IA(KNELB0), IA(KIEVOK0), IA(KIDELB0),
     *      NVARNS, NUMNPS0, IA(KNNNS0), IA(KNSVOK0), IA(KIDNS0),
     *      NVARSS, NUMESS0, IA(KNESS0), IA(KSSVOK0), IA(KIDSS0),
     &      TIME,
     *      A(KVARGL), A(KVARNP), A(KVAREL), A(KVARNS), A(KVARSS), *120)
           
           if (inod2el .gt. 0) then
             ioff = 0
             inoff = 0
             do i=1, nvarnp
               if (ia(kinod2el+i-1) .gt. 0) then
                 call elementize(a(kvarnp+inoff), a(kelmtz+ioff),
     *             nelblk, ia(knelb), ia(knlnk), ia(klink))
                 ioff = ioff + numel
              end if
              inoff = inoff + numnp
             end do
           end if

           IF (RENEL .OR. DELNP .OR. DELEL) THEN
              if (nvarnp .gt. 0) then
                 CALL MAPVAR(NUMNP0, NUMNP, NVARNP, IA(KMAPN),
     &                A(KVARNP), A(KVARSC))
              end if
              if (nvarel .gt. 0) then
C     ... Pass in the NEW element block parameters (number elements/block),
C     number element blocks, and truth table.              
                 CALL MAPEV(NUMEL0, NUMEL, NVAREL, IA(KMAPL),
     &                A(KVAREL), A(KVARSC))
              end if

              if (inod2el .gt. 0) then
                 CALL MAPEV(NUMEL0, NUMEL, inod2el, IA(KMAPL),
     &                A(kelmtz), A(KVARSC))
              end if

           END IF

           CALL DBOSTE (NDBOUT, ISTEP,
     &          NVARGL, NVARNP, NUMNP, NVAREL, INOD2EL, NELBLK,
     $          IA(KNELB), IA(KIEVOK), IA(KIDELB),
     $          NVARNS, NUMNPS0, IA(KNNNS0), IA(KNSVOK0), IA(KIDNS0),
     $          NVARSS, NUMESS0, IA(KNESS0), IA(KSSVOK0), IA(KIDSS0),
     $          TIME, 
     $          A(KVARGL), A(KVARNP), A(KVAREL), A(KVARNS), A(KVARSS),
     $          A(KELMTZ))
           
           WRITE (*, 10000) ISTEP, TIME
10000      FORMAT (' ', I8, ' time steps processed.  Time = ',1PE10.3)
        end if
 110  continue
      
 120  CONTINUE
      WRITE (SCRSTR, '(I9)', IOSTAT=K) NSTEPS
      CALL SQZSTR (SCRSTR, LSTR)
      WRITE (*, 10010) SCRSTR(:LSTR)
10010 FORMAT (/, 4X, A,
     &     ' time steps have been written to the output database')
      
      GO TO 50
 40   CONTINUE
      CALL MEMERR
      GOTO 50

 50   CONTINUE
      call exclos(ndbin, ierr)
      call exclos(ndbout, ierr)

 60   CONTINUE
      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

 70   FORMAT (/
     &     14X,' GGGGG   RRRRRR   EEEEEEE  PPPPPP    OOOOO    SSSSSS'/
     &     14X,'GG   GG  RR   RR  EE       PP   PP  OO   OO  SS     '/
     &     14X,'GG       RR   RR  EE       PP   PP  OO   OO  SS     '/
     &     14X,'GG       RRRRRR   EEEEE    PPPPPP   OO   OO   SSSSS '/
     &     14X,'GG  GGG  RRRRR    EE       PP       OO   OO       SS'/
     &     14X,'GG   GG  RR  RR   EE       PP       OO   OO       SS'/
     &     14X,' GGGGG   RR   RR  EEEEEEE  PP        OOOOO   SSSSSS ')
 900  FORMAT(/,'WARNING: ',A,i5,' is a different size in the output',
     $     /,9x,'database than in the input database.  If there are',
     $     /,9x,'variables on this sideset, they will be transferred',
     $     /,9x,'incorrectly. Contact gdsjaar@sandia.gov',
     $     /,9x,'if you need this capability.')
      END

      SUBROUTINE INIMAP(LEN, MAP)
      INTEGER MAP(*)
      DO 10 I=1, LEN
        MAP(I) = I
 10   CONTINUE
      END
      subroutine exgqaw(ndb, qarec, ierr)
      include 'gp_params.blk'
      character*(mxstln) qarec(4, *)
      call exgqa(ndb, qarec, ierr)
      return
      end
      subroutine exginw(ndb, info, ierr)
      include 'gp_params.blk'
      character*(mxlnln) info(*)
      call exginf(ndb, info, ierr)
      return
      end

      subroutine dbott(ndbout, type, nblk, nvar, iextra, ievok, itmp)
      character*1 type
      logical IEVOK(nblk,nvar)
      INTEGER ITMP(NVAR+iextra,nblk)
      
      do 20 i=1, nvar
        do 10 ielb = 1, nblk
          if (ievok(ielb,i)) then
            itmp(i,ielb) = 1
          else
            itmp(i,ielb) = 0
          end if
 10     continue
 20   continue
      
      do i=1, iextra
        do ielb = 1, nblk
          itmp(i+nvar, ielb) = 1
        end do
      end do
      if (type .eq. 'E') then
        call expvtt(ndbout, nblk, nvar+iextra, itmp, ierr)
      else if (type .eq. 'M') then
        call expnstt(ndbout, nblk, nvar+iextra, itmp, ierr)
      else if (type .eq. 'S') then
        call expsstt(ndbout, nblk, nvar+iextra, itmp, ierr)
      end if
      return 
      end

      subroutine centnam(names, isevok, nelblk, nvarel, ndim)
      include 'gp_namlen.blk'
      character*(maxnam) names(*)
      integer isevok(nelblk, nvarel)
      
      names(1)(:maxnam) = 'centroid_x'
      if (ndim .ge. 2) names(2)(:maxnam) = 'centroid_y'
      if (ndim .eq. 3) names(3)(:maxnam) = 'centroid_z'
      
      do 20 i=1, nelblk
        do 10 j=1, nvarel
          isevok(i,j) = 1
 10     continue
 20   continue

      return
      end

      subroutine eltznam(evnames, names, ixnv, nvarnp,
     *  inod2el, ixev, nvarel)
      include 'gp_namlen.blk'
      character*(maxnam) names(*), evnames(*)
      integer inod2el(*)

      do i=1, nvarel
        evnames(i) = names(ixev+i-1)
      end do

      isum = 0
      do i=1, nvarnp
        if (inod2el(i) .gt. 0) then
          isum = isum + 1
          evnames(nvarel+isum)(:2) = 'n_'
          evnames(nvarel+isum)(3:maxnam) = names(ixnv+i-1)(:maxnam-3)
        end if
      end do

      return
      end

      SUBROUTINE PRTT(NELBLK, NVAREL, truth)
      LOGICAL TRUTH(NELBLK, NVAREL)
      DO 10 I=1, NELBLK
         WRITE (*,*) (TRUTH(I,IVAR),IVAR=1,NVAREL)
 10   CONTINUE
      RETURN
      END

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

      subroutine getcon(ndb, ndim, nameco, ierr)
      include 'gp_namlen.blk'

      character*(maxnam) nameco(*)
      CALL INISTR (NDIM, ' ', nameco)
      call exgcon (ndb, nameco, ierr)
      return 
      end

      subroutine getnam(ndb, itype, isiz, names)
      include 'gp_namlen.blk'
      character*(maxnam) names(*)
      
      call exgnams(ndb, itype, isiz, names, ierr)
      return 
      end

      subroutine putnam(ndb, itype, isiz, names)
      include 'gp_namlen.blk'
      character*(maxnam) names(*)
      
      call expnams(ndb, itype, isiz, names, ierr)
      return 
      end

