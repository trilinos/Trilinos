C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      PROGRAM EX2EX1V2
C=======================================================================

C   --*** EX2EX1V2 *** EXODUS II to EXODUS I translator
C   --
C   --EX2EX1V2 reads the EXODUS II V2.02 and V2.03
C   --regular and history files and writes an EXODUS I database file.
C   --
C   --Expects the output database on unit 11.

      include 'exodusII.inc'
      INCLUDE 'argparse.inc'

      CHARACTER*8 QAINFO(6)
      PARAMETER (MAXQA = 100, MAXINF = 100)
c      CHARACTER*32 QAREC(4,MAXQA)
c      CHARACTER*80 INFO(MAXINF)

C ... Names read in are 32-characters long
      CHARACTER*(mxstln) MAMECO(6)
      CHARACTER*(mxstln) MAMES(256)
C ... Names written out are 8-characters long, truncate with no warning
      CHARACTER*8 NAMECO(6)
      CHARACTER*8 NAMELB(256)
      CHARACTER*8 NAMES(256)

      CHARACTER*80 TITLE

      DIMENSION A(1), ia(1)
C      --A - the dynamic memory base array
      equivalence (a(1), ia(1))
      CHARACTER*1 c(1)
      CHARACTER*8 cdummy

      CHARACTER*5 STRA, STRB
      CHARACTER*8 STR8
      character*2048 netfil, ndbfil, errmsg
      character*(mxstln) name
      LOGICAL WHOTIM
      real wtime, htime
      integer cpuws, iows
      LOGICAL MDEBUG

      data (qainfo(i), i=1,3) / 'ex2ex1v2', '20110616', 'v 2.08  ' /
      data cpuws, iows /0,0/

      CALL STRTUP (QAINFO)

      CALL BANNER (0, QAINFO,
     &   'EXODUS II TO EXODUS I DATABASE'//
     &   ' TRANSLATOR',' ', ' ')
      call exinq (netid, EXLBVR, idummy, exlibversion, name, nerr)
      write(*,'(A,F6.3)')'ExodusII Library version ',
     1          exlibversion

      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      MDEBUG = .false.
      if (MDEBUG) then
         call mlist()
      end if

c       make netCDF and exodus errors not show up

c      call ncpopt (0)
      call exopts (0,ierr)

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "ex2ex1v2 exo2_file exo1_file"')
        GOTO 140
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "ex2ex1v2 exo2_file exo1_file"')
        GOTO 140
      end if

c       open the netcdf file

      net = 11
      CALL get_argument(1,netfil, lnam)

      netid = EXOPEN(netfil(1:lnam), EXREAD, cpuws, iows, vers, nerr)
      if (nerr .lt. 0) then
        errmsg = 'Database "'//netfil(:lnam)//'" does not exist.'
        CALL PRTERR ('FATAL', errmsg(:lenstr(errmsg)))
        call exerr('ex2ex1v2', errmsg, exlmsg)
        goto 140
      endif

      write(*,*) 'Input file name: ',netfil(1:lnam)
      call exinq (netid, EXVERS, idummy, exversion, name, nerr)
      write(*,'(A,F6.3)')
     & 'This database was created by ExodusII version ', exversion

C       open the output database and write the initial variables

      NDB = 20
      CALL get_argument(2,ndbfil, lnam)
      open(unit=ndb, file=ndbfil(:lnam), form='unformatted',
     *     status='unknown', iostat=ierr)
       IF (IERR .NE. 0) THEN
         errmsg = 'Error opening output file "'//ndbfil(:lnam)//'".'
         CALL PRTERR ('FATAL', errmsg(:LENSTR(errmsg)))
         GOTO 140
      END IF
      write(*,*) 'Output file name: ',ndbfil(1:lnam)

c       get initialization parameters from regular netcdf file

      CALL EXGINI (netid, title, ndim, numnp, numel,
     &       nelblk, numnps, numess, nerr)
      if (nerr .lt. 0) then
        call exerr('ex2ex1v2', 'Error from exgini', exlmsg)
        goto 140
      endif

c       get the length of the node sets node list

      if (numnps .gt. 0) then
         CALL EXINQ (netid, EXNSNL, lnpsnl, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exqini', exlmsg)
           goto 140
         endif
      else
         lnpsnl = 0
      endif

      if (numess .gt. 0) then

c       get the length of the side sets node list

        CALL EXINQ (netid, EXSSNL, lessnl, dummy, cdummy, nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exqini', exlmsg)
           goto 140
        endif

c       get the length of the side sets distribution factor list

         CALL EXINQ (netid, EXSSDF, lessdl, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exqini', exlmsg)
           goto 140
         endif

c       get the length of the side sets element list

         CALL EXINQ (netid, EXSSEL, lessel, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exqini', exlmsg)
           goto 140
         endif
      else
         lessnl = 0
         lessel = 0
         lessdl = 0
      endif

c       write the initialization information to the EXODUS 1.0 database

      CALL DBOINI (NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL)

      CALL DBPINI ('TIS', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL,
     &   IDUM, IDUM, IDUM, IDUM)

C   --Read the coordinates

      CALL MDRSRV ('XN', KXN, NUMNP)
      CALL MDRSRV ('YN', KYN, NUMNP)
      IF (NDIM .GE. 3) THEN
        CALL MDRSRV ('ZN', KZN, NUMNP)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 130
c        write(*,*)' ************************* NDIM: ',ndim
        CALL EXGCOR(netid, a(kxn), a(kyn), a(kzn), nerr)
         if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgcor', exlmsg)
           goto 140
         endif

        CALL DBOXYZ (NDB, NDIM, NUMNP, A(KXN), A(KYN), A(KZN))

        CALL MDDEL ('XN')
        CALL MDDEL ('YN')
        CALL MDDEL ('ZN')
      ELSE
        CALL EXGCOR(netid, a(kxn), a(kyn), dummy, nerr)
         if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgcor', exlmsg)
           goto 140
         endif

        CALL DBOXYZ (NDB, NDIM, NUMNP, A(KXN), A(KYN), dummy)

        CALL MDDEL ('XN')
        CALL MDDEL ('YN')
      ENDIF

C   --Read the element order map

      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL EXGMAP (netid, a(KMAPEL), nerr)
      if (nerr .ne. 0) then
         if (nerr .eq. 17) then

C   -- no element order map in the EXODUS II file; create a dummy one
            do 30 i=1,numel
               ia(kmapel+i-1) = i
30          continue
         else
            goto 140
         endif
      endif

      CALL DBOMAP (NDB, NUMEL, A(KMAPEL))

      CALL MDDEL ('MAPEL')

c       Read in the element block ID array

      call MDRSRV ('IDELB', kidelb, nelblk)
      call exgebi (netid, a(kidelb), nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgebi', exlmsg)
         goto 140
      endif

C   --Read the element blocks

      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
      CALL MDRSRV ('LINK', KLINK, 0)
      CALL MDRSRV ('ATRIB', KATRIB, 0)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      nel = 0
      DO 50 IELB = 1, NELBLK

         CALL EXGELB (netid, a(kidelb+ielb-1), name,
     &        a(knelb+ielb-1), numlnk, numatr, nerr)
         if (nerr .lt. 0) then
            call exerr('ex2ex1v2', 'Error from exgelb', exlmsg)
            goto 140
         endif
         namelb(ielb) = name(:8)

         call getin (ia(knelb+ielb-1),num)
         if (numlnk .gt. 0) then
           CALL MDLONG ('LINK', KLINK, num*numlnk)
           CALL EXGELC (netid, a(kidelb+ielb-1),
     &         a(klink), nerr)
           if (nerr .lt. 0) then
              call exerr('ex2ex1v2', 'Error from exgelc', exlmsg)
              goto 140
           endif
         end if

         if (numatr .gt. 0) then
           CALL MDLONG ('ATRIB', KATRIB, num*numatr)
           CALL EXGEAT (netid, a(kidelb+ielb-1), a(katrib), nerr)
           if (nerr .lt. 0) then
              call exerr('ex2ex1v2', 'Error from exgeat', exlmsg)
              goto 140
           endif
         end if

         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 130

         CALL DBOELB (NDB, IELB, IELB,
     &      a(kidelb+ielb-1), A(KNELB+IELB-1), NUMLNK, NUMATR,
     &      A(KLINK), A(KATRIB))

         nel=nel+num
         CALL MDLONG ('LINK', klink, 0)
         CALL MDLONG ('ATRIB', katrib, 0)

50    CONTINUE

      CALL MDDEL ('LINK')
      CALL MDDEL ('ATRIB')

      IF (NEL .NE. NUMEL) THEN
         CALL INTSTR (1, 0, NEL, STRA, LSTRA)
         CALL INTSTR (1, 0, NUMEL, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'NUMBER OF ELEMENTS IN BLOCK = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

C   --Read the node sets

      CALL MDRSRV ('IDNPS',  KIDNS, NUMNPS)     ! Node set ids array
      CALL MDRSRV ('NNNPS', KNNNS, NUMNPS)      ! Node set node count array
      CALL MDRSRV ('NDNPS', KNDNS, NUMNPS)      ! Node set df count array
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)    ! Node set nodes index array
      CALL MDRSRV ('IXDNPS', KIXDNS, NUMNPS)    ! Node set df index array
      CALL MDRSRV ('LSTNPS', KLSTNS, LNPSNL)    ! Node set node list array
      CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)    ! Node set df list array
      CALL MDRSRV ('XFACNP', KXFACN, LNPSNL)    ! Expanded df list array
      CALL MDSTAT (NERR, MEM)

      if (numnps .gt. 0) then
         call exgcns (netid, a(kidns), a(knnns), a(kndns), a(kixnns),
     &                a(kixdns), a(klstns), a(kfacns), nerr)
         if (nerr .lt. 0) then
            call exerr('ex2ex1v2', 'Error from exgcns', exlmsg)
            goto 140
         endif
      endif

C     Massage node sets distribution factors to include '1' for node sets
C       without Dfs by walking KNDNS array, checking for 0, and filling where
C       necessary.

      do 64 i=0, numnps-1
        if (ia(kndns+i) .eq. 0) then
          do 60 ii=0, ia(knnns+i)-1
            a(kxfacn+ia(kixnns+i)-1+ii) = 1.0! Force unity distribution factor
60        continue
        else
          do 62 ii=0, ia(kndns+i)-1
            a(kxfacn+ia(kixnns+i)-1+ii) = a(kfacns+ia(kixdns+i)-1+ii)
62        continue
        endif
64      continue

      CALL DBONPS (NDB, NUMNPS, LNPSNL,
     &   A(KIDNS), A(KNNNS), A(KIXNNS), A(KLSTNS), A(KXFACN))

      CALL MDDEL ('IDNPS')
      CALL MDDEL ('NNNPS')
      CALL MDDEL ('NDNPS')
      CALL MDDEL ('IXNNPS')
      CALL MDDEL ('IXDNPS')
      CALL MDDEL ('LSTNPS')
      CALL MDDEL ('FACNPS')
      CALL MDDEL ('XFACNP')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

C   --Read the side sets

      CALL MDRSRV ('IDESS', KIDSS, NUMESS)      ! side set id array
c     write(*,*)'side set id array size: ',numess
      CALL MDRSRV ('NEESS', KNESS, NUMESS)      ! number of ss elems array
c     write(*,*)'number of side set elements  array size: ',numess
      CALL MDRSRV ('NDESS', KNDSS, NUMESS)      ! number of dist factors array
c     write(*,*)'number of dist factors array size: ',numess
      CALL MDRSRV ('NNESS', KNNSS, NUMESS)      ! number of nodes array
c     write(*,*)'number of side set nodes array size: ',numess
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)    ! index into elements array
c     write(*,*)'index into side set elements array size: ',numess
      CALL MDRSRV ('IXDESS', KIXDSS, NUMESS)    ! index into dist factors array
c     write(*,*)'index into side set dist factors  array size: ',numess
      CALL MDRSRV ('IXNESS', KIXNSS, NUMESS)    ! index into nodes array
c     write(*,*)'index into side set nodes array size: ',numess
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)    ! element list
c     write(*,*)'side set element list array size: ',lessel
      CALL MDRSRV ('LTNESS', KLTNSS, LESSNL)    ! node list (21 is max possible)
c     write(*,*)'side set node list array size: ',lessnl
      CALL MDRSRV ('LTNNSS', KLTNNS, LESSEL)    ! node count array
c     write(*,*)'side set node count array size: ',lessel
      CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)    ! side list
c     write(*,*)'side set side list array size: ',lessel
      CALL MDRSRV ('FACESS', KFACSS, LESSDL)    ! dist factors list
c     write(*,*)'side set dist factors list array size: ',lessdl
      CALL MDRSRV ('XFACES', KXFACS, LESSNL)    ! dist factors list(w/all DF)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      if (numess .gt. 0) then
        call exgcss (netid, a(kidss), a(kness), a(kndss),
     &                a(kixess), a(kixdss),
     &                a(kltess), a(kltsss), a(kfacss), nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgcss', exlmsg)
           goto 140
        endif

C Convert sides to nodes

        isoff = 0               ! offset into element list for current side set
        nodcnt = 0              ! node count for current side set
        do 104 i=0,numess-1     ! loop through ss elem blks

          ia(kixnss+i)=nodcnt+1                         ! update index array

          call exgsp(netid,ia(kidss+i),nsess,ndess,nerr)! get num of sides & df
          if (nerr .lt. 0) then
             call exerr('ex2ex1v2', 'Error from exgsp', exlmsg)
             goto 140
          endif
c         write(*,*)'SS ID: ',ia(kidss+i)
c         write(*,*)' # of sides: ',nsess
c         write(*,*)' # of dist factors: ',ndess

          call exgssn(netid,ia(kidss+i),a(kltnns+isoff),
     &               a(kltnss+nodcnt),nerr)             ! get side set nodes
          if (nerr .lt. 0) then
             call exerr('ex2ex1v2', 'Error from exgssn', exlmsg)
             goto 140
          endif
          nness = 0
          do 102 ii=0,nsess-1                           ! sum node counts to
            nness=nness+ia(kltnns+isoff+ii)             ! calculate next index
102       continue
c         write(*,*)' # of nodes: ',nness
          ia(knnss+i)=nness
          nodcnt=nodcnt+nness
          isoff=isoff+nsess
104     continue
      endif

C     Massage side sets distribution factors to include '1' for side sets
C       without Dfs by walking KNDSS array, checking for 0, and filling where
C       necessary.

      do 110 i=0, numess-1
        if (ia(kndss+i) .eq. 0) then
          do 106 ii=0, ia(knnss+i)-1
            a(kxfacs+ia(kixnss+i)-1+ii) = 1.0! Force unity distribution factor
106       continue
        else
          do 108 ii=0, ia(knnss+i)-1
            a(kxfacs+ia(kixnss+i)-1+ii) = a(kfacss+ia(kixdss+i)-1+ii)
108       continue
        endif
110     continue

      CALL DBOESS (NDB, NUMESS, LESSEL, LESSNL,
     &   A(KIDSS), A(KNESS), A(KNNSS), A(KIXESS), A(KIXNSS),
     &   A(KLTESS), A(KLTNSS), A(KXFACS))

      CALL MDDEL ('IDESS')
      CALL MDDEL ('NEESS')
      CALL MDDEL ('NDESS')
      CALL MDDEL ('NNESS')
      CALL MDDEL ('IXEESS')
      CALL MDDEL ('IXDESS')
      CALL MDDEL ('IXNESS')
      CALL MDDEL ('LTEESS')
      CALL MDDEL ('LTNESS')
      CALL MDDEL ('LTNNSS')
      CALL MDDEL ('LTSESS')
      CALL MDDEL ('FACESS')
      CALL MDDEL ('XFACES')

C   --Read the QA records

      nqarec = 0
      call exinq (netid, EXQA, nqarec, r, name, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exinq', exlmsg)
         goto 140
      endif

      if (nqarec .gt. 0 .and. nqarec .le. MAXQA) then
        call mcrsrv('QARECS', kqarec, 4*nqarec*8)
        call mcrsrv('QATMP', kqatmp, 4*nqarec*mxstln)
        call mcstat(nerr, mem)
        if (nerr .ne. 0) goto 130
      else
        kqarec = 1
      end if
      if (nqarec .gt. MAXQA) nqarec = 0

      ninfo = 0
      call exinq (netid, EXINFO, ninfo, r, name, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exinq', exlmsg)
         goto 140
      endif

      if (ninfo .gt. 0 .and. ninfo .le. MAXINF) then
        call mcrsrv('INFO', kinfo, ninfo*mxlnln)
        call mcstat(nerr, mem)
        if (nerr .ne. 0) goto 130
      else
        kinfo = 1
      end if
      if (ninfo .gt. MAXINF) ninfo = 0

      call rdqain (netid, nqarec, c(kqatmp), ninfo, c(kinfo))

      if (nqarec .gt. 0)
     &    call resize (nqarec, c(kqarec), c(kqatmp))

      IF (NQAREC .GE. 0) THEN
         CALL DBOQA (NDB, NQAREC, c(kqarec), NINFO, c(kinfo))
      END IF

C   --Read in the number of element variable names

      call exgvp (netid, 'e', nvarel, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgvp', exlmsg)
         goto 140
      endif

C   --Read in the number of global variable names

      call exgvp (netid, 'g', nvargl, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgvp', exlmsg)
         goto 140
      endif

C   --Read in the number of nodal variable names

      call exgvp (netid, 'n', nvarnp, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgvp', exlmsg)
         goto 140
      endif

      nvarhi = 0

      call mdrsrv ('ISEVOK', kievok, nvarel*nelblk)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

c       read in the element variable truth table

      if (nvarel .gt. 0) then
        call exgvtt (netid, nelblk, nvarel, a(kievok), nerr)
        if (nerr .gt. 0) then
          if (nvarel .gt. 0) then
            write (*,'(4x,"must have element variable truth table")')
            goto 140
          endif
        endif
        if (nerr .lt. 0) then
          call exerr('ex2ex1v2', 'Error from exgvtt', exlmsg)
          goto 140
        endif
      end if

c       read in the element variable names

      ixev = 1
      if (nvarel .gt. 0) then
        call exgvan (netid, 'e', nvarel,mames(ixev), nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgvan', exlmsg)
           goto 140
        endif
      end if

c       read in the global variable names

      ixgv = ixev + nvarel
      if (nvargl .gt. 0) then
        call exgvan (netid, 'g', nvargl,mames(ixgv), nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgvan', exlmsg)
           goto 140
        endif
      end if

c       read in the nodal variable names

      ixnv = ixgv + nvargl
      if (nvarnp .gt. 0) then
        call exgvan (netid, 'n', nvarnp, mames(ixnv), nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgvan', exlmsg)
           goto 140
        endif
      end if

c       read in the history variable names

      ixhv = ixnv + nvarnp

c       read coordinate names

      call exgcon (netid, mameco, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgcon', exlmsg)
         goto 140
      endif

      CALL DBPINI ('V', NTXT, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL,
     &      NVARHI, NVARGL, NVARNP, NVAREL)

      do 111 i=1, ndim
        nameco(i) = mameco(i)(:8)
 111  continue
      do 112 i=1, (nvarhi+nvargl+nvarnp+nvarel)
        names(i) = mames(i)(:8)
 112  continue

      CALL DBONAM (NDB, NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
     &   nameco, namelb,
     &   names(ixhv), names(ixgv), names(ixnv), names(ixev),
     &   A(KIEVOK))

      CALL MDRSRV ('VARHI', KVARHI, NVARHI)
      CALL MDRSRV ('VARGL', KVARGL, NVARGL)
      CALL MDRSRV ('VARNP', KVARNP, NVARNP * NUMNP)
      CALL MDRSRV ('VAREL', KVAREL, NVAREL * NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

c       read in the number of history time steps and the number of
C       whole time steps

      call exinq (netid, EXTIMS, ntime, s, name, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exqini', exlmsg)
         goto 140
      endif
      if (ntime .eq. 0) then
        write(errmsg,'("GENISIS file - no time steps written")')
        call exerr('ex2ex1v2', errmsg, EXPMSG)
        goto 140
      endif
      numstp = ntime

c       read the time step information

      istep = 0
      call exgtim(netid, istep+1, wtime, nerr)
      if (nerr .lt. 0) then
         call exerr('ex2ex1v2', 'Error from exgtim', exlmsg)
         goto 140
      endif

      do 300 ihstep=1,numstp

        write (*,'(4x,"processing time step ", i4)') ihstep

c         get history information

        whotim = .true.
        call exgtim(netid, ihstep, wtime, nerr)
        if (nerr .lt. 0) then
           call exerr('ex2ex1v2', 'Error from exgtim', exlmsg)
           goto 140
        endif
        htime = wtime

c          If a whole time step, do global, nodal, and element
c          variables for the time step.

        if ((whotim) .or. (wtime .eq. htime)) then

          whotim =.true.
          istep = istep + 1

c           get the global variable values

          if( nvargl .gt. 0) then
            call exggv (netid, istep, nvargl, a(kvargl), nerr)
            if (nerr .lt. 0) then
               call exerr('ex2ex1v2', 'Error from exggv', exlmsg)
               goto 140
            endif
          end if

c           get the nodal variable values

          do 210 j=1, nvarnp
            call exgnv (netid, istep, j, numnp,
     &         a(kvarnp+(j-1)*numnp), nerr)
            if (nerr .lt. 0) then
               call exerr('ex2ex1v2', 'Error from exgnv', exlmsg)
               goto 140
            endif
210       continue

c           get element variable values

          if (nvarel .gt. 0) then
            ielo=0
            do 250 k = 1,nelblk
              l=(k-1)*nvarel
              do 240 j=1, nvarel

c                If truth table indicates element values are available
c                for the element variable, get the values for the
c                element variable.

                if(a(kievok+l +j-1) .ne. 0) then
                  call exgev (netid, istep, j, a(kidelb+k-1),
     &                 a(knelb+k-1), a(kvarel+ielo), nerr)
                  if (nerr .lt. 0) then
                     call exerr('ex2ex1v2', 'Error from exgev', exlmsg)
                     goto 140
                  endif
                  call getin (ia(knelb+k-1),num)
                  ielo = ielo+num
                end if
240           continue
250         continue
          end if
        else
          whotim=.false.
        end if
        CALL DBOSTE (NDB, ihstep, NVARHI, NVARGL, NVARNP, NUMNP,
     &      NVAREL, NELBLK, a(knelb), a(kievok),
     &      HTIME, WHOTIM, A(KVARHI), A(KVARGL), A(KVARNP),
     &      A(KVAREL))

300   continue

      call MDDEL ('IDELB')
      CALL MDDEL ('VARHI')
      CALL MDDEL ('VARGL')
      CALL MDDEL ('VARNP')
      CALL MDDEL ('VAREL')
      CALL MDDEL ('NUMELB')

      CALL INTSTR (1, 0, IHSTEP-1, STR8, LSTR)
      WRITE (*, 10010) STR8(:LSTR)
10010  FORMAT (/, 4X, A,
     &   ' time steps have been written to the database')

      GOTO 140

  130 CONTINUE
      CALL MEMERR
      GOTO 140

  140 CONTINUE

c       close all files

      CLOSE (NDB, IOSTAT=IDUM)

      if (netid .ge. 0 ) call exclos (netid, ierr)

      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))

      END

      subroutine mlist()
      call mdlist(6)
      return
      end

      subroutine rdqain (ndb, nqarec, qarec, ninfo, info)
      include 'exodusII.inc'
      integer ndb
      character*(32) qarec(4,nqarec)
      character*(80) info(ninfo)

      if (nqarec .gt. 0) then
         call exgqa (ndb, qarec, nerr)
         if (nerr .lt. 0) then
            call exerr('ex2ex1v2', 'Error from exgqa', exlmsg)
         endif
      endif
      if (ninfo .gt. 0) then
         call exginf (ndb, info, nerr)
         if (nerr .lt. 0) then
            call exerr('ex2ex1v2', 'Error from exginf', exlmsg)
         endif
      endif

      return
      end

C=======================================================================
      SUBROUTINE RESIZE (NQAREC, QAREC, QATMP)
C=======================================================================
C   --
C   --RESIZE - resizes the qa records from length 32 to 8
C   --
C   --Parameters:
C   --   NQAREC - IN - the number of QA records
C   --   QAREC  - IN - the QA records containing size = 8
C   --   QATMP  - IN - the QA records containing size = 32

      INTEGER NQAREC
      CHARACTER*8 QAREC(4,NQAREC)
      CHARACTER*32 QATMP(4,NQAREC)

      IF (NQAREC .GT. 0) THEN
         DO 50 I = 1, NQAREC
            DO 75 J = 1, 4
               QAREC(J,I) = QATMP(J,I)(:8)
 75         CONTINUE
 50      CONTINUE
      END IF

      RETURN
      END

      subroutine getin (num1,num2)
      integer num1(*)
      num2=num1(1)
      return
      end
