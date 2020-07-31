C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      PROGRAM EX1EX2V2
C=======================================================================

C   --*** EX1EX2V2 *** EXODUS I to EXODUS II translator
C   --
C   --EX1EX2V2 reads EXODUS I database and writes an EXODUS II V2.03 database
C   --

      include 'exodusII.inc'
      INCLUDE 'argparse.inc'

      CHARACTER*8 QAINFO(6)

      CHARACTER*80 TITLE

      PARAMETER (MAXQA = 100, MAXINF = 100)
      PARAMETER (MAXDIM=6, MAXELB=512, MAXVAR=512)
      CHARACTER*8 QAREC(4,MAXQA)
      CHARACTER*80 INFO(MAXINF)
      CHARACTER*8 NAMECO(MAXDIM)
      CHARACTER*8 NAMELB(MAXELB)
      CHARACTER*8 NAMES(MAXVAR)
      character*2048 exofil, netfil, scratch
      character*8 name

      integer cpuws,wsout

      DIMENSION A(1)
      DIMENSION IA(1)
      EQUIVALENCE (A(1), IA(1))
C      --A - the dynamic memory base array

      CHARACTER*8 STR8
      LOGICAL EXODUS
      LOGICAL WHOTIM

      data (qainfo(i), i=1,3) / 'ex1ex2v2', '20110616', 'v 2.11  ' /
      data iin,iout/5,6/
      data nsteps /0/
      data cpuws,wsout /0,0/

      CALL STRTUP (QAINFO)

      CALL BANNER (0, QAINFO,
     &   'EXODUS I TO EXODUS II FILE TRANSLATOR',
     &   ' ', ' ')
      call exinq (netid, EXLBVR, idummy, exlibversion, name, nerr)
      write(*,'(A,F6.3)')' ExodusII Library version ',
     1          exlibversion

C   --Open the input and output files

      NDB = 11
      NET = 20

c       make netCDF and exodus errors not show up

      call exopts(0,ierr)

C .. Get filename from command line.  If not specified, emit error message
      NARG = argument_count()
      if (narg .lt. 2) then
        CALL PRTERR ('FATAL', 'Filenames not specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "ex1ex2v2 exo1_file exo2_file"')
        GOTO 140
      else if (narg .gt. 2) then
        CALL PRTERR ('FATAL', 'Too many arguments specified.')
        CALL PRTERR ('FATAL',
     *    'Syntax is: "ex1ex2v2 exo1_file exo2_file"')
        GOTO 140
      end if

      CALL get_argument(1,exofil, lnam)
      write(*,*)'Input filename: ',exofil(1:lnam)
      open(unit=ndb, file=exofil(:lnam), form='unformatted',
     *     status='old', iostat=ierr)
       IF (IERR .NE. 0) THEN
         SCRATCH = 'Database "'//exofil(:lnam)//'" does not exist.'
         CALL PRTERR ('FATAL', SCRATCH(:LENSTR(SCRATCH)))
         GOTO 140
      END IF

      CALL get_argument(2,netfil, lnam)
      write(*,*)'Output filename: ',netfil(1:lnam)

      wsout = 8
      write(*,*)'Output word size: ',wsout
      CALL MDINIT (A)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130
      CALL MDFILL(-1)

C   --Read the initial variables from exodusI database

      CALL DBIINI (NDB, '*', NVERS, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, *150)

c     create the a netcdf file
      idexo = excre (netfil(1:lnam), EXCLOB, cpuws, wsout, ierr)
      if (ierr .lt. 0) then
        call exerr('ex1ex2v2','Error from excre', EXLMSG)
        go to 140
      end if

c     write initial variables to netcdf file
      call expini (idexo, title, ndim, numnp, numel, nelblk, numnps,
     &    numess, ierr)
      if (ierr .lt. 0) then
        call exerr ('ex1ex2v2',' Error from expini', EXLMSG)
        goto 150
      end if

      CALL DBPINI ('NTIS', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL,
     &   IDUM, IDUM, IDUM, IDUM)
      WRITE(*,*)

C   --Read the coordinates from the exodusI database

      CALL MDRSRV ('XN', KXN, NUMNP)
      CALL MDRSRV ('YN', KYN, NUMNP)
      IF (NDIM .GE. 3) THEN
        CALL MDRSRV ('ZN', KZN, NUMNP)
      ELSE
        KZN = 1
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL DBIXYZ (NDB, '*', NDIM, NUMNP, A(KXN), A(KYN), A(KZN), *150)

c     write the coordinates to the regular netcdf file

      call expcor (idexo, a(kxn), a(kyn), a(kzn), ierr)

      if (ierr .lt. 0) then
        call exerr ('ex1ex2v2','Error from expcor', EXLMSG)
        goto 150
      end if
      CALL MDDEL ('XN')
      CALL MDDEL ('YN')
      IF (NDIM .GE. 3) CALL MDDEL ('ZN')
C   --Read the element order map from the exodusI database

      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL DBIMAP (NDB, '*', NUMEL, IA(KMAPEL), *150)

c     write the element order map to the regular netcdf file

      call expmap (idexo, ia(kmapel), ierr)
      if (ierr .lt. 0) then
        call exerr ('ex1ex2v2','Error from expmap', EXLMSG)
        goto 150
      end if
      CALL MDDEL ('MAPEL')

C   --Read the element blocks

      CALL MDRSRV ('NUMELB', KNELB, NELBLK)
      CALL MDRSRV ('LINK', KLINK, 0)
      CALL MDRSRV ('ATRIB', KATRIB, 0)
      call MDRSRV ('IDELB', KIDELB, NELBLK)
      call MDRSRV ('NUMLNK', KNMLNK, NELBLK)
      call MDRSRV ('NUMATR', KNMATR, NELBLK)
      CALL MDSTAT (IERR, MEM)
      IF (IERR .GT. 0) goto 150

      CALL DBIELB (NDB, '*', 1, nelblk, IA(KIDELB), IA(KNELB),
     *  IA(KNMLNK), IA(KNMATR), A, KLINK, KATRIB, *150)
      CALL MDSTAT (IERR, MEM)
      IF (IERR .GT. 0) goto 150

C   --Read the nodal points sets

      CALL MDRSRV ('IDNPS',  KIDNS, NUMNPS)
      CALL MDRSRV ('NNNPS', KNNNS, NUMNPS)
      CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
      CALL MDRSRV ('LSTNPS', KLSTNS, LNPSNL)
      CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL DBINPS (NDB, '*', NUMNPS, LNPSNL,
     &   IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLSTNS), A(KFACNS), *150)

C   --Read the element side sets

      CALL MDRSRV ('IDESS', KIDSS, NUMESS)
      CALL MDRSRV ('NEESS', KNESS, NUMESS)
      CALL MDRSRV ('NNESS', KNNSS, NUMESS)
      CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
      CALL MDRSRV ('IXNESS', KIXNSS, NUMESS)
      CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
      CALL MDRSRV ('LTNESS', KLTNSS, LESSNL)
      CALL MDRSRV ('FACESS', KFACSS, LESSNL)
      CALL MDRSRV ('SACESS', KLTSSS, LESSEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      CALL DBIESS (NDB, '*', NUMESS, LESSEL, LESSNL,
     &  IA(KIDSS), IA(KNESS), IA(KNNSS), IA(KIXESS), IA(KIXNSS),
     &  IA(KLTESS), IA(KLTNSS), A(KFACSS), *150)

C   --Read the QA and info records
C ... Exodus set to .FALSE. if end of file during this read
      CALL DBIQA (NDB, '*', MAXQA, MAXINF, NQAREC, QAREC, NINFO, INFO,
     &   EXODUS, *150)

c**********************************************************************
C   --Read the database names

      if (exodus) then
C ... Exodus set to .FALSE. if end of file during this read
        CALL DBINAM (NDB, '*', NDIM, NELBLK,
     &    NNDIM, NNELB, NVARHI, NVARGL, NVARNP, NVAREL,
     &    NAMECO, NAMELB, NAMES, IXHV, IXGV, IXNV, IXEV,
     &    A, KIEVOK, EXODUS, *150)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 130

        if (.not. exodus) then
          nvarhi = 0
          nvargl = 0
          nvarnp = 0
          nvarel = 0
        end if
        CALL DBPINI ('V', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &    NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL,
     &    NVARHI, NVARGL, NVARNP, NVAREL)
      else
C ... Assign dummy coordinate names...
        NAMECO(1) = 'X'
        NAMECO(2) = 'Y'
        NAMECO(3) = 'Z'
C ... Try to infer the element block names also
        do 50 ielb = 0, nelblk-1
          namelb(ielb+1) = 'UNKNOWN'
          npe = ia(knmlnk+ielb)
          if (ndim .eq. 2) then
            if (npe .eq. 4) then
              namelb(ielb+1) = 'QUAD'
            else if (npe .eq. 8) then
              namelb(ielb+1) = 'QUAD8'
            else if (npe .eq. 3) then
              namelb(ielb+1) = 'TRIANGLE'
            end if
          else if (ndim .eq. 3) then
            if (npe .eq. 4) then
              if (ia(knmatr+ielb) .eq. 0) then
                namelb(ielb+1) = 'TETRA'
              else
                namelb(ielb+1) = 'SHELL'
              end if
            else if (npe .eq. 8) then
              if (ia(knmatr+ielb) .eq. 0) then
                namelb(ielb+1) = 'HEX'
              else
                namelb(ielb+1) = 'SHELL8'
              end if
            else if (npe .eq. 6) then
              namelb(ielb+1) = 'WEDGE'
            end if
          end if
 50     continue
      END IF

c*********
      ioff = 0
      DO 100 IELB = 1, NELBLK

c          write element block parameters to the netcdf file

         call expelb (IDEXO, IA(KIDELB+IELB-1), namelb(IELB),
     1      IA(KNELB+IELB-1),
     2      IA(KNMLNK+IELB-1), IA(KNMATR+IELB-1), IERR)
         IF (IERR .lt. 0) THEN
                CALL exerr('ex1ex2v2','Error from expelb',EXLMSG)
         ENDIF

c          write block attributes to the netcdf file

         IF (IA(KNMATR+IELB-1) .GT. 0) THEN
           call expeat (IDEXO, IA(KIDELB+IELB-1), A(KATRIB+ioff), IERR)
           IF (IERR .lt. 0) THEN
                CALL exerr ('rdelb','Error from expeat', EXLMSG)
           ENDIF
         end if

         ioff = ioff + ia(knmatr + ielb-1) * ia(knelb + ielb-1)

c         CALL MDLONG ('LINK', KLINK, 0)
c         CALL MDLONG ('ATRIB', KATRIB, 0)
  100 CONTINUE

      iptr = klink
      do 101 ielb = 1, nelblk

c          write the element block connectivity to the netcdf file
c            skipping null element blocks

        if (IA(KNELB+IELB-1) .eq. 0) then
          write(*,*)'Null element block: ',ielb
        else
          call expelc (idexo, ia(kidelb+ielb-1), ia(iptr), ierr)
          if (ierr .lt. 0) then
            call exerr ('ex1ex2v2','Error from expelc', exlmsg)
            goto 150
          end if
        end if

        iptr = iptr + ( ia(knmlnk+ielb-1) * ia(knelb+ielb-1) )

101   continue

c     write out the nodal point sets to the regular netcdf file
c       Note: For exodus I data, dist factors always exist.

      if (numnps .gt. 0) then
        call expcns (idexo, ia(kidns), ia(knnns), ia(knnns),
     *    ia(kixnns), ia(kixnns), ia(klstns), a(kfacns), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expcns', exlmsg)
          goto 150
        end if
      endif

c     write element side sets

c       Note: Exodus II V2.0 represents a major change for side sets:
c               They are represented as side IDs - not node IDs and
c               must be translated.
      if (numess .gt. 0) then
        call excn2s (idexo, ia(kness), ia(knnss), ia(kixess),
     1          ia(kixnss), ia(kltess), ia(kltnss), ia(kltsss), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from excn2s', exlmsg)
          goto 150
        end if

        call expcss (idexo, ia(kidss), ia(kness), ia(knnss), ia(kixess),
     &      ia(kixnss), ia(kltess), ia(kltsss), a(kfacss), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expcss', exlmsg)
          goto 150
        end if
      endif

      call mddel ('LINK')
      call mddel ('NUMLNK')
      CALL MDDEL ('ATRIB')
      CALL MDDEL ('NUMATR')

      CALL MDDEL ('IDNPS')
      CALL MDDEL ('NNNPS')
      CALL MDDEL ('IXNNPS')
      CALL MDDEL ('LSTNPS')
      CALL MDDEL ('FACNPS')

      CALL MDDEL ('IDESS')
      CALL MDDEL ('NEESS')
      CALL MDDEL ('NNESS')
      CALL MDDEL ('IXEESS')
      CALL MDDEL ('IXNESS')
      CALL MDDEL ('LTEESS')
      CALL MDDEL ('LTNESS')
      CALL MDDEL ('FACESS')
      CALL MDDEL ('SACESS')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

c      write the QA records

      IF (NQAREC .GT. 0) then
        call expqa (idexo, NQAREC, QAREC, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expqa', exlmsg)
          goto 150
        end if
      end if

c       write the info records

      if (NINFO .gt. 0) then
        call expinf (idexo, ninfo, info, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expinf', exlmsg)
          goto 150
        end if
      end if

c**********************************************************************

c        write coordinate names

      call expcon (idexo, nameco, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expcon', exlmsg)
          goto 150
        end if

      if (.not. EXODUS) goto 150

c       write the number of global variables

      if (nvargl .gt. 0) then
        call expvp (idexo, 'G', nvargl, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvp', exlmsg)
          goto 140
        end if

c       write the global variable names

        call expvan (idexo, 'G', nvargl, names(ixgv), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvan', exlmsg)
          goto 140
        end if
      end if

c       write the number of nodal variables

      if (nvarnp .gt. 0) then
        call expvp (idexo, 'N', nvarnp, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvp', exlmsg)
          goto 140
        end if

c       write the nodal variable names

        call expvan (idexo, 'N', nvarnp, names(ixnv), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvan', exlmsg)
          goto 140
        end if
      end if

c       write the number of element variables

      if (nvarel .gt. 0) then
        call expvp (idexo, 'E', nvarel, ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvp', exlmsg)
          goto 140
        end if

c       write the element variable names

        call expvan (idexo, 'E', nvarel, names(ixev), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from exvan', exlmsg)
          goto 140
        end if
      end if

c       write the element variable truth table

      call mdrsrv ('ebids', kebids, nelblk)
      call exgebi (idexo, ia(kebids), ierr)
      if (ierr .lt. 0) then
        call exerr ('ex1ex2v2','Error from exgebi', exlmsg)
        goto 140
      end if

      if (nvarel .gt. 0) then
        call expvtt(idexo, nelblk, nvarel, ia(kievok), ierr)
        if (ierr .lt. 0) then
          call exerr ('ex1ex2v2','Error from expvtt', exlmsg)
          goto 140
        end if
      end if

      call mddel ('ebids')

      IF (.NOT. EXODUS) GOTO 140

C   --Read the database time steps

      CALL MDRSRV ('VARHI', KVARHI, NVARHI)
      CALL MDRSRV ('VARGL', KVARGL, NVARGL)
      CALL MDRSRV ('VARNP', KVARNP, NVARNP * NUMNP)
      CALL MDRSRV ('VAREL', KVAREL, NVAREL * NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 130

      nwstep = 0
      nhstep = 0

  110 CONTINUE
      IF (.TRUE.) THEN

         CALL DBISTE (NDB, '*', NhSTEP+1,
     &      NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK,
     &      ia(knelb), IA(KIEVOK), TIME, WHOTIM, A(KVARHI),
     &      A(KVARGL), A(KVARNP), A(KVAREL), *120)

         nhstep = nhstep + 1
         if (whotim) then
           nwstep = nwstep+1
           write (*,'(4x, "processing whole time step ", i4)') nwstep

c            write global variables

           if (nvargl .gt. 0) then
             call expgv (idexo, nwstep, nvargl, a(kvargl), ierr)
             if (ierr .lt. 0) then
               call exerr ('ex1ex2v2','Error from expgv', exlmsg)
               goto 140
             end if
           end if

c            write nodal variable values

           if (nvarnp .gt. 0) then
             do 111 i= 1,nvarnp
               call expnv (idexo, nwstep, i, numnp,
     &           a(kvarnp+(i-1)*numnp), ierr)
               if (ierr .lt. 0) then
                 call exerr ('ex1ex2v2','Error from expnv', exlmsg)
                 goto 140
               end if
111          continue
           end if

c            write element variable values

           if (nvarel .gt. 0) then
             call putev (idexo, nwstep, nelblk, nvarel,
     &         ia(knelb), a(kvarel), ia(kidelb), ia(kievok), ierr)
               if (ierr .lt. 0) then
                 call exerr ('ex1ex2v2','Error from putev', exlmsg)
                 goto 140
               end if
           end if

c            write whole time step

           call exptim (idexo, nwstep, time, ierr)
           if (ierr .lt. 0) then
             call exerr ('ex1ex2v2','Error from exptim', exlmsg)
             goto 140
           end if

         end if
         GOTO 110
      END IF

  120 CONTINUE

      call mddel ('IDELB')

      WRITE (STR8, '(I8)', IOSTAT=K) NwSTEP
      CALL SQZSTR (STR8, LSTR)
      WRITE (*, 10010) STR8(:LSTR)
10010  FORMAT (/, 4X, A,
     &   ' time steps have been written to the file')

      GOTO 140

  130 CONTINUE
      CALL MEMERR
      CALL MDSTAT (NERR, MEM)
      WRITE(*,10020) QAINFO(1), MEM
10020 FORMAT(1x,A8,' currently using ', I10, ' words of memory')

  140 CONTINUE

  150 call exclos (idexo, ierr)

      IF (NDB .NE. 0) CLOSE (NDB, IOSTAT=K)

      call addlog (QAINFO(1)(:lenstr(QAINFO(1))))
      CALL WRAPUP (QAINFO(1))
      STOP
      END

