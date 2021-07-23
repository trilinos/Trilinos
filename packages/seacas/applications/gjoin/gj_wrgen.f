C Copyright(C) 1999-2021 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRGEN (A,IA, FILNAM, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSDL,
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, KFACSS,
     &   kltsss, NQAREC, QAREC, NINFO, INFREC, NAMELB, L64BIT, NC4,
     $   NAMBK, NAMNS, NAMSS, *)
C=======================================================================

C   --*** WRGEN *** (GJOIN) Writes the GENESIS database
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --WRGEN writes the GENESIS database.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   FILNAM - IN - the database filename
C   --   TITLE - IN - the database title
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   NUMEL - IN - the number of elements
C   --   NELBLK - IN - the number of element blocks
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN - the length of the nodal point sets node list
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the element side sets element list
C   --   LESSNL - IN - the length of the element side sets node list
C   --   KXN, KYN, KZN - IN - index of XN, YN, ZN; nodal coordinates
C   --   KMAPEL - IN - index of MAPEL; the element order map
C   --   KIDELB - IN - index of IDELB; the element block IDs for each block
C   --   KNELB - IN - index of NUMELB; the number of elements in each block
C   --   KNLNK - IN - index of NUMLNK; the number of nodes per element
C   --      in each block
C   --   KNATR - IN - index of NUMATR; the number of attributes in each block
C   --   KLINK - IN - index of LINK; the connectivity for each block
C   --   KATRIB - IN - index of ATRIB; the attributes for each block
C ... Nodesets
C   --   KIDNS - IN - index of IDNPS; the nodal point set ID for each set
C   --   KNNNS - IN - index of NNNPS; the number of nodes for each set
C   --   KIXNNS - IN - index of IXNNPS; the index of the first node for each set
C   --   KLTNNS - IN - index of LTNNPS; the nodes for all sets
C   --   KFACNS - IN - index of FACNPS; the distribution factors for all sets
C .... Sidesets
C   --   KIDSS - IN - index of IDESS; the element side set ID for each set
C   --   KNESS - IN - index of NEESS; the number of elements for each set
C   --   KIXESS - IN - index of IXEESS; the index of the first element
C   --      for each set
C   --   KLTESS - IN - index of LTEESS; the elements for all sets
C   --   KFACSS - IN - index of FACESS; the distribution factors for all sets
C   --   NAMELB - IN - the names of the element blocks
C   --   L64BIT - IN - true if use 64-bit integer output database

      include 'exodusII.inc'
      include 'gj_params.blk'
      include 'gj_namlen.blk'

      DIMENSION A(*), IA(*)
      CHARACTER*(*) FILNAM
      character*(MXLNLN) title
      character*(MXSTLN) qarec(4,MAXQA)
      character*(MXLNLN) infrec(MAXINF)
      character*(MXSTLN) nameco(6), namelb(*)
      character*(namlen) nambk(*), namns(*), namss(*)
      LOGICAL            l64bit, NC4

C      --QAREC - the QA records
C      --INFREC - the information records
      integer cpuws,wsout

      DATA NAMECO /'X', 'Y', 'Z', 'RX', 'RY', 'RZ'/

C     Make netCDF and exodus errors not show up

      call exopts (0,ierr)
      wsout = iowdsz()
      write(*,*)'Output word size: ',wsout

C   --Open the database
C     Create the netcdf file

      CALL SQZSTR(FILNAM, LNAM)
      LNAM = LENSTR (FILNAM)

      cpuws = 0
      MODE = EX_CLOBBER
      if (l64bit) then
        MODE = MODE + EX_ALL_INT64_DB + EX_ALL_INT64_API
      end if
      if (nc4) then
        MODE = MODE + EX_NETCDF4
      end if
      idexo = excre (filnam(:lnam), MODE, cpuws, wsout, ierr)
      if (ierr .lt. 0) then
         call exerr('gjoin2', 'Error from excre', exlmsg)
         go to 150
      endif

C   --Write the QA records

      if (nqarec .gt. 0) then
         call expqa (idexo, nqarec, qarec, ierr)
         if (ierr .lt. 0) then
            call exerr ('gjoin2', 'Error from expqa', exlmsg)
            goto 150
         endif
      endif

C   --Write the info records

      if (ninfo .gt. 0) then
         call expinf (idexo, ninfo, infrec, ierr)
         if (ierr .lt. 0) then
            call exerr ('gjoin2', 'Error from expinf', exlmsg)
            goto 150
         endif
      endif

      lnpsdl = 0
      if (numnps .gt. 0) then
         call mdrsrv('NSDF', kansdf, numnps)
         do 90 i=1, numnps
            ia(kansdf+i-1) = ia(knnns+i-1)
            lnpsdl = lnpsdl + ia(kansdf+i-1)
 90      continue
      end if

      IDUM = 0
      CALL DBPINI ('NTIS', idexo, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, lnpsdl, NUMESS, LESSEL, LESSDL,
     &   IDUM, IDUM, IDUM, FILNAM(:LNAM))

C   --Write the initial variables

      call expini (idexo, title, ndim, numnp, numel, nelblk, numnps,
     &    numess, ierr)
      if (ierr .lt. 0) then
         call exerr ('gjoin2', ' Error from expini', exlmsg)
         goto 150
      endif

C   --Write the coordinates

      call expcor (idexo, a(kxn), a(kyn), a(kzn), ierr)
      if (ierr .lt. 0) then
         call exerr ('gjoin2', 'Error from expcor', exlmsg)
         goto 150
      endif

      call expcon (idexo, nameco, ierr)
      if (ierr .lt. 0) then
         call exerr ('gjoin2', 'Error from expcon', exlmsg)
         goto 150
      endif

C   --Write out the nodal point sets

      if (numnps .gt. 0) then
         call expcns (idexo, ia(kidns), ia(knnns), ia(kansdf),
     &                ia(kixnns), ia(kixnns), ia(kltnns),
     &                a(kfacns), ierr)
         if (ierr .lt. 0) then
            call exerr ('gjoin2', 'Error from expcns', exlmsg)
            goto 150
         endif
         call mddel('NSDF')
         call expnams(idexo, 2, numnps, namns, ierr)
      endif

C   --Write element side sets
      if (numess .gt. 0) then
         call expcss (idexo, ia(kidss), ia(kness), ia(kndss),
     &    ia(kixess), ia(kixdss), ia(kltess),
     &    ia(kltsss), a(kfacss), ierr)
        if (ierr .lt. 0) then
          call exerr ('gjoin2', 'Error from expcss', exlmsg)
          goto 150
        endif
         call expnams(idexo, 3, numess, namss, ierr)
      endif

C   --Write the element blocks

C        Write concatenated element block parameters
      call expclb (idexo, ia(kidelb), namelb,
     &  ia(knelb), ia(knlnk), ia(knatr), .FALSE., ierr)
      if (ierr .lt. 0) then
        call exerr('gjoin2', 'Error from expclb', exlmsg)
        goto 150
      endif

      ioff = katrib
      iptr = klink
      do 100 ielb = 1, nelblk

C        Write block attributes

         if (ia(knatr+ielb-1) .gt. 0) then
            call expeat (idexo, ia(kidelb+ielb-1), a(ioff), ierr)
            if (ierr .lt. 0) then
               call exerr ('gjoin2', 'Error from expeat', exlmsg)
               goto 150
            endif
         endif

C        Write the element block connectivity,
C        skipping null element blocks

         if (ia(knelb+ielb-1) .eq. 0) then
            write(*,*)'Null element block: ',ielb
         else
            call expelc (idexo, ia(kidelb+ielb-1), ia(iptr), ierr)
            if (ierr .lt. 0) then
               call exerr ('gjoin2', 'Error from expelc', exlmsg)
               goto 150
            endif
         endif

         ioff = ioff + ( ia(knatr+ielb-1) * ia(knelb+ielb-1) )
         iptr = iptr + ( ia(knlnk+ielb-1) * ia(knelb+ielb-1) )
 100  continue
      call expnams(idexo, 1, nelblk, nambk, ierr)

 150  call exclos (idexo, ierr)

      RETURN
      END
