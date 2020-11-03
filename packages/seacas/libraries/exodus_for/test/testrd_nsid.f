C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      program testrd

c This is a test program for the Fortran binding of the EXODUS II
c database read routines

c       09/07/93 V.R. Yarberry - Modified for API 2.00
      implicit none

      include 'exodusII.inc'

      integer iin, iout, ierr, ioff
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer num_node_sets
      integer num_side_sets
      integer i, j, k, elem_map(100), connect(100), nnpe(10)
      integer ids(10)
      integer num_elem_in_block(10), num_nodes_per_elem(10)
      integer num_attr(10)
      integer cpu_ws,io_ws, mod_sz

      real x(100), y(100), z(100)
      real vers

      character*(MXSTLN) coord_names(3)
      character*(MXLNLN) titl
      character*(MXSTLN) eltype(10)

      data iin /5/, iout /6/

c open EXODUS II files

      cpu_ws = 0
      io_ws = 0

      exoid = exopen ("test-nsided.exo", EXREAD, cpu_ws, io_ws,
     *  vers, ierr)
      write (iout, '(/"after exopen, error = ",i3)')
     1                  ierr

      write (iout, '("test-nsided.exo is an EXODUSII file; version ",
     1                f4.2)') vers
      write (iout, '("  I/O word size",i2)') io_ws

      mod_sz = exlgmd(exoid)
      write (iout, '("  Model Size",i2)') mod_sz

c read database parameters

      call exgini (exoid, titl, num_dim, num_nodes, num_elem,
     1             num_elem_blk, num_node_sets, num_side_sets, ierr)
      write (iout, '(/"after exgini, error = ", i3)' ) ierr

      write (iout, '("database parameters:"/
     1               "title = ", a81 /
     2               "num_dim = ", i3 /
     3               "num_nodes = ", i3 /
     4               "num_elem = ", i3 /
     5               "num_elem_blk = ", i3 /
     6               "num_node_sets = ", i3 /
     7               "num_side_sets = ", i3)')
     8               titl,num_dim, num_nodes, num_elem,
     9               num_elem_blk,num_node_sets, num_side_sets

c read nodal coordinates values and names from database

      call exgcor (exoid, x, y, z, ierr)
      write (iout, '(/"after exgcor, error = ", i3)' ) ierr

      write (iout, '("x coords = ")')
      do 10 i = 1, num_nodes
         write (iout, '(f5.1)') x(i)
10    continue

      write (iout, '("y coords = ")')
      do 20 i = 1, num_nodes
         write (iout, '(f5.1)') y(i)
20    continue

      if (num_dim .gt. 2) then
      write (iout, '("z coords = ")')
      do 22 i = 1, num_nodes
         write (iout, '(f5.1)') z(i)
22    continue
      endif

      call exgcon (exoid, coord_names, ierr)
      write (iout, '(/"after exgcon, error = ", i3)' ) ierr

      write (iout, '("x coord name = ", a9)') coord_names(1)
      write (iout, '("y coord name = ", a9)') coord_names(2)

c read element order map

      call exgmap (exoid, elem_map, ierr)
      write (iout, '(/"after exgmap, error = ", i3)' ) ierr

      do 30 i = 1, num_elem
         write (iout, '("elem_map(",i1,") = ", i1)') i, elem_map(i)
30    continue

c read element block parameters

      call exgebi (exoid, ids, ierr)
      write (iout, '(/"after exgebi, error = ", i3)' ) ierr

      do 40 i = 1, num_elem_blk

         call exgelb (exoid, ids(i), eltype(i), num_elem_in_block(i),
     1                num_nodes_per_elem(i), num_attr(i), ierr)
         write (iout, '(/"after exgelb, error = ", i3)' ) ierr

         write (iout, '("element block id = ", i2,/
     1                  "element type = ", a9,/
     2                  "num_elem_in_block = ", i2,/
     3                  "num_nodes_per_elem = ", i2,/
     4                  "num_attr = ", i2)')
     5                  ids(i), eltype(i), num_elem_in_block(i),
     6                  num_nodes_per_elem(i), num_attr(i)

40    continue

c read element connectivity

      do 60 i = 1, num_elem_blk

        call exgelc (exoid, ids(i), connect, ierr)
        write (iout, '(/"after exgelc, error = ", i3)' ) ierr

        if (eltype(i) .eq. 'nsided' .or. eltype(i) .eq. 'NSIDED') then
          call exgecpp(exoid, EXEBLK, ids(i), nnpe, ierr)
          write (iout, '(/"after exgecpp, error = ", i3)' ) ierr

          write (iout, '("connect array for elem block ", i2)') ids(i)

          ioff = 0
          do j = 1, num_elem_in_block(i)
            write (iout, 100) j, nnpe(j), (connect(ioff+k),k=1,nnpe(j))
            ioff = ioff + nnpe(j)
          end do

        end if
60    continue

C ... Test coordinate frames
      call getfrm(exoid)

      call exclos (exoid, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr
 100  format(' Element ',I3,', Nodes/Element = ',I3,' -- ',20I3)
      stop
      end

      subroutine getfrm(exoid)
      implicit none
      include 'exodusII.inc'

      integer iout

      integer exoid, ierr, i, j
      integer numfrm;   ! Assumed to be 3 for remaining dimensions
      integer cfids(3), tags(3)
      character*32 ctag
      real    coord(27)

      data iout /6/

      numfrm = exinqi(exoid, EXNCF)
      write (iout,
     1  '(/"There are ",i2," coordinate frames")')
     2  numfrm

      numfrm = 3
      call exgfrm(exoid, numfrm, cfids, coord, tags, ierr);
      write (6,'("after exgfrm, error = ", i4)') ierr

! NOTE: These values may not be sensical; just used for testing.
      do i=0,2
        if (tags(i+1) .eq. EXCFREC) then
          ctag = 'RECTANGULAR'
        else if (tags(i+1) .eq. EXCFCYL) then
          ctag = 'CYLINDRICAL'
        else if (tags(i+1) .eq. EXCFSPH) then
          ctag = 'SPHERICAL'
        end if

        write (iout, 100) cfids(i+1), ctag, (COORD(9*i+j),j=1,9)
      end do

 100  format(i5, 3x, A, 9F10.2)
      return
      end
