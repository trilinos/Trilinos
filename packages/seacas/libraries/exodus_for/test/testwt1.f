C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      program testwt1

c This is a test program for the Fortran binding of the EXODUS II
c database write routines.

c       history -
c       Original L.A. Schoof
c       02/25/93 V.R. Yarberry - Added error checks for file creation.
c       03/04/93 V.R. Yarberry - Fixed bug in expvtt test, ebids was not passed
c       08/31/93 VRY - updated to match API version 2.00

      include 'exodusII.inc'

      integer iin, iout
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer num_elem_in_block(10), num_nodes_per_elem(10),numattr(10)
      integer num_node_sets, num_side_sets
      integer i, j, k, m, elem_map(10), node_map(100), connect(10)
      integer node_list(100), elem_list(100), side_list(100)
      integer ebids(10),ids(10), num_nodes_per_set(10)
      integer num_elem_per_set(10), num_df_per_set(10)
      integer df_ind(10), node_ind(10), elem_ind(10)
      integer num_qa_rec, num_info
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer truth_tab(3,7)
      integer whole_time_step, num_time_steps
      integer cpu_word_size, io_word_size
      integer prop_array(2)

      real glob_var_vals(100), nodal_var_vals(100)
      real time_value, elem_var_vals(100)
      real x(100), y(100), z(100)
      real attrib(100), dist_fact(100)

      character*(MXSTLN) coord_names(3)
      character*(MXSTLN) cname
      character*(MXSTLN) var_names(3)
      character*(MXSTLN) qa_record(4,2)
      character*(MXLNLN) inform(3)
      character*(MXSTLN) prop_names(2)
      character*(MXSTLN) attrib_names(3)
      character*(MXSTLN) blk_names(7)

      data iin /5/, iout /6/

      cpu_word_size = 0
      io_word_size = 0

c  create EXODUS II files

      exoid = excre ("test.exo",
     1               EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo, id: ", i4)') exoid
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout,'("after excre, error = ", i4)') ierr

c  initialize file with parameters

      num_dim = 3
      num_nodes = 28
      num_elem = 8
      num_elem_blk = 7
      num_node_sets = 2
      num_side_sets = 5
c     Uncomment the following line to test NULL side sets
c     num_side_sets = 6

      call expini (exoid, "This is testwt1", num_dim, num_nodes,
     1             num_elem, num_elem_blk, num_node_sets,
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i4)' ) ierr

c  write nodal coordinates values and names to database

c  Quad #1
      x(1) = 0.0
      x(2) = 1.0
      x(3) = 1.0
      x(4) = 0.0

      y(1) = 0.0
      y(2) = 0.0
      y(3) = 1.0
      y(4) = 1.0

      z(1) = 0.0
      z(2) = 0.0
      z(3) = 0.0
      z(4) = 0.0

c  Quad #2
      x(5) = 1.0
      x(6) = 2.0
      x(7) = 2.0
      x(8) = 1.0

      y(5) = 0.0
      y(6) = 0.0
      y(7) = 1.0
      y(8) = 1.0

      z(5) = 0.0
      z(6) = 0.0
      z(7) = 0.0
      z(8) = 0.0

c  Hex #1
      x(9)  =  0.0
      x(10) = 10.0
      x(11) = 10.0
      x(12) =  1.0
      x(13) =  1.0
      x(14) = 10.0
      x(15) = 10.0
      x(16) =  1.0

      y(9)  =  0.0
      y(10) =  0.0
      y(11) =  0.0
      y(12) =  0.0
      y(13) = 10.0
      y(14) = 10.0
      y(15) = 10.0
      y(16) = 10.0

      z(9)  =  0.0
      z(10) =  0.0
      z(11) =-10.0
      z(12) =-10.0
      z(13) =  0.0
      z(14) =  0.0
      z(15) =-10.0
      z(16) =-10.0

c  Tetra #1
      x(17) =  0.0
      x(18) =  1.0
      x(19) = 10.0
      x(20) =  7.0

      y(17) =  0.0
      y(18) =  0.0
      y(19) =  0.0
      y(20) =  5.0

      z(17) =  0.0
      z(18) =  5.0
      z(19) =  2.0
      z(20) =  3.0

c Circle #1
      x(21) = 100.0
      y(21) = 100.0
      z(21) =   0.0

c  Sphere #1
      x(22) = 50.0
      y(22) = 50.0
      z(22) = 20.0

c  Wedge #1
      x(23) =  3.0
      x(24) =  6.0
      x(25) =  0.0
      x(26) =  3.0
      x(27) =  6.0
      x(28) =  0.0

      y(23) =  0.0
      y(24) =  0.0
      y(25) =  0.0
      y(26) =  2.0
      y(27) =  2.0
      y(28) =  2.0

      z(23) =  6.0
      z(24) =  0.0
      z(25) =  0.0
      z(26) =  6.0
      z(27) =  2.0
      z(28) =  0.0
      call expcor (exoid, x, y, z, ierr)
      write (iout, '("after expcor, error = ", i4)' ) ierr

      coord_names(1) = "xcoor"
      coord_names(2) = "ycoor"
      coord_names(3) = "zcoor"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i4)' ) ierr

c write node and element map parameters

      n_node_maps = 1
      n_elem_maps = 2

      call expmp (exoid, n_node_maps, n_elem_maps, ierr)
      write (iout, '("after expmp, error = ", i4)' ) ierr

c write element map properties

      prop_names(1) = "ORDER"
      prop_names(2) = "NUMBER"
      call exppn(exoid,EXEMAP,2,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr

c write element order map

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      id = 111
      call expem (exoid, id, elem_map, ierr)
      write (iout, '("after expem, error = ", i4)' ) ierr

      call expp(exoid, EXEMAP, id, "ORDER", 1, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

c write element numbering map

      id = 222
C write map an element at a time...
      do 11 i = 1, num_elem
        elem_map(i) = i*2
        call exppem (exoid, id, i, 1, elem_map(i), ierr)
        write (iout, '("after exppem, error = ", i4)' ) ierr
11    continue

      call expp(exoid, EXEMAP, id, "NUMBER", 1, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

c write node map properties

      prop_names(1) = "NUMBER"
      call exppn(exoid,EXNMAP,1,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr

c write node numbering map

      do 13 i = 1, num_nodes
        node_map(i) = i*3
 13   continue

      id = 333
      call expnm (exoid, id, node_map, ierr)
      write (iout, '("after expnm, error = ", i4)' ) ierr

      call expp(exoid, EXNMAP, id, "NUMBER", 1, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

c write element block parameters

      num_elem_in_block(1) = 1
      num_elem_in_block(2) = 2
      num_elem_in_block(3) = 1
      num_elem_in_block(4) = 1
      num_elem_in_block(5) = 1
      num_elem_in_block(6) = 1
      num_elem_in_block(7) = 1

      num_nodes_per_elem(1) = 4
      num_nodes_per_elem(2) = 4
      num_nodes_per_elem(3) = 8
      num_nodes_per_elem(4) = 4
      num_nodes_per_elem(5) = 1
      num_nodes_per_elem(6) = 1
      num_nodes_per_elem(7) = 6

      ebids(1) = 10
      ebids(2) = 11
      ebids(3) = 12
      ebids(4) = 13
      ebids(5) = 14
      ebids(6) = 15
      ebids(7) = 16

      numattr(1) = 3
      numattr(2) = 3
      numattr(3) = 3
      numattr(4) = 3
      numattr(5) = 3
      numattr(6) = 3
      numattr(7) = 3

      blk_names(1) = "e_block_i"
      blk_names(2) = "e_block_ii"
      blk_names(3) = "e_block_iii"
      blk_names(4) = "e_block_iv"
      blk_names(5) = "e_block_v"
      blk_names(6) = "e_block_vi"
      blk_names(7) = "e_block_vii"

      cname = "quad"
      call expelb (exoid,ebids(1),cname,num_elem_in_block(1),
     1          num_nodes_per_elem(1),numattr(1),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(1),blk_names(1), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      call expelb (exoid,ebids(2),cname,num_elem_in_block(2),
     1          num_nodes_per_elem(2),numattr(2),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(2),blk_names(2), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      cname = "hex"
      call expelb (exoid,ebids(3),cname,num_elem_in_block(3),
     1          num_nodes_per_elem(3),numattr(3),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(3),blk_names(3), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      cname = "tetra"
      call expelb (exoid,ebids(4),cname,num_elem_in_block(4),
     1          num_nodes_per_elem(4),numattr(4),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(4),blk_names(4), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      cname = "circle"
      call expelb (exoid,ebids(5),cname,num_elem_in_block(5),
     1          num_nodes_per_elem(5),numattr(5),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(5),blk_names(5), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      cname = "sphere"
      call expelb (exoid,ebids(6),cname,num_elem_in_block(6),
     1          num_nodes_per_elem(6),numattr(6),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(6),blk_names(6), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

      cname = "wedge"
      call expelb (exoid,ebids(7),cname,num_elem_in_block(7),
     1          num_nodes_per_elem(7),numattr(7),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expnam (exoid,EXEBLK,ebids(7),blk_names(7), ierr)
      write (iout, '("after expnam, error = ", i4)' ) ierr

c  write element block properties

      prop_names(1) = "MATL"
      prop_names(2) = "DENSITY"
      call exppn(exoid,EXEBLK,2,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr

      call expp(exoid, EXEBLK, ebids(1), "MATL", 10, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(2), "MATL", 20, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(3), "MATL", 30, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(4), "MATL", 40, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(5), "MATL", 50, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(6), "MATL", 60, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(7), "MATL", 70, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

c write element connectivity

      connect(1) = 1
      connect(2) = 2
      connect(3) = 3
      connect(4) = 4

      call expelc (exoid, ebids(1), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 1
      connect(2) = 2
      connect(3) = 3
      connect(4) = 4
      connect(5) = 5
      connect(6) = 6
      connect(7) = 7
      connect(8) = 8

      call expelc (exoid, ebids(2), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) =  9
      connect(2) = 10
      connect(3) = 11
      connect(4) = 12
      connect(5) = 13
      connect(6) = 14
      connect(7) = 15
      connect(8) = 16

      call expelc (exoid, ebids(3), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 17
      connect(2) = 18
      connect(3) = 19
      connect(4) = 20

      call expelc (exoid, ebids(4), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 21

      call expelc (exoid, ebids(5), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 22

      call expelc (exoid, ebids(6), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 23
      connect(2) = 24
      connect(3) = 25
      connect(4) = 26
      connect(5) = 27
      connect(6) = 28

      call expelc (exoid, ebids(7), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

c write element block attributes

      attrib(1) = 1.0  !  block 1
      attrib(2) = 2.0
      attrib(3) = 3.0
      attrib(4) = 1.11 !  block 2, element 1
      attrib(5) = 2.11
      attrib(6) = 3.11
      attrib(7) = 1.12 !  block 2, element 2
      attrib(8) = 2.12
      attrib(9) = 3.12
      attrib(10) = 1.2  !  block 3
      attrib(11) = 2.2
      attrib(12) = 3.2
      attrib(13) = 1.3 !  block 4
      attrib(14) = 2.3
      attrib(15) = 3.3
      attrib(16) = 1.4 !  block 5
      attrib(17) = 2.4
      attrib(18) = 3.4
      attrib(19) = 1.5 !  block 6
      attrib(20) = 2.5
      attrib(21) = 3.5
      attrib(22) = 1.6 !  block 7
      attrib(23) = 2.6
      attrib(24) = 3.6

      call expeat (exoid, ebids(1), attrib(1), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(2), attrib(4), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(3), attrib(10), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(4), attrib(13), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(5), attrib(16), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(6), attrib(19), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(7), attrib(22), ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib_names(1) = "attribute_1"
      attrib_names(2) = "attribute_2"
      attrib_names(3) = "attribute_3"
      do i=1, num_elem_blk
        call expean (exoid, ebids(i), numattr(i), attrib_names, ierr)
        write (iout, '("after expean, error = ", i4)' ) ierr
      end do

c write individual node sets

      node_list(1) = 100
      node_list(2) = 101
      node_list(3) = 102
      node_list(4) = 103
      node_list(5) = 104

      dist_fact(1) = 1.0
      dist_fact(2) = 2.0
      dist_fact(3) = 3.0
      dist_fact(4) = 4.0
      dist_fact(5) = 5.0

c     call expnp (exoid, 20, 5, 5, ierr)
c     write (iout, '("after expnp, error = ", i4)' ) ierr
c     call expns (exoid, 20, node_list, ierr)
c     write (iout, '("after expns, error = ", i4)' ) ierr
c     call expnsd (exoid, 20, dist_fact, ierr)
c     write (iout, '("after expnsd, error = ", i4)' ) ierr

      node_list(1) = 200
      node_list(2) = 201
      node_list(3) = 202

      dist_fact(1) = 1.1
      dist_fact(2) = 2.1
      dist_fact(3) = 3.1

c     call expnp (exoid, 21, 3, 3, ierr)
c     write (iout, '("after expnp, error = ", i4)' ) ierr
c     call expns (exoid, 21, node_list, ierr)
c     write (iout, '("after expns, error = ", i4)' ) ierr
c     call expnsd (exoid, 21, dist_fact, ierr)
c     write (iout, '("after expnsd, error = ", i4)' ) ierr

c write concatenated node sets; this produces the same information as
c the above code which writes individual node sets

      ids(1) = 20
      ids(2) = 21

      num_nodes_per_set(1) = 5
      num_nodes_per_set(2) = 3

      num_df_per_set(1) = 5
      num_df_per_set(2) = 3

      node_ind(1) = 1
      node_ind(2) = 6

      df_ind(1) = 1
      df_ind(2) = 6

      node_list(1) = 100
      node_list(2) = 101
      node_list(3) = 102
      node_list(4) = 103
      node_list(5) = 104
      node_list(6) = 200
      node_list(7) = 201
      node_list(8) = 202

      dist_fact(1) = 1.0
      dist_fact(2) = 2.0
      dist_fact(3) = 3.0
      dist_fact(4) = 4.0
      dist_fact(5) = 5.0
      dist_fact(6) = 1.1
      dist_fact(7) = 2.1
      dist_fact(8) = 3.1

      call expcns (exoid, ids, num_nodes_per_set, num_df_per_set,
     1        node_ind, df_ind, node_list, dist_fact, ierr)
      write (iout, '("after expcns, error = ", i4)' ) ierr

c     write node set properties

      prop_names(1) = "FACE"
      call expp(exoid, EXNSET, 20, prop_names(1), 4, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      call expp(exoid, EXNSET, 21, prop_names(1), 5, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      prop_array(1) = 1000
      prop_array(2) = 2000

      prop_names(1) = "VELOCITY"
      call exppa(exoid, EXNSET, prop_names(1), prop_array, ierr)
      write (iout, '("after exppa, error = ", i4)' ) ierr

c write individual side sets

      elem_list(1) = 11
      elem_list(2) = 12

      side_list(1) = 1
      side_list(2) = 2

      dist_fact(1) = 30.0
      dist_fact(2) = 30.1
      dist_fact(3) = 30.2
      dist_fact(4) = 30.3

c     call expsp (exoid, 30, 2, 4, ierr)
c     write (iout, '("after expsp, error = ", i4)' ) ierr

c     call expss (exoid, 30, elem_list, side_list, ierr)
c     write (iout, '("after expss, error = ", i4)' ) ierr

c     call expssd (exoid, 30, dist_fact, ierr)
c     write (iout, '("after expssd, error = ", i4)' ) ierr

      elem_list(1) = 13
      elem_list(2) = 14

      side_list(1) = 3
      side_list(2) = 4

      dist_fact(1) = 31.0
      dist_fact(2) = 31.1
      dist_fact(3) = 31.2
      dist_fact(4) = 31.3

c     call expsp (exoid, 31, 2, 4, ierr)
c     write (iout, '("after expsp, error = ", i4)' ) ierr

c     call expss (exoid, 31, elem_list, side_list, ierr)
c     write (iout, '("after expss, error = ", i4)' ) ierr

c     call expssd (exoid, 31, dist_fact, ierr)
c     write (iout, '("after expssd, error = ", i4)' ) ierr

c write concatenated side sets; this produces the same information as
c the above code which writes individual side sets

      ids(1) = 30
      ids(2) = 31
      ids(3) = 32
      ids(4) = 33
      ids(5) = 34
      ids(6) = 35

c     side set #1 - quad
      node_list(1) = 8
      node_list(2) = 5
      node_list(3) = 6
      node_list(4) = 7

c     side set #2 - quad/hex, spanning 2 element types
      node_list(5) = 2
      node_list(6) = 3
      node_list(7) = 7
      node_list(8) = 8

c     side set #3 - hex
      node_list(9)  =  9
      node_list(10) = 12
      node_list(11) = 11
      node_list(12) = 10

      node_list(13) = 11
      node_list(14) = 12
      node_list(15) = 16
      node_list(16) = 15

      node_list(17) = 16
      node_list(18) = 15
      node_list(19) = 11
      node_list(20) = 12

      node_list(21) = 10
      node_list(22) = 11
      node_list(23) = 15
      node_list(24) = 14

      node_list(25) = 13
      node_list(26) = 16
      node_list(27) = 12
      node_list(28) =  9

      node_list(29) = 14
      node_list(30) = 13
      node_list(31) =  9
      node_list(32) = 10

      node_list(33) = 16
      node_list(34) = 13
      node_list(35) = 14
      node_list(36) = 15

c     side set #4 - Tetra
      node_list(37) = 17
      node_list(38) = 18
      node_list(39) = 20

      node_list(40) = 18
      node_list(41) = 19
      node_list(42) = 20

      node_list(43) = 20
      node_list(44) = 19
      node_list(45) = 17

      node_list(46) = 19
      node_list(47) = 18
      node_list(48) = 17

c     side set #5 - Circle/Sphere
      node_list(49) = 21
      node_list(50) = 22

c     side set #6 - Wedges
      node_list(51) = 27
      node_list(52) = 26
      node_list(53) = 23
      node_list(54) = 24

      node_list(55) = 28
      node_list(56) = 27
      node_list(57) = 24
      node_list(58) = 25

      node_list(59) = 28
      node_list(60) = 25
      node_list(61) = 23
      node_list(62) = 26

      node_list(63) = 25
      node_list(64) = 24
      node_list(65) = 23

      node_list(66) = 26
      node_list(67) = 27
      node_list(68) = 28

      num_elem_per_set(1) = 2
      num_elem_per_set(2) = 2
      num_elem_per_set(3) = 7
      num_elem_per_set(4) = 4
      num_elem_per_set(5) = 2
      num_elem_per_set(6) = 5
c     Uncomment following line to test NULL side sets
c     num_elem_per_set(6) = 0

      num_nodes_per_set(1) = 4
      num_nodes_per_set(2) = 4
      num_nodes_per_set(3) = 28
      num_nodes_per_set(4) = 12
      num_nodes_per_set(5) =  2
      num_nodes_per_set(6) = 18

      elem_ind(1) = 1
      elem_ind(2) = 3
      elem_ind(3) = 5
      elem_ind(4) = 12
      elem_ind(5) = 16
      elem_ind(6) = 18

      node_ind(1) = 1
      node_ind(2) = 5
      node_ind(3) = 9
      node_ind(4) = 37
      node_ind(5) = 48
      node_ind(6) = 50

      elem_list(1) = 3
      elem_list(2) = 3
      elem_list(3) = 1
      elem_list(4) = 3
      elem_list(5) = 4
      elem_list(6) = 4
      elem_list(7) = 4
      elem_list(8) = 4
      elem_list(9) = 4
      elem_list(10) = 4
      elem_list(11) = 4
      elem_list(12) = 5
      elem_list(13) = 5
      elem_list(14) = 5
      elem_list(15) = 5
      elem_list(16) = 6
      elem_list(17) = 7
      elem_list(18) = 8
      elem_list(19) = 8
      elem_list(20) = 8
      elem_list(21) = 8
      elem_list(22) = 8

c     side_list(1) = 1
c     side_list(2) = 2
c     side_list(3) = 3
c     side_list(4) = 4

      call excn2s(exoid, num_elem_per_set, num_nodes_per_set, elem_ind,
     1          node_ind, elem_list, node_list, side_list, ierr)
      write (iout, '("after excn2s, error = ", i4)' ) ierr

      num_df_per_set(1) = 4
      num_df_per_set(2) = 4
      num_df_per_set(3) = 0
      num_df_per_set(4) = 0
      num_df_per_set(5) = 0
      num_df_per_set(6) = 0

      df_ind(1) = 1
      df_ind(2) = 5
      df_ind(3) = 9
      df_ind(4) = 9
      df_ind(5) = 9
      df_ind(6) = 9

      dist_fact(1) = 30.0
      dist_fact(2) = 30.1
      dist_fact(3) = 30.2
      dist_fact(4) = 30.3
      dist_fact(5) = 31.0
      dist_fact(6) = 31.1
      dist_fact(7) = 31.2
      dist_fact(8) = 31.3

      call expcss (exoid, ids, num_elem_per_set, num_df_per_set,
     1             elem_ind, df_ind, elem_list, side_list, dist_fact,
     2             ierr)
      write (iout, '("after expcss, error = ", i4)' ) ierr

      prop_names(1) = "COLOR"
      call expp(exoid, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      call expp(exoid, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

c write QA records

      num_qa_rec = 2

      qa_record(1,1) = "TESTWT1 fortran version"
      qa_record(2,1) = "testwt1"
      qa_record(3,1) = "03/16/94"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "07/07/93"
      qa_record(4,2) = "16:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i4)' ) ierr

c write information records

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i4)' ) ierr

c write results variables parameters and names

      num_glo_vars = 1

      var_names(1) = "glo vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvnm (exoid, "g", 1, var_names(1), ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_nod_vars = 2

      var_names(1) = "nod_var0"
      var_names(2) = "nod_var1"

      call expvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_ele_vars = 3

      var_names(1) = "ele_var0"
      var_names(2) = "ele_var1"
      var_names(3) = "ele_var2"

      call expvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

c write element variable truth table

      k = 0

      do 30 i = 1,num_elem_blk
         do 20 j = 1,num_ele_vars
            truth_tab(j,i) = 1
20       continue
30    continue

      truth_tab(1,3) = 0

c     call expvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ierr)
c     write (iout, '("after expvtt, error = ", i4)' ) ierr

c for each time step, write the analysis results;
c the code below fills the arrays glob_var_vals,
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays

      whole_time_step = 1
      num_time_steps = 10

      do 110 i = 1, num_time_steps
        time_value = real(i)/100.

c write time value

        call exptim (exoid, whole_time_step, time_value, ierr)
        write (iout, '("after exptim, error = ", i4)' ) ierr

c write global variables

        do 50 j = 1, num_glo_vars
          glob_var_vals(j) = real(j+1) * time_value
50      continue

        call expgv (exoid, whole_time_step, num_glo_vars,
     1              glob_var_vals, ierr)
        write (iout, '("after expgv, error = ", i4)' ) ierr

c write nodal variables

        do 70 k = 1, num_nod_vars
          do 60 j = 1, num_nodes

            nodal_var_vals(j) = real(k) + (real(j) * time_value)

60        continue

          call expnv (exoid, whole_time_step, k, num_nodes,
     1                nodal_var_vals, ierr)
          write (iout, '("after expnv, error = ", i4)' ) ierr

70      continue

c write element variables

        do 100 k = 1, num_ele_vars
          do 90 j = 1, num_elem_blk
            do 80 m = 1, num_elem_in_block(j)

              elem_var_vals(m) = real(k+1) + real(j+1) +
     1                          (real(m)*time_value)
c             write(iout,*)'elem_var_val(',m,'): ',elem_var_vals(m)

80          continue

            if (k .eq. 1 .and. j .eq. 3) then
                continue        ! skip element block 3, variable 1
            else
              call expev (exoid, whole_time_step, k, ebids(j),
     1                  num_elem_in_block(j), elem_var_vals, ierr)
              write (iout, '("after expev, error = ", i4)' ) ierr
            endif

90        continue
100     continue

        whole_time_step = whole_time_step + 1

c update the data file; this should be done at the end of every time
c step to ensure that no data is lost if the analysis dies

        call exupda (exoid, ierr)
        write (iout, '("after exupda, error = ", i4)' ) ierr

110   continue

c close the EXODUS files

      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      stop
      end
