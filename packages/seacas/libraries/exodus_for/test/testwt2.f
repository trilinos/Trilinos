C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      program testwt2

c This is a test program for the Fortran binding of the EXODUS II
c database write routines. It tests multiple simultaneous output files.

c     09/07/93  V.R. Yarberry - Revised for 2.00 API

      include 'exodusII.inc'

      integer iin, iout
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer exoid2, num_dim2, num_nodes2, num_elem2, num_elem_blk2
      integer num_elem_in_block(10), num_node_sets
      integer num_elem_in_block2(10), num_node_sets2
      integer num_side_sets, num_nodes_per_elem(10), numattr(10)
      integer num_side_sets2, num_nodes_per_elem2(10), numattr2(10)
      integer i, j, k, m, elem_map(5), connect(10)
      integer elem_map2(5), connect2(10)
      integer node_list(100), elem_list(100), side_list(100)
      integer node_list2(100), elem_list2(100), side_list2(100)
      integer ebids(10),ids(10),num_nodes_per_set(10)
      integer num_elem_per_set(10), num_df_per_set(10)
      integer ebids2(10)
      integer df_ind(10), node_ind(10), elem_ind(10)
      integer num_qa_rec, num_info
      integer num_qa_rec2,num_info2
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer num_glo_vars2, num_nod_vars2, num_ele_vars2
      integer truth_tab(3,5)
      integer whole_time_step, num_time_steps
      integer cpu_word_size, io_word_size
      integer prop_array(2)

      real glob_var_vals(100), nodal_var_vals(100)
      real time_value, elem_var_vals(100)
      real time_value2
      real x(100), y(100), z(100)
      real x2(100), y2(100), z2(100)
      real attrib(100), dist_fact(100)
      real attrib2(100), dist_fact2(100)

      character*(MXLNLN) title
      character*(MXLNLN) title2
      character*(MXSTLN) coord_names(3)
      character*(MXSTLN) coord_names2(3)
      character*(MXSTLN) cname
      character*(MXSTLN) cname2
      character*(MXSTLN) var_names(3)
      character*(MXSTLN) var_names2(3)
      character*(MXSTLN) qa_record(4,2)
      character*(MXSTLN) qa_record2(4,2)
      character*(MXLNLN) inform(3)
      character*(MXLNLN) inform2(3)
      character*(MXSTLN) prop_names(2)

      data iin /5/, iout /6/

c  create EXODUS II files

      cpu_word_size = 0
      io_word_size = 4

c     first create a "regular" file that contains everything except
c     history variable info

      exoid = excre ("test.exo",
     1               EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo,id: ",i4,", err=",i3)')
     1           exoid, ierr
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout, '("after excre, error = ", i4)' ) ierr

      exoid2= excre ("test2.exo",
     1               EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test2.exo,id: ",i4,", err=",i3)')
     1           exoid2, ierr
      write (iout, '("after excre (2), error = ", i4)' ) ierr

c  initialize file with parameters

      title = "This is test 2"
      num_dim = 3
      num_nodes = 26
      num_elem = 5
      num_elem_blk = 5
      num_node_sets = 2
      num_side_sets = 5

      call expini (exoid, title, num_dim, num_nodes,
     1             num_elem, num_elem_blk, num_node_sets,
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i4)' ) ierr

      title2 = "This is test 2"
      num_dim2 = 3
      num_nodes2 = 26
      num_elem2 = 5
      num_elem_blk2 = 5
      num_node_sets2 = 2
      num_side_sets2 = 5

      call expini (exoid2, title2, num_dim2, num_nodes2,
     1             num_elem2, num_elem_blk2, num_node_sets2,
     2             num_side_sets2, ierr)

      write (iout, '("after expini (2), error = ", i4)' ) ierr

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

c  Wedge #1
      x(21) =  3.0
      x(22) =  6.0
      x(23) =  0.0
      x(24) =  3.0
      x(25) =  6.0
      x(26) =  0.0

      y(21) =  0.0
      y(22) =  0.0
      y(23) =  0.0
      y(24) =  2.0
      y(25) =  2.0
      y(26) =  2.0

      z(21) =  6.0
      z(22) =  0.0
      z(23) =  0.0
      z(24) =  6.0
      z(25) =  2.0
      z(26) =  0.0

      call expcor (exoid, x, y, z, ierr)
      write (iout, '("after expcor, error = ", i4)' ) ierr

c  Quad #1
      x2(1) = 0.0
      x2(2) = 1.0
      x2(3) = 1.0
      x2(4) = 0.0

      y2(1) = 0.0
      y2(2) = 0.0
      y2(3) = 1.0
      y2(4) = 1.0

      z2(1) = 0.0
      z2(2) = 0.0
      z2(3) = 0.0
      z2(4) = 0.0

c  Quad #2
      x2(5) = 1.0
      x2(6) = 2.0
      x2(7) = 2.0
      x2(8) = 1.0

      y2(5) = 0.0
      y2(6) = 0.0
      y2(7) = 1.0
      y2(8) = 1.0

      z2(5) = 0.0
      z2(6) = 0.0
      z2(7) = 0.0
      z2(8) = 0.0

c  Hex #1
      x2(9)  =  0.0
      x2(10) = 10.0
      x2(11) = 10.0
      x2(12) =  1.0
      x2(13) =  1.0
      x2(14) = 10.0
      x2(15) = 10.0
      x2(16) =  1.0

      y2(9)  =  0.0
      y2(10) =  0.0
      y2(11) =  0.0
      y2(12) =  0.0
      y2(13) = 10.0
      y2(14) = 10.0
      y2(15) = 10.0
      y2(16) = 10.0

      z2(9)  =  0.0
      z2(10) =  0.0
      z2(11) =-10.0
      z2(12) =-10.0
      z2(13) =  0.0
      z2(14) =  0.0
      z2(15) =-10.0
      z2(16) =-10.0

c  Tetra #1
      x2(17) =  0.0
      x2(18) =  1.0
      x2(19) = 10.0
      x2(20) =  7.0

      y2(17) =  0.0
      y2(18) =  0.0
      y2(19) =  0.0
      y2(20) =  5.0

      z2(17) =  0.0
      z2(18) =  5.0
      z2(19) =  2.0
      z2(20) =  3.0

c  Wedge #1
      x2(21) =  3.0
      x2(22) =  6.0
      x2(23) =  0.0
      x2(24) =  3.0
      x2(25) =  6.0
      x2(26) =  0.0

      y2(21) =  0.0
      y2(22) =  0.0
      y2(23) =  0.0
      y2(24) =  2.0
      y2(25) =  2.0
      y2(26) =  2.0

      z2(21) =  6.0
      z2(22) =  0.0
      z2(23) =  0.0
      z2(24) =  6.0
      z2(25) =  2.0
      z2(26) =  0.0

      call expcor (exoid2, x2, y2, z2, ierr)
      write (iout, '("after expcor (2), error = ", i4)' ) ierr

      coord_names(1) = "xcoor"
      coord_names(2) = "ycoor"
      coord_names(3) = "zcoor"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i4)' ) ierr

      coord_names2(1) = "xcoor"
      coord_names2(2) = "ycoor"
      coord_names2(3) = "zcoor"

      call expcon (exoid2, coord_names2, ierr)
      write (iout, '("after expcon (2), error = ", i4)' ) ierr

c write element order map

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      call expmap (exoid, elem_map, ierr)
      write (iout, '("after expmap, error = ", i4)' ) ierr

      do 12 i = 1, num_elem2
         elem_map2(i) = i
12    continue

      call expmap (exoid2, elem_map2, ierr)
      write (iout, '("after expmap (2), error = ", i4)' ) ierr

c write element block parameters

      num_elem_in_block(1) = 1
      num_elem_in_block(2) = 1
      num_elem_in_block(3) = 1
      num_elem_in_block(4) = 1
      num_elem_in_block(5) = 1

      num_nodes_per_elem(1) = 4
      num_nodes_per_elem(2) = 4
      num_nodes_per_elem(3) = 8
      num_nodes_per_elem(4) = 4
      num_nodes_per_elem(5) = 6

      ebids(1) = 10
      ebids(2) = 11
      ebids(3) = 12
      ebids(4) = 13
      ebids(5) = 14

      numattr(1) = 1
      numattr(2) = 1
      numattr(3) = 1
      numattr(4) = 1
      numattr(5) = 1

      cname = "quad"

      call expelb (exoid,ebids(1),cname,num_elem_in_block(1),
     1          num_nodes_per_elem(1),numattr(1),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expelb (exoid,ebids(2),cname,num_elem_in_block(2),
     1          num_nodes_per_elem(2),numattr(2),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      cname = "hex"
      call expelb (exoid,ebids(3),cname,num_elem_in_block(3),
     1          num_nodes_per_elem(3),numattr(3),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      cname = "tetra"
      call expelb (exoid,ebids(4),cname,num_elem_in_block(4),
     1          num_nodes_per_elem(4),numattr(4),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      cname = "wedge"
      call expelb (exoid,ebids(5),cname,num_elem_in_block(5),
     1          num_nodes_per_elem(5),numattr(5),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      num_elem_in_block2(1) = 1
      num_elem_in_block2(2) = 1
      num_elem_in_block2(3) = 1
      num_elem_in_block2(4) = 1
      num_elem_in_block2(5) = 1

      num_nodes_per_elem2(1) = 4
      num_nodes_per_elem2(2) = 4
      num_nodes_per_elem2(3) = 8
      num_nodes_per_elem2(4) = 4
      num_nodes_per_elem2(5) = 6

      ebids2(1) = 10
      ebids2(2) = 11
      ebids2(3) = 12
      ebids2(4) = 13
      ebids2(5) = 14

      numattr2(1) = 1
      numattr2(2) = 1
      numattr2(3) = 1
      numattr2(4) = 1
      numattr2(5) = 1

      cname2 = "quad"

      call expelb(exoid2,ebids2(1),cname2,num_elem_in_block2(1),
     1          num_nodes_per_elem2(1),numattr2(1),ierr)
      write (iout, '("after expelb (2), error = ", i4)' ) ierr

      call expelb(exoid2,ebids2(2),cname2,num_elem_in_block2(2),
     1          num_nodes_per_elem2(2),numattr2(2),ierr)
      write (iout, '("after expelb (2), error = ", i4)' ) ierr

      cname2 = "hex"
      call expelb(exoid2,ebids2(3),cname2,num_elem_in_block2(3),
     1          num_nodes_per_elem(3),numattr(3),ierr)
      write (iout, '("after expelb (2), error = ", i4)' ) ierr

      cname2 = "tetra"
      call expelb(exoid2,ebids2(4),cname2,num_elem_in_block2(4),
     1          num_nodes_per_elem2(4),numattr2(4),ierr)
      write (iout, '("after expelb (2), error = ", i4)' ) ierr

      cname2 = "wedge"
      call expelb(exoid2,ebids2(5),cname2,num_elem_in_block2(5),
     1          num_nodes_per_elem2(5),numattr2(5),ierr)
      write (iout, '("after expelb (2), error = ", i4)' ) ierr

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

      call exppn(exoid2,EXEBLK,2,prop_names,ierr)
      write (iout, '("after exppn (2), error = ", i4)' ) ierr

      call expp(exoid2, EXEBLK, ebids(1), "MATL", 100, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr
      call expp(exoid2, EXEBLK, ebids(2), "MATL", 200, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr
      call expp(exoid2, EXEBLK, ebids(3), "MATL", 300, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr
      call expp(exoid2, EXEBLK, ebids(4), "MATL", 400, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr
      call expp(exoid2, EXEBLK, ebids(5), "MATL", 500, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr

c write element connectivity

      connect(1) = 1
      connect(2) = 2
      connect(3) = 3
      connect(4) = 4

      call expelc (exoid, ebids(1), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect(1) = 5
      connect(2) = 6
      connect(3) = 7
      connect(4) = 8

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
      connect(2) = 22
      connect(3) = 23
      connect(4) = 24
      connect(5) = 25
      connect(6) = 26

      call expelc (exoid, ebids(5), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr

      connect2(1) = 1
      connect2(2) = 2
      connect2(3) = 3
      connect2(4) = 4

      call expelc (exoid2, ebids2(1), connect2, ierr)
      write (iout, '("after expelc (2), error = ", i4)' ) ierr

      connect2(1) = 5
      connect2(2) = 6
      connect2(3) = 7
      connect2(4) = 8

      call expelc (exoid2, ebids2(2), connect2, ierr)
      write (iout, '("after expelc (2), error = ", i4)' ) ierr

      connect2(1) =  9
      connect2(2) = 10
      connect2(3) = 11
      connect2(4) = 12
      connect2(5) = 13
      connect2(6) = 14
      connect2(7) = 15
      connect2(8) = 16

      call expelc (exoid2, ebids2(3), connect2, ierr)
      write (iout, '("after expelc (2), error = ", i4)' ) ierr

      connect2(1) = 17
      connect2(2) = 18
      connect2(3) = 19
      connect2(4) = 20

      call expelc (exoid2, ebids2(4), connect2, ierr)
      write (iout, '("after expelc (2), error = ", i4)' ) ierr

      connect2(1) = 21
      connect2(2) = 22
      connect2(3) = 23
      connect2(4) = 24
      connect2(5) = 25
      connect2(6) = 26

      call expelc (exoid2, ebids2(5), connect2, ierr)
      write (iout, '("after expelc (2), error = ", i4)' ) ierr

c write element block attributes

      attrib(1) = 3.14159
      call expeat (exoid, ebids(1), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib(1) = 6.14159
      call expeat (exoid, ebids(2), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(3), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(4), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      call expeat (exoid, ebids(5), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib2(1) = 3.
      call expeat (exoid2, ebids2(1), attrib2, ierr)
      write (iout, '("after expeat (2), error = ", i4)' ) ierr

      attrib2(1) = 6.
      call expeat (exoid2, ebids2(2), attrib2, ierr)
      write (iout, '("after expeat (2), error = ", i4)' ) ierr

      call expeat (exoid2, ebids2(3), attrib2, ierr)
      write (iout, '("after expeat (2), error = ", i4)' ) ierr

      call expeat (exoid2, ebids2(4), attrib2, ierr)
      write (iout, '("after expeat (2), error = ", i4)' ) ierr

      call expeat (exoid2, ebids(5), attrib2, ierr)
      write (iout, '("after expeat (2), error = ", i4)' ) ierr

c write individual node sets

      call expnp (exoid, 20, 5, 5, ierr)
      write (iout, '("after expnp, error = ", i4)' ) ierr

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

      call expns (exoid, 20, node_list, ierr)
      write (iout, '("after expns, error = ", i4)' ) ierr
      call expnsd (exoid, 20, dist_fact, ierr)
      write (iout, '("after expnsd, error = ", i4)' ) ierr

      call expnp (exoid, 21, 3, 3, ierr)
      write (iout, '("after expnp, error = ", i4)' ) ierr

      node_list(1) = 200
      node_list(2) = 201
      node_list(3) = 202

      dist_fact(1) = 1.1
      dist_fact(2) = 2.1
      dist_fact(3) = 3.1

      call expns (exoid, 21, node_list, ierr)
      write (iout, '("after expns, error = ", i4)' ) ierr
      call expnsd (exoid, 21, dist_fact, ierr)
      write (iout, '("after expnsd, error = ", i4)' ) ierr

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

C**** file 2

      call expnp (exoid2, 20, 5, 5, ierr)
      write (iout, '("after expnp (2), error = ", i4)' ) ierr

      node_list2(1) = 100
      node_list2(2) = 101
      node_list2(3) = 102
      node_list2(4) = 103
      node_list2(5) = 104

      dist_fact2(1) = 1.0
      dist_fact2(2) = 2.0
      dist_fact2(3) = 3.0
      dist_fact2(4) = 4.0
      dist_fact2(5) = 5.0

      call expns (exoid2, 20, node_list2, ierr)
      write (iout, '("after expns (2), error = ", i4)' ) ierr
      call expnsd (exoid2, 20, dist_fact2, ierr)
      write (iout, '("after expnsd (2), error = ", i4)' ) ierr

      call expnp (exoid2, 21, 3, 3, ierr)
      write (iout, '("after expnp (2), error = ", i4)' ) ierr

      node_list2(1) = 200
      node_list2(2) = 201
      node_list2(3) = 202

      dist_fact2(1) = 1.1
      dist_fact2(2) = 2.1
      dist_fact2(3) = 3.1

      call expns (exoid2, 21, node_list2, ierr)
      write (iout, '("after expns (2), error = ", i4)' ) ierr
      call expnsd (exoid2, 21, dist_fact2, ierr)
      write (iout, '("after expnsd (2), error = ", i4)' ) ierr

c write concatenated node sets; this produces the same information as
c the above code which writes individual node sets

      ids(1) = 20
      ids(2) = 21

      num_nodes_per_set(1) = 5
      num_nodes_per_set(2) = 3

      node_ind(1) = 1
      node_ind(2) = 6

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

c     call expcns (exoid, ids, num_nodes_per_set, node_ind, node_list,
c    1        dist_fact, ierr)
c     write (iout, '("after expcns, error = ", i4)' ) ierr

      prop_names(1) = "FACE"
      call expp(exoid2, EXNSET, 20, prop_names(1), 4, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr

      prop_names(1) = "FACE"
      call expp(exoid2, EXNSET, 21, prop_names(1), 5, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr

      prop_array(1) = 1000
      prop_array(2) = 2000

      prop_names(1) = "VELOCITY"
      call exppa(exoid2, EXNSET, prop_names(1), prop_array, ierr)
      write (iout, '("after exppa (2), error = ", i4)' ) ierr

c write individual side sets

c     side set #1 - quad

      elem_list(1) = 2
      elem_list(2) = 2

      side_list(1) = 4
      side_list(2) = 2

      dist_fact(1) = 30.0
      dist_fact(2) = 30.1
      dist_fact(3) = 30.2
      dist_fact(4) = 30.3

      call expsp (exoid, 30, 2, 4, ierr)
      write (iout, '("after expsp, error = ", i4)' ) ierr

      call expss (exoid, 30, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr

      call expssd (exoid, 30, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i4)' ) ierr

c     side set #2 - quad, spanning 2 elements

      elem_list(1) = 1
      elem_list(2) = 2

      side_list(1) = 2
      side_list(2) = 3

      dist_fact(1) = 31.0
      dist_fact(2) = 31.1
      dist_fact(3) = 31.2
      dist_fact(4) = 31.3

      call expsp (exoid, 31, 2, 4, ierr)
      write (iout, '("after expsp, error = ", i3)' ) ierr

      call expss (exoid, 31, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i3)' ) ierr

      call expssd (exoid, 31, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i3)' ) ierr

c     side set #3 - hex

      elem_list(1) = 3
      elem_list(2) = 3
      elem_list(3) = 3
      elem_list(4) = 3
      elem_list(5) = 3
      elem_list(6) = 3
      elem_list(7) = 3

      side_list(1) = 5
      side_list(2) = 3
      side_list(3) = 3
      side_list(4) = 2
      side_list(5) = 4
      side_list(6) = 1
      side_list(7) = 6

      call expsp (exoid, 32, 7, 0, ierr)
      write (iout, '("after expsp, error = ", i4)' ) ierr

      call expss (exoid, 32, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr

c     side set #4 - tetras

      elem_list(1) = 4
      elem_list(2) = 4
      elem_list(3) = 4
      elem_list(4) = 4

      side_list(1) = 1
      side_list(2) = 2
      side_list(3) = 3
      side_list(4) = 4

      call expsp (exoid, 33, 4, 0, ierr)
      write (iout, '("after expsp, error = ", i4)' ) ierr

      call expss (exoid, 33, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr

c     side set #5 - wedges

      elem_list(1) = 5
      elem_list(2) = 5
      elem_list(3) = 5
      elem_list(4) = 5
      elem_list(5) = 5

      side_list(1) = 1
      side_list(2) = 2
      side_list(3) = 3
      side_list(4) = 4
      side_list(5) = 5

      call expsp (exoid, 34, 5, 0, ierr)
      write (iout, '("after expsp, error = ", i4)' ) ierr

      call expss (exoid, 34, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr

c     side set #1  - quad

      elem_list2(1) = 2
      elem_list2(2) = 2

      side_list2(1) = 4
      side_list2(2) = 2

      dist_fact2(1) = 30.0
      dist_fact2(2) = 30.1
      dist_fact2(3) = 30.2
      dist_fact2(4) = 30.3

      call expsp (exoid2, 30, 2, 4, ierr)
      write (iout, '("after expsp (2), error = ", i4)' ) ierr

      call expss (exoid2, 30, elem_list2, side_list2, ierr)
      write (iout, '("after expss (2), error = ", i4)' ) ierr

      call expssd (exoid2, 30, dist_fact2, ierr)
      write (iout, '("after expssd (2), error = ", i4)' ) ierr

c     side set #2 - quad, spanning 2 elements

      elem_list2(1) = 1
      elem_list2(2) = 2

      side_list2(1) = 2
      side_list2(2) = 3

      dist_fact2(1) = 31.0
      dist_fact2(2) = 31.1
      dist_fact2(3) = 31.2
      dist_fact2(4) = 31.3

      call expsp (exoid2, 31, 2, 4, ierr)
      write (iout, '("after expsp (2), error = ", i3)' ) ierr

      call expss (exoid2, 31, elem_list2, side_list2, ierr)
      write (iout, '("after expss (2), error = ", i3)' ) ierr

      call expssd (exoid2, 31, dist_fact2, ierr)
      write (iout, '("after expssd (2), error = ", i3)' ) ierr

c     side set #3 - hex

      elem_list2(1) = 3
      elem_list2(2) = 3
      elem_list2(3) = 3
      elem_list2(4) = 3
      elem_list2(5) = 3
      elem_list2(6) = 3
      elem_list2(7) = 3

      side_list2(1) = 5
      side_list2(2) = 3
      side_list2(3) = 3
      side_list2(4) = 2
      side_list2(5) = 4
      side_list2(6) = 1
      side_list2(7) = 6

      call expsp (exoid2, 32, 7, 0, ierr)
      write (iout, '("after expsp (2), error = ", i4)' ) ierr

      call expss (exoid2, 32, elem_list2, side_list2, ierr)
      write (iout, '("after expss (2), error = ", i4)' ) ierr

c     side set #4 - tetras

      elem_list2(1) = 4
      elem_list2(2) = 4
      elem_list2(3) = 4
      elem_list2(4) = 4

      side_list2(1) = 1
      side_list2(2) = 2
      side_list2(3) = 3
      side_list2(4) = 4

      call expsp (exoid2, 33, 4, 0, ierr)
      write (iout, '("after expsp (2), error = ", i4)' ) ierr

      call expss (exoid2, 33, elem_list2, side_list2, ierr)
      write (iout, '("after expss (2), error = ", i4)' ) ierr

c     side set #5 - wedges

      elem_list2(1) = 5
      elem_list2(2) = 5
      elem_list2(3) = 5
      elem_list2(4) = 5
      elem_list2(5) = 5

      side_list2(1) = 1
      side_list2(2) = 2
      side_list2(3) = 3
      side_list2(4) = 4
      side_list2(5) = 5

      call expsp (exoid2, 34, 5, 0, ierr)
      write (iout, '("after expsp (2), error = ", i4)' ) ierr

      call expss (exoid2, 34, elem_list2, side_list2, ierr)
      write (iout, '("after expss (2), error = ", i4)' ) ierr

c write concatenated side sets; this produces the same information as
c the above code which writes individual side sets

      ids(1) = 30
      ids(2) = 31
      ids(3) = 32
      ids(4) = 33
      ids(5) = 34

c     side set #1
      node_list(1) = 8
      node_list(2) = 5
      node_list(3) = 6
      node_list(4) = 7

c     side set #2
      node_list(5) = 2
      node_list(6) = 3
      node_list(7) = 7
      node_list(8) = 8

c     side set #3
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

c     side set #4
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

c     side set #5
      node_list(49) = 25
      node_list(50) = 24
      node_list(51) = 21
      node_list(52) = 22

      node_list(53) = 26
      node_list(54) = 25
      node_list(55) = 22
      node_list(56) = 23

      node_list(57) = 26
      node_list(58) = 23
      node_list(59) = 21
      node_list(60) = 24

      node_list(61) = 23
      node_list(62) = 22
      node_list(63) = 21

      node_list(64) = 24
      node_list(65) = 25
      node_list(66) = 26

      num_elem_per_set(1) = 2
      num_elem_per_set(2) = 2
      num_elem_per_set(3) = 7
      num_elem_per_set(4) = 4
      num_elem_per_set(5) = 5

      num_nodes_per_set(1) = 4
      num_nodes_per_set(2) = 4
      num_nodes_per_set(3) = 28
      num_nodes_per_set(4) = 12
      num_nodes_per_set(5) = 20

      elem_ind(1) = 1
      elem_ind(2) = 3
      elem_ind(3) = 5
      elem_ind(4) = 12
      elem_ind(5) = 16

      node_ind(1) = 1
      node_ind(2) = 5
      node_ind(3) = 9
      node_ind(4) = 37
      node_ind(5) = 48

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
      elem_list(17) = 6
      elem_list(18) = 6
      elem_list(19) = 6
      elem_list(20) = 6

c     side_list(1) = 1
c     side_list(2) = 2
c     side_list(3) = 3
c     side_list(4) = 4

c     call excn2s(exoid, num_elem_per_set, num_nodes_per_set, elem_ind,
c    1          node_ind, elem_list, node_list, side_list, ierr)
c     write (iout, '("after excn2s, error = ", i4)' ) ierr

      num_df_per_set(1) = 4
      num_df_per_set(2) = 4
      num_df_per_set(3) = 0
      num_df_per_set(4) = 0
      num_df_per_set(5) = 0

      df_ind(1) = 1
      df_ind(2) = 5

      dist_fact(1) = 30.0
      dist_fact(2) = 30.1
      dist_fact(3) = 30.2
      dist_fact(4) = 30.3
      dist_fact(5) = 31.0
      dist_fact(6) = 31.1
      dist_fact(7) = 31.2
      dist_fact(8) = 31.3

c     call expcss (exoid, ids, num_elem_per_set, num_df_per_set,
c    1             elem_ind, df_ind, elem_list, side_list, dist_fact,
c    2             ierr)
c     write (iout, '("after expcss, error = ", i4)' ) ierr

c     call expcss (exoid2, ids, num_elem_per_set, num_df_per_set,
c    1             elem_ind, df_ind, elem_list, side_list, dist_fact,
c    2             ierr)
c     write (iout, '("after expcss (2), error = ", i4)' ) ierr

      prop_names(1) = "COLOR"
      call expp(exoid, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      call expp(exoid, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      prop_names(1) = "COLOR"
      call expp(exoid2, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr

      call expp(exoid2, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp (2), error = ", i4)' ) ierr

c write QA records

      num_qa_rec = 2

      qa_record(1,1) = "TESTWT2 fortran version"
      qa_record(2,1) = "testwt2"
      qa_record(3,1) = "07/07/93"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "07/07/93"
      qa_record(4,2) = "16:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i4)' ) ierr

      num_qa_rec2 = 2

      qa_record2(1,1) = "TESTWT2 fortran version"
      qa_record2(2,1) = "testwt2"
      qa_record2(3,1) = "07/07/93"
      qa_record2(4,1) = "15:41:33"
      qa_record2(1,2) = "FASTQ"
      qa_record2(2,2) = "fastq"
      qa_record2(3,2) = "07/07/93"
      qa_record2(4,2) = "16:41:33"

      call expqa (exoid2, num_qa_rec2, qa_record2, ierr)
      write (iout, '("after expqa (2), error = ", i4)' ) ierr

c write information records

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i4)' ) ierr

      num_info2 = 3

      inform2(1) = "This is the first information record."
      inform2(2) = "This is the second information record."
      inform2(3) = "This is the third information record."

      call expinf (exoid2, num_info2, inform2, ierr)
      write (iout, '("after expinf (2), error = ", i4)' ) ierr

c write results variables parameters and names

      num_glo_vars = 1

      var_names(1) = "glo_vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "g", num_glo_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_glo_vars2 = 1

      var_names2(1) = "glo_vars"

      call expvp (exoid2, "g", num_glo_vars2, ierr)
      write (iout, '("after expvp (2), error = ", i4)' ) ierr
      call expvan (exoid2, "g", num_glo_vars2, var_names2, ierr)
      write (iout, '("after expvan (2), error = ", i4)' ) ierr

      num_nod_vars = 2

      var_names(1) = "nod_var0"
      var_names(2) = "nod_var1"

      call expvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_nod_vars2 = 2

      var_names2(1) = "nod_var0"
      var_names2(2) = "nod_var1"

      call expvp (exoid2, "n", num_nod_vars2, ierr)
      write (iout, '("after expvp (2), error = ", i4)' ) ierr
      call expvan (exoid2, "n", num_nod_vars2, var_names2, ierr)
      write (iout, '("after expvan (2), error = ", i4)' ) ierr

      num_ele_vars = 3

      var_names(1) = "ele_var0"
      var_names(2) = "ele_var1"
      var_names(3) = "ele_var2"

      call expvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_ele_vars2 = 3

      var_names2(1) = "ele_var0"
      var_names2(2) = "ele_var1"
      var_names2(3) = "ele_var2"

      call expvp (exoid2, "e", num_ele_vars2, ierr)
      write (iout, '("after expvp (2), error = ", i4)' ) ierr
      call expvan (exoid2, "e", num_ele_vars2, var_names2, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

c write element variable truth table

      k = 0

      do 30 i = 1,num_elem_blk
         do 20 j = 1,num_ele_vars
            truth_tab(j,i) = 1
20       continue
30    continue

      call exgebi (exoid, ebids, ierr)
      write (iout, '("after exgebi, error = ", i4)' ) ierr
      call exgebi (exoid2, ebids2, ierr)
      write (iout, '("after exgebi (2), error = ", i4)' ) ierr
      call expvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ierr)
      write (iout, '("after expvtt, error = ", i4)' ) ierr
      call expvtt (exoid2, num_elem_blk, num_ele_vars, truth_tab, ierr)
      write (iout, '("after expvtt, error = ", i4)' ) ierr

c for each time step, write the analysis results;
c the code below fills the arrays glob_var_vals,
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays

      whole_time_step = 1
      num_time_steps = 10

      do 110 i = 1, num_time_steps
        time_value = real(i)/100
        time_value2 = real(i)/100

c write time value to regular file

        call exptim (exoid, whole_time_step, time_value, ierr)
        write (iout, '("after exptim, error = ", i4)' ) ierr

        call exptim (exoid2, whole_time_step, time_value2, ierr)
        write (iout, '("after exptim (2), error = ", i4)' ) ierr

c write global variables

        do 50 j = 1, num_glo_vars
          glob_var_vals(j) = real(j+1) * time_value
50      continue

        call expgv (exoid, whole_time_step, num_glo_vars,
     1              glob_var_vals, ierr)
        write (iout, '("after expgv, error = ", i4)' ) ierr

        call expgv (exoid2, whole_time_step, num_glo_vars,
     1              glob_var_vals, ierr)
        write (iout, '("after expgv (2), error = ", i4)' ) ierr

c write nodal variables

        do 70 k = 1, num_nod_vars
          do 60 j = 1, num_nodes

            nodal_var_vals(j) = real(k) + (real(j) * time_value)

60        continue

          call expnv (exoid, whole_time_step, k, num_nodes,
     1                nodal_var_vals, ierr)
          write (iout, '("after expnv, error = ", i4)' ) ierr

          call expnv (exoid2, whole_time_step, k, num_nodes,
     1                nodal_var_vals, ierr)
          write (iout, '("after expnv (2), error = ", i4)' ) ierr

70      continue

c write element variables

        do 100 k = 1, num_ele_vars
          do 90 j = 1, num_elem_blk
            do 80 m = 1, num_elem_in_block(j)

              elem_var_vals(m) = real(k+1) + real(j+1) +
     1                          (real(m)*time_value)

80          continue

            call expev (exoid, whole_time_step, k, ebids(j),
     1                  num_elem_in_block(j), elem_var_vals, ierr)
            write (iout, '("after expev, error = ", i4)' ) ierr
            call expev (exoid2, whole_time_step, k, ebids(j),
     1                  num_elem_in_block(j), elem_var_vals, ierr)
            write (iout, '("after expev (2), error = ", i4)' ) ierr

90        continue
100     continue

        whole_time_step = whole_time_step + 1

c update the data file; this should be done at the end of every time
c step to ensure that no data is lost if the analysis dies

        call exupda (exoid, ierr)
        write (iout, '("after exupda, error = ", i4)' ) ierr
        call exupda (exoid2, ierr)
        write (iout, '("after exupda (2), error = ", i4)' ) ierr

110   continue

c close the EXODUS files

      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      call exclos (exoid2, ierr)
      write (iout, '("after exclos (2), error = ", i4)' ) ierr

      stop
      end

