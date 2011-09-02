      program testwtm
c
c This is a test program for the Fortran binding of the EXODUS II
c database write routines. It tests multiple simultaneous output files.
c
c     09/07/93	V.R. Yarberry - Revised for 2.00 API

      include 'exodusII.inc'

      integer iin, iout
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer exoidm(10),num_dim2,num_nodes2,num_elem2,num_elem_blk2
      integer num_elem_in_block(2), num_node_sets
      integer num_elem_in_block2(2), num_node_sets2
      integer num_side_sets
      integer num_side_sets2
      integer nexofiles
      integer i, j, k, m, elem_map(2), connect(4) 
      integer elem_map2(2), connect2(4) 
      integer node_list(10), elem_list(10), side_list(10)
      integer node_list2(10), elem_list2(10), side_list2(10)
      integer ebids(2),ids(2), num_nodes_per_set(2), num_elem_per_set(2)
      integer ebids2(2)
      integer num_df_per_set(2)
      integer df_ind(2), node_ind(2), elem_ind(2), num_qa_rec, num_info
      integer num_qa_rec2,num_info2
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer num_glo_vars2, num_nod_vars2, num_ele_vars2
      integer truth_tab(3,2)
      integer whole_time_step, num_time_steps
      integer cpu_word_size, io_word_size
      integer prop_array(2)

      real glob_var_vals(10), nodal_var_vals(8) 
      real time_value, elem_var_vals(20)
      real time_value2
      real x(8), y(8), dummy(1)
      real x2(8), y2(8)
      real attrib(1), dist_fact(8)
      real attrib2(1), dist_fact2(8)

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
      character*(MXSTLN) exofname

      data iin /5/, iout /6/, nexofiles /5/

c
c  create EXODUS II files 
c
      cpu_word_size = 0
      io_word_size = 4
c
      exoid = excre ("test.exo",
     1               EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo,id: ",i4,", err=",i3)')
     1           exoid, ierr
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout, '("after excre, error = ", i4)' ) ierr

      do 1000 i=1,nexofiles
        write(exofname,'("test",i1,".exo")')i
        exoidm(i)= excre (exofname, 
     1               EXCLOB, cpu_word_size, io_word_size, ierr)
        write (iout,
     1    '("after excre for test",i1,".exo,id: ",i4,", err=",i3)')
     2	  i, exoidm(i), ierr
        write (iout, '("after excre (",i1,"), error = ", i4)' )
     1		i, ierr
1000  continue

c
c  initialize file with parameters
c

      title = "This is test m"
      num_dim = 2
      num_nodes = 8
      num_elem = 2
      num_elem_blk = 2
      num_node_sets = 2
      num_side_sets = 2

      call expini (exoid, title, num_dim, num_nodes, 
     1             num_elem, num_elem_blk, num_node_sets, 
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i4)' ) ierr
      
      title2 = "This is test m"
      num_dim2 = 2
      num_nodes2 = 8
      num_elem2 = 2
      num_elem_blk2 = 2
      num_node_sets2 = 2
      num_side_sets2 = 2

      do 1001 i=1,nexofiles
        call expini (exoidm(i), title2, num_dim2, num_nodes2, 
     1             num_elem2, num_elem_blk2, num_node_sets2, 
     2             num_side_sets2, ierr)

        write (iout, '("after expini (",i1,"), error = ", i4)' )
     1		i, ierr
1001  continue


c
c  write nodal coordinates values and names to database
c

      x(1) = 0.0 
      x(2) = 1.0 
      x(3) = 1.0 
      x(4) = 0.0 
      x(5) = 1.0 
      x(6) = 2.0 
      x(7) = 2.0 
      x(8) = 1.0
      y(1) = 0.0 
      y(2) = 0.0 
      y(3) = 1.0 
      y(4) = 1.0 
      y(5) = 0.0 
      y(6) = 0.0 
      y(7) = 1.0 
      y(8) = 1.0

      call expcor (exoid, x, y, dummy, ierr)
      write (iout, '("after expcor, error = ", i4)' ) ierr

      x2(1) = 0.0 
      x2(2) = 1.0 
      x2(3) = 1.0 
      x2(4) = 0.0 
      x2(5) = 1.0 
      x2(6) = 2.0 
      x2(7) = 2.0 
      x2(8) = 1.0
      y2(1) = 0.0 
      y2(2) = 0.0 
      y2(3) = 1.0 
      y2(4) = 1.0 
      y2(5) = 0.0 
      y2(6) = 0.0 
      y2(7) = 1.0 
      y2(8) = 1.0

      do 1002 i=1,nexofiles
        call expcor (exoidm(i), x2, y2, dummy, ierr)
        write (iout, '("after expcor (",i1,"), error = ", i4)')
     1		i, ierr
1002  continue

      coord_names(1) = "xcoor"
      coord_names(2) = "ycoor"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i4)' ) ierr

      coord_names2(1) = "xcoor"
      coord_names2(2) = "ycoor"

      do 1003 i=1,nexofiles
        call expcon (exoidm(i), coord_names2, ierr)
        write (iout, '("after expcon (",i1,"), error = ", i4)')
     1		i, ierr
1003  continue


c
c write element order map
c

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      call expmap (exoid, elem_map, ierr)
      write (iout, '("after expmap, error = ", i4)' ) ierr

      do 12 i = 1, num_elem2
         elem_map2(i) = i
12    continue

      do 1004 i=1,nexofiles
        call expmap (exoidm(i), elem_map2, ierr)
        write (iout, '("after expmap (",i1,"), error = ", i4)')
     1		i, ierr
1004  continue

c
c write element block parameters
c

      num_elem_in_block(1) = 1
      num_elem_in_block(2) = 1

      ebids(1) = 10
      ebids(2) = 11

      cname = "quad"

      call expelb (exoid,ebids(1),cname,num_elem_in_block(1)
     1		,4,1,ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expelb (exoid,ebids(2),cname,num_elem_in_block(2),
     1		 4,1,ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      num_elem_in_block2(1) = 1
      num_elem_in_block2(2) = 1

      ebids2(1) = 10
      ebids2(2) = 11

      cname2 = "quad2"

      do 1005 i=1,nexofiles
        call expelb(exoidm(i),ebids2(1),cname2,num_elem_in_block2(1),
     1		4,1,ierr)
        write (iout, '("after expelb (",i1,"), error = ", i4)')
     1		i, ierr

        call expelb(exoidm(i),ebids2(2),cname2,num_elem_in_block2(2),
     1		4,1,ierr)
        write (iout, '("after expelb (",i1,"), error = ", i4)')
     1		i, ierr
1005  continue

c  write element block properties

      prop_names(1) = "MATL"
      prop_names(2) = "DENSITY"
      call exppn(exoid,EXEBLK,2,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr

      call expp(exoid, EXEBLK, ebids(1), "MATL", 10, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(2), "MATL", 20, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      do 1006 i=1,nexofiles
        call exppn(exoidm(i),EXEBLK,2,prop_names,ierr)
        write (iout, '("after exppn (",i1,"), error = ", i4)')
     1		i, ierr

        call expp(exoidm(i), EXEBLK, ebids(1), "MATL", 10, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr
        call expp(exoidm(i), EXEBLK, ebids(2), "MATL", 20, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr
1006  continue

c
c write element connectivity
c

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

      connect2(1) = 1
      connect2(2) = 2 
      connect2(3) = 3 
      connect2(4) = 4

      do 1007 i=1,nexofiles
        call expelc (exoidm(i), ebids2(1), connect2, ierr)
        write (iout, '("after expelc (",i1,"), error = ", i4)')
     1		i, ierr
1007  continue

      connect2(1) = 5
      connect2(2) = 6 
      connect2(3) = 7 
      connect2(4) = 8

      do 1008 i=1,nexofiles
        call expelc (exoidm(i), ebids2(2), connect2, ierr)
        write (iout, '("after expelc (",i1,"), error = ", i4)')
     1		i, ierr
1008  continue

c
c write element block attributes
c

      attrib(1) = 3.14159
      call expeat (exoid, ebids(1), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib(1) = 6.14159
      call expeat (exoid, ebids(2), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib2(1) = 3.
      do 1009 i=1,nexofiles
        call expeat (exoidm(i), ebids2(1), attrib2, ierr)
        write (iout, '("after expeat (",i1,"), error = ", i4)')
     1		i, ierr
1009  continue

      attrib2(1) = 6.
      do 1010 i=1,nexofiles
        call expeat (exoidm(i), ebids2(2), attrib2, ierr)
        write (iout, '("after expeat (",i1,"), error = ", i4)')
     1		i, ierr
1010  continue

c
c write individual node sets
c

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

      node_list2(1) = 2100 
      node_list2(2) = 2101 
      node_list2(3) = 2102 
      node_list2(4) = 2103 
      node_list2(5) = 2104 

      dist_fact2(1) = 21.0 
      dist_fact2(2) = 22.0 
      dist_fact2(3) = 23.0
      dist_fact2(4) = 24.0 
      dist_fact2(5) = 25.0

      do 1011 i=1,nexofiles
        call expnp (exoidm(i), 20, 5, 5, ierr)
        write (iout, '("after expnp (",i1,"), error = ", i4)')
     1		i, ierr

        call expns (exoidm(i), 20, node_list, ierr)
        write (iout, '("after expns (",i1,"), error = ", i4)')
     1		i, ierr
        call expnsd (exoidm(i), 20, dist_fact, ierr)
        write (iout, '("after expnsd (",i1,"), error = ", i4)')
     1		i, ierr

        call expnp (exoidm(i), 21, 3, 3, ierr)
        write (iout, '("after expnp (",i1,"), error = ", i4)')
     1		i, ierr
1011  continue

      node_list2(1) = 2200 
      node_list2(2) = 2201 
      node_list2(3) = 2202 
   
      dist_fact2(1) = 21.1 
      dist_fact2(2) = 22.1 
      dist_fact2(3) = 23.1

      do 1012 i=1,nexofiles
        call expns (exoidm(i), 21, node_list, ierr)
        write (iout, '("after expns (",i1,"), error = ", i4)')
     1		i, ierr
        call expnsd (exoidm(i), 21, dist_fact, ierr)
        write (iout, '("after expnsd (",i1,"), error = ", i4)')
     1		i, ierr
1012  continue

c
c write concatenated node sets; this produces the same information as
c the above code which writes individual node sets
c

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
c

      do 1013 i=1,nexofiles
        prop_names(1) = "FACE"
        call expp(exoidm(i), EXNSET, 20, prop_names(1), 4, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr

        call expp(exoidm(i), EXNSET, 21, prop_names(1), 5, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr

        prop_array(1) = 1000
        prop_array(2) = 2000

        prop_names(1) = "VELOCITY"
        call exppa(exoidm(i), EXNSET, prop_names(1), prop_array, ierr)
        write (iout, '("after exppa (",i1,"), error = ", i4)')
     1		i, ierr
1013  continue

c write individual side sets
c

      elem_list(1) = 11
      elem_list(2) = 12

      node_list(1) = 1 
      node_list(2) = 2 
      node_list(3) = 3 
      node_list(4) = 4

      dist_fact(1) = 30.0
      dist_fact(2) = 30.1
      dist_fact(3) = 30.2
      dist_fact(4) = 30.3

      call expsp (exoid, 30, 2, 4, ierr)
      write (iout, '("after expsp, error = ", i4)' ) ierr

      call expss (exoid, 30, elem_list, node_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr

      call expssd (exoid, 30, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i4)' ) ierr

      elem_list(1) = 13
      elem_list(2) = 14

      side_list(1) = 3
      side_list(2) = 4

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


      elem_list2(1) = 11
      elem_list2(2) = 12

      node_list2(1) = 1 
      node_list2(2) = 2 
      node_list2(3) = 3 
      node_list2(4) = 4

      dist_fact2(1) = 1.1 
      dist_fact2(2) = 2.1 
      dist_fact2(3) = 3.1
      dist_fact2(4) = 4.1

      do 1014 i=1,nexofiles
        call expsp (exoidm(i), 30, 2, 4, ierr)
        write (iout, '("after expsp (",i1,"), error = ", i4)')
     1		i, ierr

        call expss (exoidm(i), 30, elem_list2, node_list2, ierr)
        write (iout, '("after expss (",i1,"), error = ", i4)')
     1		i, ierr

        call expssd (exoidm(i), 30, dist_fact2, ierr)
        write (iout, '("after expssd (",i1,"), error = ", i4)')
     1		i, ierr
1014  continue

      elem_list2(1) = 13
      elem_list2(2) = 14

      side_list2(1) = 3
      side_list2(2) = 4

      dist_fact2(1) = 31.0
      dist_fact2(2) = 31.1
      dist_fact2(3) = 31.2
      dist_fact2(4) = 31.3

      do 1015 i=1,nexofiles
        call expsp (exoidm(i), 31, 2, 4, ierr)
        write (iout, '("after expsp (",i1,"), error = ", i3)')
     1		i, ierr

        call expss (exoidm(i), 31, elem_list2, side_list2, ierr)
        write (iout, '("after expss (",i1,"), error = ", i3)')
     1		i, ierr

        call expssd (exoidm(i), 31, dist_fact2, ierr)
        write (iout, '("after expssd (",i1,"), error = ", i3)')
     1		i, ierr
1015  continue

c
c write concatenated side sets; this produces the same information as
c the above code which writes individual side sets
c

      ids(1) = 30
      ids(2) = 31

      num_elem_per_set(1) = 2
      num_elem_per_set(2) = 2

      num_df_per_set(1) = 4
      num_df_per_set(2) = 4

      elem_ind(1) = 1
      elem_ind(2) = 3

      df_ind(1) = 1
      df_ind(2) = 5

      elem_list(1) = 11
      elem_list(2) = 12
      elem_list(3) = 13
      elem_list(4) = 14

      side_list(1) = 1
      side_list(2) = 2
      side_list(3) = 3
      side_list(4) = 4

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

c     call expcss (exoidm(i), ids, num_elem_per_set, num_df_per_set,
c    1             elem_ind, df_ind, elem_list, side_list, dist_fact,
c    2             ierr)
c     write (iout, '("after expcss (",i1,"), error = ", i4)' ) ierr

      prop_names(1) = "COLOR"
      call expp(exoid, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      call expp(exoid, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr


      do 1016 i=1,nexofiles
        prop_names(1) = "COLOR"
        call expp(exoidm(i), EXSSET, 30, prop_names(1), 100, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr

        call expp(exoidm(i), EXSSET, 31, prop_names(1), 101, ierr)
        write (iout, '("after expp (",i1,"), error = ", i4)')
     1		i, ierr
1016  continue


c
c write QA records
c

      num_qa_rec = 2

      qa_record(1,1) = "TESTWTM fortran version"
      qa_record(2,1) = "testwtm"
      qa_record(3,1) = "07/07/93"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "07/07/93"
      qa_record(4,2) = "16:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i4)' ) ierr

      num_qa_rec2 = 2

      qa_record2(1,1) = "TESTWTM fortran version"
      qa_record2(2,1) = "testwtm"
      qa_record2(3,1) = "07/07/93"
      qa_record2(4,1) = "15:41:33"
      qa_record2(1,2) = "FASTQ"
      qa_record2(2,2) = "fastq"
      qa_record2(3,2) = "07/07/93"
      qa_record2(4,2) = "16:41:33"

      do 1017 i=1,nexofiles
        call expqa (exoidm(i), num_qa_rec2, qa_record2, ierr)
        write (iout, '("after expqa (",i1,"), error = ", i4)')
     1		i, ierr
1017  continue


c
c write information records
c

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i4)' ) ierr

      num_info2 = 3

      inform2(1) = "This is the first info record."
      inform2(2) = "This is the second info record."
      inform2(3) = "This is the third info record."

      do 1018 i=1,nexofiles
        call expinf (exoidm(i), num_info2, inform2, ierr)
        write (iout, '("after expinf (",i1,"), error = ", i4)')
     1		i, ierr
1018  continue

c write results variables parameters and names

      num_glo_vars = 1
  
      var_names(1) = "glo_vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "g", num_glo_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_glo_vars2 = 1
  
      var_names2(1) = "glovars2"

      do 1019 i=1,nexofiles
        call expvp (exoidm(i), "g", num_glo_vars2, ierr)
        write (iout, '("after expvp (",i1,"), error = ", i4)')
     1		i, ierr
        call expvan (exoidm(i), "g", num_glo_vars2, var_names2, ierr)
        write (iout, '("after expvan (",i1,"), error = ", i4)')
     1		i, ierr
1019  continue

      num_nod_vars = 2

      var_names(1) = "nod_var0"
      var_names(2) = "nod_var1"

      call expvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_nod_vars2 = 2

      var_names2(1) = "nodvar20"
      var_names2(2) = "nodvar21"

      do 1020 i=1,nexofiles
        call expvp (exoidm(i), "n", num_nod_vars2, ierr)
        write (iout, '("after expvp (",i1,"), error = ", i4)')
     1		i, ierr
        call expvan (exoidm(i), "n", num_nod_vars2, var_names2, ierr)
        write (iout, '("after expvan (",i1,"), error = ", i4)')
     1		i, ierr
1020  continue
   
      num_ele_vars = 3

      var_names(1) = "ele_var0"
      var_names(2) = "ele_var1"
      var_names(3) = "ele_var2"

      call expvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr

      num_ele_vars2 = 3

      var_names2(1) = "elevar20"
      var_names2(2) = "elevar21"
      var_names2(3) = "elevar22"

      do 1021 i=1,nexofiles
        call expvp (exoidm(i), "e", num_ele_vars2, ierr)
        write (iout, '("after expvp (",i1,"), error = ", i4)')
     1		i, ierr
        call expvan (exoidm(i), "e", num_ele_vars2, var_names2, ierr)
        write (iout, '("after expvan (",i1,"), error = ", i4)')
     1		i, ierr
1021  continue
c
c write element variable truth table
c

      k = 0

      do 30 i = 1,num_elem_blk
         do 20 j = 1,num_ele_vars
            truth_tab(j,i) = 1
20       continue
30    continue

      call exgebi (exoid, ebids, ierr)
      write (iout, '("after exgebi, error = ", i4)' ) ierr
      call expvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ierr)
      write (iout, '("after expvtt, error = ", i4)' ) ierr

      do 1022 i=1,nexofiles
        call exgebi (exoidm(i), ebids2, ierr)
        write (iout, '("after exgebi (",i1,"), error = ", i4)')
     1		i, ierr
        call expvtt (exoidm(i),num_elem_blk,num_ele_vars,truth_tab,ierr)
        write (iout, '("after expvtt (",i1,"), error = ", i4)')
     1		i, ierr
1022  continue
c
c for each time step, write the analysis results;
c the code below fills the arrays glob_var_vals, 
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays
c

      whole_time_step = 1
      num_time_steps = 10

      do 110 iii = 1, num_time_steps
        time_value = real(iii)/100
        time_value2 = real(iii)/100
c
c write time value to regular file
c

        call exptim (exoid, whole_time_step, time_value, ierr)
        write (iout, '("after exptim, error = ", i4)' ) ierr

        do 1023 i=1,nexofiles
          call exptim (exoidm(i), whole_time_step, time_value2, ierr)
          write (iout, '("after exptim (",i1,"), error = ", i4)')
     1		i, ierr
1023    continue

c
c write global variables
c

        do 50 j = 1, num_glo_vars
          glob_var_vals(j) = real(j+1) * time_value
50      continue

        call expgv (exoid, whole_time_step, num_glo_vars, 
     1              glob_var_vals, ierr)
        write (iout, '("after expgv, error = ", i4)' ) ierr

        do 1024 i=1,nexofiles
          call expgv (exoidm(i), whole_time_step, num_glo_vars, 
     1              glob_var_vals, ierr)
          write (iout, '("after expgv (",i1,"), error = ", i4)')
     1		i, ierr
1024    continue

c
c write nodal variables
c

        do 70 k = 1, num_nod_vars
          do 60 j = 1, num_nodes

            nodal_var_vals(j) = real(k) + (real(j) * time_value)

60        continue

          call expnv (exoid, whole_time_step, k, num_nodes, 
     1                nodal_var_vals, ierr)
          write (iout, '("after expnv, error = ", i4)' ) ierr

          do 1025 i=1,nexofiles
            call expnv (exoidm(i), whole_time_step, k, num_nodes, 
     1                nodal_var_vals, ierr)
            write (iout, '("after expnv (",i1,"), error = ", i4)')
     1		i, ierr
1025      continue

70      continue

c
c write element variables
c

        do 100 k = 1, num_ele_vars
          do 90 j = 1, num_elem_blk
            do 80 m = 1, num_elem_in_block(j)

              elem_var_vals(m) = real(k+1) + real(j+1) + 
     1                          (real(m)*time_value)

80          continue

            call expev (exoid, whole_time_step, k, ebids(j), 
     1                  num_elem_in_block(j), elem_var_vals, ierr)
            write (iout, '("after expev, error = ", i4)' ) ierr
            do 1026 i=1,nexofiles
              call expev (exoidm(i), whole_time_step, k, ebids(j), 
     1                  num_elem_in_block(j), elem_var_vals, ierr)
              write (iout, '("after expev (",i1,"), error = ", i4)')
     1		i, ierr
1026        continue

90        continue
100     continue

        whole_time_step = whole_time_step + 1

c
c update the data file; this should be done at the end of every time 
c step to ensure that no data is lost if the analysis dies
c
        call exupda (exoid, ierr)
        write (iout, '("after exupda, error = ", i4)' ) ierr
        do 1027 i=1,nexofiles
          call exupda (exoidm(i), ierr)
          write (iout, '("after exupda (",i1,"), error = ", i4)')
     1		i, ierr
1027    continue

110   continue

c
c close the EXODUS files
c
      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      do 1028 i=1,nexofiles
        call exclos (exoidm(i), ierr)
        write (iout, '("after exclos (",i1,"), error = ", i4)')
     1		i, ierr
1028  continue

      stop
      end

