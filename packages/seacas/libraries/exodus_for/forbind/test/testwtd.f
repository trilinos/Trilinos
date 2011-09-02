      program testwtd
c
c This is a test program for the Fortran binding of the EXODUS II
c database write routines using double precision reals.
c

c	history - 
c	Original L.A. Schoof
c	02/25/93 V.R. Yarberry - Added error checks for file creation.
c	03/04/93 V.R. Yarberry - Fixed bug in expvtt test, ebids was not passed 
c	08/31/93 VRY - updated to match API version 2.00
c
      include 'exodusII.inc'

      integer iin, iout
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer num_elem_in_block(2), num_node_sets
      integer num_side_sets
      integer i, j, k, m, elem_map(2), connect(4) 
      integer node_list(10), elem_list(10), side_list(10)
      integer ebids(2),ids(2), num_nodes_per_set(2), num_elem_per_set(2)
      integer num_df_per_set(2)
      integer df_ind(2), node_ind(2), elem_ind(2), num_qa_rec, num_info
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer truth_tab(3,2)
      integer whole_time_step, num_time_steps
      integer cpu_word_size, io_word_size
      integer prop_array(2)

      real*8 glob_var_vals(10), nodal_var_vals(8) 
      real*8 time_value, elem_var_vals(20)
      real*8 x(8), y(8), dummy(1)
      real*8 attrib(1), dist_fact(8)

      character*(MXSTLN) coord_names(3)
      character*(MXSTLN) cname
      character*(MXSTLN) var_names(3)
      character*(MXSTLN) qa_record(4,2)
      character*(MXLNLN) inform(3)
      character*(MXSTLN) prop_names(2)

      logical whole

      data iin /5/, iout /6/

      cpu_word_size = 8
      io_word_size = 8
c
c  create EXODUS II files 
c
      exoid = excre ("test.exo",
     1	 	     EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo,id: ",i4,", err=",i3)')
     1           exoid, ierr
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout,'("after excre, error = ", i4)') ierr
c
c  initialize file with parameters
c

      num_dim = 2
      num_nodes = 8
      num_elem = 2
      num_elem_blk = 2
      num_node_sets = 2
      num_side_sets = 2

      call expini (exoid, "This is a test", num_dim, num_nodes, 
     1             num_elem, num_elem_blk, num_node_sets, 
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i4)' ) ierr

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

      coord_names(1) = "xcoor"
      coord_names(2) = "ycoor"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i4)' ) ierr


c
c write element order map
c

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      call expmap (exoid, elem_map, ierr)
      write (iout, '("after expmap, error = ", i4)' ) ierr

c
c write element block parameters
c

      num_elem_in_block(1) = 1
      num_elem_in_block(2) = 1

      ebids(1) = 10
      ebids(2) = 11

      cname = "quad"

      call expelb (exoid,ebids(1),cname,num_elem_in_block(1),4,1,ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

      call expelb (exoid,ebids(2),cname,num_elem_in_block(2),4,1,ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr

c  write element block properties

      prop_names(1) = "MATL"
      prop_names(2) = "DENSITY"
      call exppn(exoid,EXEBLK,2,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr

      call expp(exoid, EXEBLK, ebids(1), "MATL", 10, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      call expp(exoid, EXEBLK, ebids(2), "MATL", 20, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

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

c
c write element block attributes
c

      attrib(1) = 3.14159
      call expeat (exoid, ebids(1), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

      attrib(1) = 6.14159
      call expeat (exoid, ebids(2), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr

c
c write individual node sets
c

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

c
c write concatenated node sets; this produces the same information as
c the above code which writes individual node sets
c

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

c
c write individual side sets
c

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

      call expcss (exoid, ids, num_elem_per_set, num_df_per_set,
     1             elem_ind, df_ind, elem_list, side_list, dist_fact,
     2             ierr)
      write (iout, '("after expcss, error = ", i4)' ) ierr

      prop_names(1) = "COLOR"
      call expp(exoid, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr

      call expp(exoid, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
c
c
c write QA records
c

      num_qa_rec = 2

      qa_record(1,1) = "TESTWTD fortran version"
      qa_record(2,1) = "testwtd"
      qa_record(3,1) = "07/07/93"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "07/07/93"
      qa_record(4,2) = "16:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i4)' ) ierr


c
c write information records
c

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i4)' ) ierr


c write results variables parameters and names

      num_glo_vars = 1
  
      var_names(1) = "glo_vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      call expvan (exoid, "g", num_glo_vars, var_names, ierr)
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

c
c write element variable truth table
c

      k = 0

      do 30 i = 1,num_elem_blk
         do 20 j = 1,num_ele_vars
            truth_tab(j,i) = 1
20       continue
30    continue
      call expvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ierr)
      write (iout, '("after expvtt, error = ", i4)' ) ierr

c
c for each time step, write the analysis results;
c the code below fills the arrays hist_var_vals, glob_var_vals, 
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays
c

      whole = .true.
      hist_time_step = 1
      whole_time_step = 1
      num_time_steps = 10

      do 110 i = 1, num_time_steps
        time_value = dble(i)/100
c
c write time value
c

        call exptim (exoid, whole_time_step, time_value, ierr)
        write (iout, '("after exptim, error = ", i4)' ) ierr

c
c write global variables
c

        do 50 j = 1, num_glo_vars
          glob_var_vals(j) = real(j+1) * time_value
50      continue

        call expgv (exoid, whole_time_step, num_glo_vars, 
     1              glob_var_vals, ierr)
        write (iout, '("after expgv, error = ", i4)' ) ierr

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

70      continue

c
c write element variables
c

        do 100 k = 1, num_ele_vars
          do 90 j = 1, num_elem_blk
            do 80 m = 1, num_elem_in_block(j)

              elem_var_vals(m) = real(k+1) + real(j+1) + 
     1                          (real(m)*time_value)
c             write(iout,*)'elem_var_val(',m,'): ',elem_var_vals(m)

80          continue

            call expev (exoid, whole_time_step, k, ebids(j), 
     1                  num_elem_in_block(j), elem_var_vals, ierr)
            write (iout, '("after expev, error = ", i4)' ) ierr

90        continue
100     continue

        whole_time_step = whole_time_step + 1

c
c update the data file; this should be done at the end of every time 
c step to ensure that no data is lost if the analysis dies
c
        call exupda (exoid, ierr)
        write (iout, '("after exupda, error = ", i4)' ) ierr

110   continue

c
c close the EXODUS files
c
      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      stop
      end
