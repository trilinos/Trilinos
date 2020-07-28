C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      program testwt3

c This is a test program for the Fortran binding of the EXODUS II
c database write routines. This test writes GENISIS (geometry)
c data to the history file.

c     08/10/93	V.R. Yarberry - Updated for use with 2.01 API

      include 'exodus_app.inc'

      integer iin, iout
      integer exoid, exoidh, num_dim, num_nodes, num_elem, num_elem_blk
      integer num_elem_in_block(2), num_node_sets
      integer num_side_sets, error
      integer i, j, k, m, elem_map(2), connect(4)
      integer node_list(10), elem_list(10)
      integer ebids(2),ids(2), num_nodes_per_set(2), num_elem_per_set(1)
      integer node_ind(2), elem_ind(1), num_qa_rec, num_info
      integer num_his_vars, num_glo_vars, num_nod_vars, num_ele_vars
      integer truth_tab(3,2)
      integer hist_time_step, whole_time_step, num_time_steps
      integer cpu_word_size, io_word_size

      real hist_var_vals(10), glob_var_vals(10), nodal_var_vals(8)
      real time_value, elem_var_vals(20)
      real x(8), y(8), dummy(1)
      real attrib(1), dist_fact(8)

      character*(MXLNLN) title
      character*(MXSTLN) coord_names(3)
      character*(MXSTLN) cname
      character*(MXSTLN) var_names(3)
      character*(MXSTLN) qa_record(4,2)
      character*(MXLNLN) inform(3)

      logical whole

      data iin /5/, iout /6/

c  create EXODUS II files

      cpu_word_size = 4
      io_word_size = 4

c     first create a "regular" file that contains everything except
c     history variable info

      exoid = excre ("test.exo",
     1               "r", EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo, id: ", i3)') exoid
      write (iout,'("after excre, error = ", i3)') ierr

c     create a "history" file if you will output history variables

      exoidh = excre ("testh.exo",
     1               "h", EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for testh.exo, id: ", i3)') exoidh
      write (iout,'("after excre, error = ", i3)') ierr

c  initialize file with parameters

      title = "This is test 3 - genisis data in history file"
      num_dim = 2
      num_nodes = 8
      num_elem = 2
      num_elem_blk = 2
      num_node_sets = 2
      num_side_sets = 1

      call expini (exoid, title, num_dim, num_nodes,
     1             num_elem, num_elem_blk, num_node_sets,
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i3)' ) ierr

      call expini (exoidh, title, num_dim, num_nodes,
     1             num_elem, num_elem_blk, num_node_sets,
     2             num_side_sets, ierr)

      write (iout, '("after expini (h), error = ", i3)' ) ierr

c  write nodal coordinates values and names to database

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
      write (iout, '("after expcor, error = ", i3)' ) ierr

      call expcor (exoidh, x, y, dummy, ierr)
      write (iout, '("after expcor (h), error = ", i3)' ) ierr

      coord_names(1) = "xcoorjun"
      coord_names(2) = "ycoorjun"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i3)' ) ierr

      call expcon (exoidh, coord_names, ierr)
      write (iout, '("after expcon (h), error = ", i3)' ) ierr

c write element order map

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      call expmap (exoid, elem_map, ierr)
      write (iout, '("after expmap, error = ", i3)' ) ierr

      call expmap (exoidh, elem_map, ierr)
      write (iout, '("after expmap (h), error = ", i3)' ) ierr

c write element block parameters

      num_elem_in_block(1) = 1
      num_elem_in_block(2) = 1

      ebids(1) = 10
      ebids(2) = 11

      cname = "quadjunk"

      call expelb (exoid, ebids(1), cname, num_elem_in_block(1),
     1		 4,1,ierr)
      write (iout, '("after expelb, error = ", i3)' ) ierr

      call expelb (exoid, ebids(2), cname, num_elem_in_block(2),
     1		4,1,ierr)
      write (iout, '("after expelb, error = ", i3)' ) ierr

      call expelb (exoidh, ebids(1), cname, num_elem_in_block(1),
     1		4,1,ierr)
      write (iout, '("after expelb (h), error = ", i3)' ) ierr

      call expelb (exoidh, ebids(2), cname, num_elem_in_block(2),
     1		 4,1,ierr)
      write (iout, '("after expelbi(h), error = ", i3)' ) ierr

c write element connectivity

      connect(1) = 1
      connect(2) = 2
      connect(3) = 3
      connect(4) = 4

      call expelc (exoid, ebids(1), connect, ierr)
      write (iout, '("after expelc, error = ", i3)' ) ierr

      call expelc (exoidh, ebids(1), connect, ierr)
      write (iout, '("after expelci (h), error = ", i3)' ) ierr

      connect(1) = 5
      connect(2) = 6
      connect(3) = 7
      connect(4) = 8

      call expelc (exoid, ebids(2), connect, ierr)
      write (iout, '("after expelc, error = ", i3)' ) ierr

      call expelc (exoidh, ebids(2), connect, ierr)
      write (iout, '("after expelc (h), error = ", i3)' ) ierr

c write element block attributes

      attrib(1) = 3.14159
      call expeat (exoid, ebids(1), attrib, ierr)
      write (iout, '("after expeat, error = ", i3)' ) ierr

      call expeat (exoidh, ebids(1), attrib, ierr)
      write (iout, '("after expeat (h), error = ", i3)' ) ierr

      attrib(1) = 6.14159
      call expeat (exoid, ebids(2), attrib, ierr)
      write (iout, '("after expeat, error = ", i3)' ) ierr

      call expeat (exoidh, ebids(2), attrib, ierr)
      write (iout, '("after expeat (h), error = ", i3)' ) ierr

c write individual node sets

      call expnp (exoid, 20, 5, ierr)
      write (iout, '("after expnp, error = ", i3)' ) ierr

      call expnp (exoidh, 20, 5, ierr)
      write (iout, '("after expnp (h), error = ", i3)' ) ierr

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

      call expns (exoid, 20, node_list, dist_fact, ierr)
      write (iout, '("after expns, error = ", i3)' ) ierr

      call expns (exoidh, 20, node_list, dist_fact, ierr)
      write (iout, '("after expns (h), error = ", i3)' ) ierr

      call expnp (exoid, 21, 3, ierr)
      write (iout, '("after expnp, error = ", i3)' ) ierr

      call expnp (exoidh, 21, 3, ierr)
      write (iout, '("after expnp (h), error = ", i3)' ) ierr

      node_list(1) = 200
      node_list(2) = 201
      node_list(3) = 202

      dist_fact(1) = 1.1
      dist_fact(2) = 2.1
      dist_fact(3) = 3.1

      call expns (exoid, 21, node_list, dist_fact, ierr)
      write (iout, '("after expns, error = ", i3)' ) ierr

      call expns (exoidh, 21, node_list, dist_fact, ierr)
      write (iout, '("after expns (h), error = ", i3)' ) ierr

c write concatenated node sets; this produces the same information as
c the above code which writes individual node sets

c     ids(1) = 20
c     ids(2) = 21

c     num_nodes_per_set(1) = 5
c     num_nodes_per_set(2) = 3

c     node_ind(1) = 1
c     node_ind(2) = 6

c     node_list(1) = 100
c     node_list(2) = 101
c     node_list(3) = 102
c     node_list(4) = 103
c     node_list(5) = 104
c     node_list(6) = 200
c     node_list(7) = 201
c     node_list(8) = 202

c     dist_fact(1) = 1.0
c     dist_fact(2) = 2.0
c     dist_fact(3) = 3.0
c     dist_fact(4) = 4.0
c     dist_fact(5) = 5.0
c     dist_fact(6) = 1.1
c     dist_fact(7) = 2.1
c     dist_fact(8) = 3.1

c     call expcns (exoid, ids, num_nodes_per_set, node_ind, node_list,
c    1        dist_fact, ierr)
c     write (iout, '("after expcns, error = ", i3)' ) ierr

c write individual side sets

      call expsp (exoid, 30, 2, 4, ierr)
      write (iout, '("after expsp, error = ", i3)' ) ierr

      call expsp (exoidh, 30, 2, 4, ierr)
      write (iout, '("after expsp (h), error = ", i3)' ) ierr

      elem_list(1) = 1
      elem_list(2) = 2

      node_list(1) = 1
      node_list(2) = 2
      node_list(3) = 3
      node_list(4) = 4

      dist_fact(1) = 0.0
      dist_fact(2) = 0.0
      dist_fact(3) = 0.0
      dist_fact(4) = 0.0

      call expss (exoid, 30, elem_list, node_list, ierr)
      write (iout, '("after expss, error = ", i3)' ) ierr

      call expssd (exoid, 30, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i3)' ) ierr

      call expss (exoidh, 30, elem_list, node_list, ierr)
      write (iout, '("after expss (h), error = ", i3)' ) ierr

      call expssd (exoidh, 30, dist_fact, ierr)
      write (iout, '("after expssd (h), error = ", i3)' ) ierr

c write concatenated side sets; this produces the same information as
c the above code which writes individual side sets

c      ids(1) = 30

c      num_elem_per_set(1) = 2

c      num_nodes_per_set(1) = 4

c      elem_ind(1) = 1

c      node_ind(1) = 1

c      elem_list(1) = 1
c      elem_list(2) = 2

c      node_list(1) = 1
c      node_list(2) = 2
c      node_list(3) = 3
c      node_list(4) = 4

c      dist_fact(1) = 0.0
c      dist_fact(2) = 0.0
c      dist_fact(3) = 0.0
c      dist_fact(4) = 0.0

c      call expcss (exoid, ids, num_elem_per_set, num_nodes_per_set,
c     1             elem_ind, node_ind, elem_list, node_list, dist_fact,
c     2             ierr)
c      write (iout, '("after expcss, error = ", i3)' ) ierr

c write QA records

      num_qa_rec = 2

      qa_record(1,1) = "PRONTO2D"
      qa_record(2,1) = "pronto2d"
      qa_record(3,1) = "3/10/92"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "2/10/92"
      qa_record(4,2) = "11:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i3)' ) ierr

      call expqa (exoidh, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa (h), error = ", i3)' ) ierr

c write information records

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i3)' ) ierr

      call expinf (exoidh, num_info, inform, ierr)
      write (iout, '("after expinf (h), error = ", i3)' ) ierr

c write results variables parameters and names

      num_his_vars = 1

      var_names(1) = "his_vars"

      call expvp (exoidh, "h", num_his_vars, ierr)
      write (iout, '("after expvp, error = ", i3)' ) ierr
      call expvan (exoidh, "h", num_his_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i3)' ) ierr

      num_glo_vars = 1

      var_names(1) = "glo_vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i3)' ) ierr
      call expvan (exoid, "g", num_glo_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i3)' ) ierr

      num_nod_vars = 2

      var_names(1) = "nod_var0"
      var_names(2) = "nod_var1"

      call expvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '("after expvp, error = ", i3)' ) ierr
      call expvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i3)' ) ierr

      num_ele_vars = 3

      var_names(1) = "ele_var0"
      var_names(2) = "ele_var1"
      var_names(3) = "ele_var2"

      call expvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '("after expvp, error = ", i3)' ) ierr
      call expvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i3)' ) ierr

c write element variable truth table

      k = 0

      do 30 i = 1,num_elem_blk
         do 20 j = 1,num_ele_vars
            truth_tab(j,i) = 1
20       continue
30    continue

      call exgebi (exoid, ebids, ierr)
      write (iout, '("after exgebi, error = ", i3)' ) ierr
      call expvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ebids,
     &             ierr)
      write (iout, '("after expvtt, error = ", i3)' ) ierr

c for each time step, write the analysis results;
c the code below fills the arrays hist_var_vals, glob_var_vals,
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays

      whole = .true.
      hist_time_step = 1
      whole_time_step = 1
      num_time_steps = 10

      do 110 i = 1, num_time_steps
         time_value = real(i)/100

c if history time step

c write time value to history file

         call exptim (exoidh, hist_time_step, time_value, ierr)
         write (iout, '("after exptim, error = ", i3)' ) ierr

c write history variables to history file

         do 40 j = 1, num_his_vars
            hist_var_vals(j) = real(j+1) * time_value
40       continue

         call exphv (exoidh, hist_time_step, num_his_vars,
     1               hist_var_vals, ierr)
         write (iout, '("after exphv, error = ", i3)' ) ierr

         hist_time_step = hist_time_step + 1

c update the history file

         call exupda (exoidh, ierr)
         write (iout, '("after exupda, error = ", i3)' ) ierr

c if whole time step

         if (whole) then

c write time value to regular file

            call exptim (exoid, whole_time_step, time_value, ierr)
            write (iout, '("after exptim, error = ", i3)' ) ierr

c write global variables

            do 50 j = 1, num_glo_vars
               glob_var_vals(j) = real(j+1) * time_value
50          continue

            call expgv (exoid, whole_time_step, num_glo_vars,
     1                  glob_var_vals, ierr)
            write (iout, '("after expgv, error = ", i3)' ) ierr

c write nodal variables

            do 70 k = 1, num_nod_vars
               do 60 j = 1, num_nodes

                  nodal_var_vals(j) = real(k) + (real(j) * time_value)

60             continue

               call expnv (exoid, whole_time_step, k, num_nodes,
     1                     nodal_var_vals, ierr)
               write (iout, '("after expnv, error = ", i3)' ) ierr

70          continue

c write element variables

            do 100 k = 1, num_ele_vars
               do 90 j = 1, num_elem_blk
                  do 80 m = 1, num_elem_in_block(j)

                     elem_var_vals(m) = real(k+1) + real(j+1) +
     1                                  (real(m)*time_value)

80                continue

                  call expev (exoid, whole_time_step, k, ebids(j),
     1                        num_elem_in_block(j), elem_var_vals, ierr)
                  write (iout, '("after expev, error = ", i3)' ) ierr

90             continue
100         continue

            whole_time_step = whole_time_step + 1

c update the data file; this should be done at the end of every time
c step to ensure that no data is lost if the analysis dies

            call exupda (exoid, ierr)
            write (iout, '("after exupda, error = ", i3)' ) ierr

         endif

110   continue

c close the EXODUS files

      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i3)' ) ierr

      call exclos (exoidh, ierr)
      write (iout, '("after exclos, error = ", i3)' ) ierr

      stop
      end

