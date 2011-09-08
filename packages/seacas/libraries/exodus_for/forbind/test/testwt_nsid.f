      program testwt
c
c This is a test program for the Fortran binding of the EXODUS II
c database write routines.
c
      include 'exodusII.inc'

      integer iin, iout
      integer exoid, num_dim,num_nodes,elem_map(5),num_elem,num_elem_blk
      integer num_elem_in_block(10), num_nodes_per_elem(10),numattr(10)
      integer num_node_sets, num_side_sets
      integer i, j, k, m, connect(37), nnpe(10)
      integer node_list(100), elem_list(100), side_list(100)
      integer ebids(10),ids(10), num_nodes_per_set(10)
      integer num_elem_per_set(10), num_df_per_set(10)
      integer df_ind(10), node_ind(10), elem_ind(10)
      integer num_qa_rec, num_info
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer truth_tab(3,5)
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
      character*(MXSTLN) attrib_names(1)

      data iin /5/, iout /6/

      call exopts (EXABRT, ierr)
      write (iout,'("after exopts, error = ", i4)') ierr
      cpu_word_size = 0
      io_word_size = 0
c
c  create EXODUS II files 
c
      exoid = excre ("test-nsided.exo",
     1	 	     EXCLOB, cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test-nsided.exo, id: ", i4)') exoid
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout,'("after excre, error = ", i4)') ierr
c
c  initialize file with parameters
c
      num_dim = 3
      num_nodes = 33
      num_elem = 7
      num_elem_blk = 1
      num_node_sets = 0
      num_side_sets = 0

      call expini (exoid, "This is a test", num_dim, num_nodes, 
     1             num_elem, num_elem_blk, num_node_sets, 
     2             num_side_sets, ierr)

      write (iout, '("after expini, error = ", i4)' ) ierr

      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c  write nodal coordinates values and names to database
c
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

C Tetra #2 
      x(27) =  2.7
      x(28) =  6.0
      x(29) =  5.7
      x(30) =  3.7

      y(27) =  1.7
      y(28) =  1.7
      y(29) =  1.7
      y(30) =  0.0

      z(27) =  2.7
      z(28) =  3.3
      z(29) =  1.7
      z(30) =  2.3

C 3d Tri 
      x(31) =  0.0
      x(32) = 10.0
      x(33) = 10.0

      y(31) =  0.0
      y(32) =  0.0
      y(33) = 10.0

      z(31) =  0.0
      z(32) =  0.0
      z(33) = 10.0

      call expcor (exoid, x, y, z, ierr)
      write (iout, '("after expcor, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      coord_names(1) = "xcoor"
      coord_names(2) = "ycoor"
      coord_names(3) = "zcoor"

      call expcon (exoid, coord_names, ierr)
      write (iout, '("after expcon, error = ", i4)' ) ierr
      call exupda(exoid,ierr)
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif


c
c write element order map
c

      do 10 i = 1, num_elem
         elem_map(i) = i
10    continue

      call expmap (exoid, elem_map, ierr)
      write (iout, '("after expmap, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c write element block parameters
c

      num_elem_in_block(1) = 7

      num_nodes_per_elem(1) = 37 ! This is total nodes per block

      ebids(1) = 10

      numattr(1) = 0

      cname = "nsided"

      call expelb (exoid,ebids(1),cname,num_elem_in_block(1),
     1		num_nodes_per_elem(1),numattr(1),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c write element connectivity
c
      connect( 1) = 1
      connect( 2) = 2 
      connect( 3) = 3 
      connect( 4) = 4
      nnpe(1) = 4

      connect( 5) = 5
      connect( 6) = 6 
      connect( 7) = 7 
      connect( 8) = 8
      nnpe(2) = 4

      connect( 9) =  9
      connect(10) = 10
      connect(11) = 11 
      connect(12) = 12
      connect(13) = 13
      connect(14) = 14
      connect(15) = 15
      connect(16) = 16
      nnpe(3) = 8

      connect(17) = 17
      connect(18) = 18
      connect(19) = 19 
      connect(20) = 20
      nnpe(4) = 4

      connect(21) = 21
      connect(22) = 22
      connect(23) = 23
      connect(24) = 24
      connect(25) = 25
      connect(26) = 26
      nnpe(5) = 6

      connect(27) = 17
      connect(28) = 18
      connect(29) = 19
      connect(30) = 20
      connect(31) = 27
      connect(32) = 28
      connect(33) = 30
      connect(34) = 29
      nnpe(6) = 8

      connect(35) = 31
      connect(36) = 32
      connect(37) = 33;
      nnpe(7) = 3

      call expelc (exoid, ebids(1), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expecpp(exoid, EXEBLK, ebids(1), nnpe, ierr)
      write (iout, '("after expecpp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
c
c
c write QA records
c

      num_qa_rec = 2

      qa_record(1,1) = "TESTWT fortran version"
      qa_record(2,1) = "testwt"
      qa_record(3,1) = "07/07/93"
      qa_record(4,1) = "15:41:33"
      qa_record(1,2) = "FASTQ"
      qa_record(2,2) = "fastq"
      qa_record(3,2) = "07/07/93"
      qa_record(4,2) = "16:41:33"

      call expqa (exoid, num_qa_rec, qa_record, ierr)
      write (iout, '("after expqa, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif


c
c write information records
c

      num_info = 3

      inform(1) = "This is the first information record."
      inform(2) = "This is the second information record."
      inform(3) = "This is the third information record."

      call expinf (exoid, num_info, inform, ierr)
      write (iout, '("after expinf, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c ... Define and write some coordinate frames
      call putfrm(exoid)

c
c close the EXODUS files
c
      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      stop
      end

      subroutine putfrm(exoid)
      implicit none
      include 'exodusII.inc'

      integer exoid, ierr, i
      integer numfrm;   ! Assumed to be 3 for remaining dimensions
      integer cfids(3), tags(3)
      real    coord(27)

      numfrm = 3
      
      cfids(1) =   1
      cfids(2) =  11
      cfids(3) = 111

      tags(1) = EXCFREC
      tags(2) = EXCFCYL
      tags(3) = EXCFSPH

! NOTE: These values may not be sensical; just used for testing.
      do i=0,2
        COORD(9*i+1) = i+0.1
        COORD(9*i+2) = i+0.2
        COORD(9*i+3) = i+0.3
        COORD(9*i+4) = i+0.4
        COORD(9*i+5) = i+0.5
        COORD(9*i+6) = i+0.6
        COORD(9*i+7) = i+0.7
        COORD(9*i+8) = i+0.8
        COORD(9*i+9) = i+0.9
      end do
      
      call expfrm(exoid, numfrm, cfids, coord, tags, ierr);
      write (6,'("after expfrm, error = ", i4)') ierr

      return
      end
