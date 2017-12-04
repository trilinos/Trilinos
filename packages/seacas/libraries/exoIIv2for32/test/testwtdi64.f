C    Copyright (c) 2005-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.  
C    
C        * Neither the name of NTESS nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

      program testwt
      implicit none
c
c This is a test program for the Fortran binding of the EXODUS II
c database write routines.
c
      include 'exodusII.inc'

      integer*4 iout, ierr
      integer*4 exoid
      integer*8 num_dim,num_nodes,elem_map(5),num_elem
      integer*8 num_elem_blk,numattr(10)
      integer*8 num_elem_in_block(10), num_nodes_per_elem(10)
      integer*8 num_node_sets, num_side_sets
      integer*8 i, j, k, m, connect(10) 
      integer*8 node_list(100), elem_list(100), side_list(100)
      integer*8 ebids(10),ids(10), num_nodes_per_set(10)
      integer*8 num_elem_per_set(10), num_df_per_set(10)
      integer*8 df_ind(10), node_ind(10), elem_ind(10)

      integer*4 num_qa_rec, num_info
      integer*4 num_glo_vars, num_nod_vars, num_ele_vars
      integer*4 truth_tab(3,5)
      integer*4 whole_time_step, num_time_steps
      integer*4 cpu_word_size, io_word_size
      integer*8 prop_array(2)

      real*8 glob_var_vals(100), nodal_var_vals(100) 
      real*8 time_value, elem_var_vals(100)
      real*8 x(100), y(100), z(100)
      real*8 attrib(100), dist_fact(100)

      character*(MXSTLN) coord_names(3)
      character*(MXSTLN) blk_names(5)
      character*(MXSTLN) nset_names(2)
      character*(MXSTLN) sset_names(5)
      character*(MXSTLN) cname
      character*(MXSTLN) var_names(3)
      character*(MXSTLN) qa_record(4,2)
      character*(MXLNLN) inform(3)
      character*(MXSTLN) prop_names(2)
      character*(MXSTLN) attrib_names(1)

      data iout /6/

      call exopts (EXABRT, ierr)
      write (iout,'("after exopts, error = ", i4)') ierr
      cpu_word_size = 8
      io_word_size = 8
c
c  create EXODUS II files 
c
C ... All integers passed through the API and stored on DB will be 64-bit integers
      exoid = excre ("test.exo",
     1  EXCLOB+EX_ALL_INT64_DB+EX_ALL_INT64_API,
     *  cpu_word_size, io_word_size, ierr)
      write (iout,'("after excre for test.exo, id: ", i8)') exoid
      write (iout,'("  cpu word size: ",i4," io word size: ",i4)')
     1                  cpu_word_size, io_word_size
      write (iout,'("after excre, error = ", i4)') ierr
c
c  initialize file with parameters
c

      num_dim = 3
      num_nodes = 26
      num_elem = 5
      num_elem_blk = 5
      num_node_sets = 2
      num_side_sets = 5
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
     1		num_nodes_per_elem(1),numattr(1),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expelb (exoid,ebids(2),cname,num_elem_in_block(2),
     1		num_nodes_per_elem(2),numattr(2),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      cname = "hex"
      call expelb (exoid,ebids(3),cname,num_elem_in_block(3),
     1		num_nodes_per_elem(3),numattr(3),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      cname = "tetra"
      call expelb (exoid,ebids(4),cname,num_elem_in_block(4),
     1		num_nodes_per_elem(4),numattr(4),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      cname = "wedge"
      call expelb (exoid,ebids(5),cname,num_elem_in_block(5),
     1		num_nodes_per_elem(5),numattr(5),ierr)
      write (iout, '("after expelb, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      blk_names(1) = "block_a";
      blk_names(2) = "block_b";
      blk_names(3) = "block_c";
      blk_names(4) = "block_d";
      blk_names(5) = "block_e";

      call expnams(exoid, EXEBLK, num_elem_blk, blk_names, ierr)
      write (iout, '("after expnams, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      
c  write element block properties

      prop_names(1) = "MATL"
      prop_names(2) = "DENSITY"
      call exppn(exoid,EXEBLK,2,prop_names,ierr)
      write (iout, '("after exppn, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expp(exoid, EXEBLK, ebids(1), "MATL", 10, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expp(exoid, EXEBLK, ebids(2), "MATL", 20, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expp(exoid, EXEBLK, ebids(3), "MATL", 30, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expp(exoid, EXEBLK, ebids(4), "MATL", 40, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expp(exoid, EXEBLK, ebids(5), "MATL", 50, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c write element connectivity
c

      connect(1) = 1
      connect(2) = 2 
      connect(3) = 3 
      connect(4) = 4

      call expelc (exoid, ebids(1), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      connect(1) = 5
      connect(2) = 6 
      connect(3) = 7 
      connect(4) = 8

      call expelc (exoid, ebids(2), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      connect(1) = 17
      connect(2) = 18
      connect(3) = 19 
      connect(4) = 20

      call expelc (exoid, ebids(4), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      connect(1) = 21
      connect(2) = 22
      connect(3) = 23
      connect(4) = 24
      connect(5) = 25
      connect(6) = 26

      call expelc (exoid, ebids(5), connect, ierr)
      write (iout, '("after expelc, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c write element block attributes
c
      attrib(1) = 3.14159
      call expeat (exoid, ebids(1), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      attrib(1) = 6.14159
      call expeat (exoid, ebids(2), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expeat (exoid, ebids(3), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expeat (exoid, ebids(4), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expeat (exoid, ebids(5), attrib, ierr)
      write (iout, '("after expeat, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      attrib_names(1) = 'THICKNESS'
      do i=1, 5
        call expean (exoid, ebids(i), 1, attrib_names, ierr)
        write (iout, '("after expean, error = ", i4)' ) ierr
        if (ierr .ne. 0) then
          call exclos(exoid,ierr)
          call exit (0)
        endif
      end do
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

      call expnp (exoid, 20, 5, 5, ierr)
      write (iout, '("after expnp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expns (exoid, 20, node_list, ierr)
      write (iout, '("after expns, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expnsd (exoid, 20, dist_fact, ierr)
      write (iout, '("after expnsd, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      node_list(1) = 200 
      node_list(2) = 201 
      node_list(3) = 202 
   
      dist_fact(1) = 1.1 
      dist_fact(2) = 2.1 
      dist_fact(3) = 3.1

      call expnp (exoid, 21, 3, 3, ierr)
      write (iout, '("after expnp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expns (exoid, 21, node_list, ierr)
      write (iout, '("after expns, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expnsd (exoid, 21, dist_fact, ierr)
      write (iout, '("after expnsd, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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

c     call expcns (exoid, ids, num_nodes_per_set, num_df_per_set,
c    1        node_ind, df_ind, node_list, dist_fact, ierr)
c     write (iout, '("after expcns, error = ", i4)' ) ierr

      nset_names(1) = "nodeset_a1";
      nset_names(2) = "nodeset_b2";

      call expnams(exoid, EXNSET, num_node_sets, nset_names, ierr)
      write (iout, '("after expnams, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      

c     write node set properties

      prop_names(1) = "FACE"
      call expp(exoid, EXNSET, 20, prop_names(1), 4, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expp(exoid, EXNSET, 21, prop_names(1), 5, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      prop_array(1) = 1000
      prop_array(2) = 2000

      prop_names(1) = "VELOCITY"
      call exppa(exoid, EXNSET, prop_names(1), prop_array, ierr)
      write (iout, '("after exppa, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c write individual side sets
c

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expss (exoid, 30, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expssd (exoid, 30, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      write (iout, '("after expsp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expss (exoid, 31, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expssd (exoid, 31, dist_fact, ierr)
      write (iout, '("after expssd, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expss (exoid, 32, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expss (exoid, 33, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expss (exoid, 34, elem_list, side_list, ierr)
      write (iout, '("after expss, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif


c write concatenated side sets; this produces the same information as
c the above code which writes individual side sets
c

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
c    1		node_ind, elem_list, node_list, side_list, ierr)
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

      prop_names(1) = "COLOR"
      call expp(exoid, EXSSET, 30, prop_names(1), 100, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      call expp(exoid, EXSSET, 31, prop_names(1), 101, ierr)
      write (iout, '("after expp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

      sset_names(1) = "surf_first"
      sset_names(2) = "surf_second";
      sset_names(3) = "surf_third";
      sset_names(4) = "surf_fourth";
      sset_names(5) = "surf_fifth";

      call expnams(exoid, EXSSET, num_side_sets, sset_names, ierr)
      write (iout, '("after expnams, error = ", i4)' ) ierr
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

c write results variables parameters and names

      num_glo_vars = 1
  
      var_names(1) = "glo_vars"

      call expvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expvan (exoid, "g", num_glo_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif


      num_nod_vars = 2

      var_names(1) = "nod_var0"
      var_names(2) = "nod_var1"

      call expvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

   
      num_ele_vars = 3

      var_names(1) = "ele_var0"
      var_names(2) = "ele_var1"
      var_names(3) = "ele_var2"

      call expvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '("after expvp, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif
      call expvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '("after expvan, error = ", i4)' ) ierr
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

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
      if (ierr .ne. 0) then
         call exclos(exoid,ierr)
         call exit (0)
      endif

c
c for each time step, write the analysis results;
c the code below fills the arrays glob_var_vals, 
c nodal_var_vals, and elem_var_vals with values for debugging purposes;
c obviously the analysis code will populate these arrays
c

      whole_time_step = 1
      num_time_steps = 10

      do 110 i = 1, num_time_steps
        time_value = real(i)/100.
c
c write time value
c

        call exptim (exoid, whole_time_step, time_value, ierr)
        write (iout, '("after exptim, error = ", i4)' ) ierr
        if (ierr .ne. 0) then
           call exclos(exoid,ierr)
           call exit (0)
        endif

c
c write global variables
c

        do 50 j = 1, num_glo_vars
          glob_var_vals(j) = real(j+1) * time_value
50      continue

        call expgv (exoid, whole_time_step, num_glo_vars, 
     1              glob_var_vals, ierr)
        write (iout, '("after expgv, error = ", i4)' ) ierr
        if (ierr .ne. 0) then
           call exclos(exoid,ierr)
           call exit (0)
        endif

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
          if (ierr .ne. 0) then
             call exclos(exoid,ierr)
             call exit (0)
          endif

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
            if (ierr .ne. 0) then
               call exclos(exoid,ierr)
               call exit (0)
            endif

90        continue
100     continue

        whole_time_step = whole_time_step + 1

c
c update the data file; this should be done at the end of every time 
c step to ensure that no data is lost if the analysis dies
c
        call exupda (exoid, ierr)
        write (iout, '("after exupda, error = ", i4)' ) ierr
        if (ierr .ne. 0) then
           call exclos(exoid,ierr)
           call exit (0)
        endif

110   continue

c
c close the EXODUS files
c
      call exclos (exoid, ierr)
      write (iout, '("after exclos, error = ", i4)' ) ierr

      stop
      end

