      program testrd1

c
c This is a test program for the Fortran binding of the EXODUS II
c database read routines
c
c	09/07/93 V.R. Yarberry - Modified for API 2.00
      implicit none

      include 'exodusII.inc'

      integer iin, iout, ierr
      integer exoid, num_dim, num_nodes, num_elem, num_elem_blk
      integer num_node_sets
      integer num_side_sets
      integer i, j, k, elem_map(10), connect(10), node_list(100) 
      integer elem_list(100), side_list(100), ids(10), node_map(100) 
      integer num_elem_per_set(10), num_nodes_per_set(10)
      integer num_df_per_set(10)
      integer num_df_in_set, num_sides_in_set
      integer df_ind(10),node_ind(10),elem_ind(10),num_qa_rec,num_info
      integer num_glo_vars, num_nod_vars, num_ele_vars
      integer truth_tab(3,7)
      integer num_time_steps
      integer num_elem_in_block(10), num_nodes_per_elem(10)
      integer num_attr(10), node_ctr_list(10), node_ctr
      integer num_nodes_in_set, num_elem_in_set
      integer df_list_len, list_len, elem_list_len, node_list_len
      integer node_num, time_step, var_index, beg_time, end_time
      integer elem_num
      integer cpu_ws,io_ws
      integer num_props, prop_value, n_node_map, n_elem_map
      integer ival
      integer ioff

      real time_value, time_values(100), var_values(100)
      real x(100), y(100), z(100)
      real attrib(100), dist_fact(100)
      real vers, fdum

      character*(MXSTLN) coord_names(3), qa_record(4,2), var_names(3)
      character*(MXLNLN) inform(3), titl
      character typ*(MXSTLN), cdum*1
      character*(MXSTLN) prop_names(3)
      character*(MXSTLN) attrib_names(100)
      data iin /5/, iout /6/

c      call exopts (EXVRBS, ierr)
c      write (iout, '(/"after exopts, error = ", i3)' ) ierr

c
c open EXODUS II files
c

      cpu_ws = 0
      io_ws = 0

      exoid = exopen ("test.exo", EXREAD, cpu_ws, io_ws, vers, ierr)
      write (iout, '(/"after exopen, error = ",i3)')
     1			ierr

      write (iout, '("test.exo is an EXODUSII file; version ",
     1                f4.2)') vers
      write (iout, '("  I/O word size",i2)') io_ws

c
c read database parameters
c

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


c
c read nodal coordinates values and names from database
c

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

c
c read map parameters
c

      call exgmp (exoid, n_node_map, n_elem_map, ierr)
      write (iout, '(/"after exgmp, error = ", i3)' ) ierr

c
c read element maps
c
 
      call exgpa (exoid, EXEMAP, "ID", ids, ierr)
      write (iout, '(/"after exgpa, error = ", i3)' ) ierr

      do 26 j = 1,n_elem_map

         call exgem (exoid, ids(j), elem_map, ierr)
         write (iout, '(/"after exgem, error = ", i3)' ) ierr
 
         write (iout, '(/"element map id = ", i3)' ) ids(j)

         do 24 i = 1, num_elem
            write (iout, '("elem_map(",i3,") = ", i3)') i, elem_map(i)
C .. test reading partial element maps...
            call exgpem(exoid, ids(j), i, 1, ival, ierr)
            if (ival .ne. elem_map(i)) then
              write (iout, '("bad value found for map ",
     *          " element ", i3)') i, elem_map(i)
            end if
24       continue

26    continue

c
c     read element map properties
c

      call exinq (exoid, EXNEMP, num_props, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
      write (iout,
     1	'(/"There are ",i2," properties for each element map")')
     2	num_props


      call exgpn(exoid, EXEMAP, prop_names, ierr)
      write (iout, '("after exgpn, error = ", i3)' ) ierr

      do 30 i = 1, num_props
        do 28 j = 1, n_elem_map
          call exgp(exoid, EXEMAP,ids(j),prop_names(i),prop_value,ierr)
          if (ierr .eq. 0) then
            write( iout,
     1	      '("elem map ",i2," property(",i2,"): ",a," = ",i5)' )
     2	      j, i, prop_names(i), prop_value
          else
	    write (iout, '(/"after exgp, error = ", i3)' ) ierr
          endif
28      continue
30    continue

c
c read node maps
c
 
      call exgpa (exoid, EXNMAP, "ID", ids, ierr)
      write (iout, '(/"after exgpa, error = ", i3)' ) ierr

      do 34 j = 1,n_node_map

         call exgnm (exoid, ids(j), node_map, ierr)
         write (iout, '(/"after exgnm, error = ", i3)' ) ierr
 
         write (iout, '(/"node map id = ", i3)' ) ids(j)

         do 32 i = 1, num_nodes
            write (iout, '("node_map(",i3,") = ", i3)') i, node_map(i)
32       continue

34    continue

c
c     read node map properties
c

      call exinq (exoid, EXNNMP, num_props, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
      write (iout,
     1	'(/"There are ",i2," properties for each node map")')
     2	num_props


      call exgpn(exoid, EXNMAP, prop_names, ierr)
      write (iout, '("after exgpn, error = ", i3)' ) ierr

      do 38 i = 1, num_props
        do 36 j = 1, n_node_map
          call exgp(exoid, EXNMAP,ids(j),prop_names(i),prop_value,ierr)
          if (ierr .eq. 0) then
            write( iout,
     1	      '("node map ",i2," property(",i2,"): ",a," = ",i5)' )
     2	      j, i, prop_names(i), prop_value
          else
	    write (iout, '(/"after exgp, error = ", i3)' ) ierr
          endif
36      continue
38    continue

c
c read element block parameters
c
c
      call exgebi (exoid, ids, ierr)
      write (iout, '(/"after exgebi, error = ", i3)' ) ierr

      do 40 i = 1, num_elem_blk

         call exgelb (exoid, ids(i), typ, num_elem_in_block(i),
     1                num_nodes_per_elem(i), num_attr(i), ierr)
         write (iout, '(/"after exgelb, error = ", i3)' ) ierr

         write (iout, '("element block id = ", i2,/
     1                  "element type = ", a9,/
     2                  "num_elem_in_block = ", i2,/
     3                  "num_nodes_per_elem = ", i2,/
     4                  "num_attr = ", i2)')
     5                  ids(i), typ, num_elem_in_block(i), 
     6                  num_nodes_per_elem(i), num_attr(i)

40    continue

c     read element block properties

      call exinq (exoid, EXNEBP, num_props, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
      write (iout,
     1	'(/"There are ",i2," properties for each element block")')
     2	num_props


      call exgpn(exoid, EXEBLK, prop_names, ierr)
      write (iout, '("after exgpn, error = ", i3)' ) ierr

      do 47 i = 1, num_props
        do 45 j = 1, num_elem_blk
          call exgp(exoid, EXEBLK,ids(j),prop_names(i),prop_value,ierr)
          if (ierr .eq. 0) then
            write( iout,
     1	      '("elem block ",i2," property(",i2,"): ",a," = ",i5)' )
     2	      j, i, prop_names(i), prop_value
          else
	    write (iout, '(/"after exgp, error = ", i3)' ) ierr
          endif
45      continue
47    continue

c
c read element connectivity
c

      do 60 i = 1, num_elem_blk

         call exgelc (exoid, ids(i), connect, ierr)
         write (iout, '(/"after exgelc, error = ", i3)' ) ierr
 
         write (iout, '("connect array for elem block ", i2)') ids(i)

         do 50 j = 1, num_nodes_per_elem(i)
            write (iout, '(i3)') connect(j)
50       continue

60    continue

c
c read element block attributes
c

      do 70 i = 1, num_elem_blk

         call exgeat (exoid, ids(i), attrib, ierr)
         write (iout, '(/"after exgeat, error = ", i3)' ) ierr

         call exgean (exoid, ids(i), num_attr(i), attrib_names, ierr)
         write (iout, '(/"after exgean, error = ", i3)' ) ierr
 
         write (iout,
     *     '("element block ", i2, " has ",i2," attribute(s) and ",
     *     i2, " element(s):")')
     1          ids(i), num_attr(i), num_elem_in_block(i)
         do j=1, num_attr(i)
           write (iout, 69) attrib_names(j),
     *       (attrib(k),k= j, num_attr(i)*num_elem_in_block(i),
     *       num_attr(i))
         end do
 69      format(A32," = ", 10(f6.4,2x))
70    continue

c
c read individual node sets
c
      if (num_node_sets .gt. 0) then
         call exgnsi (exoid, ids, ierr)
         write (iout, '(/"after exgnsi, error = ", i3)' ) ierr
      endif

      do 100 i = 1, num_node_sets

         call exgnp (exoid, ids(i), num_nodes_in_set, 
     1               num_df_in_set, ierr)
         write (iout, '(/"after exgnp, error = ", i3)' ) ierr

         write (iout, '(/"node set ", i2, " parameters: ",/
     2                  "num_nodes = ", i2)') ids(i), num_nodes_in_set

         call exgns (exoid, ids(i), node_list, ierr)
         write (iout, '(/"after exgns, error = ", i3)' ) ierr

         if (num_df_in_set .gt. 0) then
           call exgnsd (exoid, ids(i), dist_fact, ierr)
           write (iout, '(/"after exgnsd, error = ", i3)' ) ierr
	 endif

         write (iout, '(/"node list for node set ", i2)') ids(i)
 
         do 80 j = 1, num_nodes_in_set
            write (iout, '(i3)') node_list(j)
80       continue

	 if (num_df_in_set .gt. 0) then
           write (iout, '("dist factors for node set ", i2)') ids(i)
           do 90 j = 1, num_nodes_in_set
             write (iout, '(f5.2)') dist_fact(j)
90         continue
	 else
           write (iout, '("no dist factors for node set ", i2)') ids(i)
	 endif

100   continue

c     read node set properties

      call exinq (exoid, EXNNSP, num_props, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
      write (iout,
     1	'(/"There are ",i2," properties for each node set")')
     2	num_props


      call exgpn(exoid, EXNSET, prop_names, ierr)
      write (iout, '("after exgpn, error = ", i3)' ) ierr

      do 107 i = 1, num_props
        do 105 j = 1, num_node_sets
          call exgp(exoid,EXNSET,ids(j),prop_names(i),prop_value,ierr)
          if (ierr .eq. 0) then
            write( iout,
     1	      '("node set ",i2," property(",i2,"): ",a," = ",i5)' )
     2	      j, i, prop_names(i), prop_value
          else
	    write (iout, '(/"after exgp, error = ", i3)' ) ierr
          endif
105      continue
107    continue

c
c read concatenated node sets; this produces the same information as
c the above code which reads individual node sets
c
      call exinq (exoid, EXNODS, num_node_sets, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
 
      if (num_node_sets .gt. 0) then
         call exinq (exoid, EXNSNL, list_len, fdum, cdum, ierr)
         write(iout,'(/"after EXNSNL =",i3," exinq, error = ",i3)')
     1			list_len,ierr
 
         call exinq (exoid, EXNSDF, list_len, fdum, cdum, ierr)
         write(iout,'(/"after EXNSDF =",i3," exinq, error = ",i3)')
     1			list_len,ierr
 
         call exgcns (exoid, ids, num_nodes_per_set, num_df_per_set,
     1                node_ind, df_ind, node_list, dist_fact, ierr)
         write (iout, '(/"after exgcns, error = ", i3)' ) ierr
 
         write (iout, '(/"concatenated node set info")')
 
         write (iout, '("ids = ")')

         do 110 i = 1, num_node_sets
            write (iout, '(i3)') ids(i)
110      continue

         write (iout, '("num_nodes_per_set = ")')

         do 120 i = 1, num_node_sets
            write (iout, '(i3)') num_nodes_per_set(i)
120      continue

         write (iout, '("node_ind = ")')

         do 130 i = 1, num_node_sets
            write (iout, '(i3)') node_ind(i)
130      continue

         write (iout, '("node_list = ")')

         do 140 i = 1, list_len
            write (iout, '(i3)') node_list(i)
140      continue

         write (iout, '("dist_fact = ")')

         do 150 i = 1, list_len
            write (iout, '(f5.3)') dist_fact(i)
150      continue
      endif

c
c read individual side sets
c
      if (num_side_sets .gt. 0) then
         call exgssi (exoid, ids, ierr)
         write (iout, '(/"after exgssi, error = ", i3)' ) ierr
      endif

      do 190 i = 1, num_side_sets

         call exgsp (exoid, ids(i), num_sides_in_set, num_df_in_set,
     1               ierr)
         write (iout, '(/"after exgsp, error = ", i3)' ) ierr
 
         write (iout, '("side set ", i2, " parameters:",/
     2          "num_sides = ", i3,/
     3          "num_dist_factors = ", i3)')
     4          ids(i), num_sides_in_set, num_df_in_set
 
         call exgss (exoid, ids(i), elem_list, side_list, ierr)
         write (iout, '(/"after exgss, error = ", i3)' ) ierr

         call exgssn (exoid, ids(i), node_ctr_list, node_list, ierr)
         write (iout, '(/"after exgssn, error = ", i3)' ) ierr

	 if (num_df_in_set .gt. 0) then
           call exgssd (exoid, ids(i), dist_fact, ierr)
           write (iout, '(/"after exgssd, error = ", i3)' ) ierr
         endif

         write (iout, '(/"element list for side set ", i2)') ids(i)

         num_elem_in_set = num_sides_in_set
         do 160 j = 1, num_elem_in_set
            write (iout, '(i3)') elem_list(j)
160      continue

         write (iout, '("side list for side set ", i2)') ids(i)

         do 170 j = 1, num_sides_in_set
            write (iout, '(i3)') side_list(j)
170      continue

	 node_ctr = 0
	 write (iout, '("node list for side set ", i2)') ids(i)
	 do 178 k=1, num_elem_in_set
	   do 175 j=1, node_ctr_list(k)
	     write (iout, '(i3)') node_list(j+node_ctr)
175	   continue
	   node_ctr = node_ctr+node_ctr_list(k)
178	 continue

	 if (num_df_in_set .gt. 0) then
           write (iout, '("dist factors for side set ", i2)') ids(i)
           do 180 j = 1, num_df_in_set
             write (iout, '(f6.3)') dist_fact(j)
180        continue
	 else
           write (iout, '("no dist factors for side set ", i2)') ids(i)
         endif

190   continue

c     read side set properties

      call exinq (exoid, EXNSSP, num_props, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
      write (iout,
     1	'(/"There are ",i2," properties for each side set")')
     2	num_props


      call exgpn(exoid, EXSSET, prop_names, ierr)
      write (iout, '("after exgpn, error = ", i3)' ) ierr

      do 197 i = 1, num_props
        do 195 j = 1, num_side_sets
          call exgp(exoid, EXSSET,ids(j),prop_names(i),prop_value,ierr)
          if (ierr .eq. 0) then
            write( iout,
     1	      '("side set ",i2," property(",i2,"): ",a," = ",i5)' )
     2	      j, i, prop_names(i), prop_value
          else
	    write (iout, '(/"after exgp, error = ", i3)' ) ierr
          endif
195      continue
197    continue

      call exinq (exoid, EXSIDS, num_side_sets, fdum, cdum, ierr)
      write (iout, '(/"after exinq: EXSIDS =",i3,", error = ",i3)')
     1		num_side_sets,ierr
 
      if (num_side_sets .gt. 0) then
         call exinq (exoid, EXSSEL, elem_list_len, fdum, cdum, ierr)
         write (iout, '(/"after exinq: EXSSEL =",i3,", error = ",i3)')
     1		elem_list_len,ierr
 
         call exinq (exoid, EXSSNL, node_list_len, fdum, cdum, ierr)
         write (iout, '(/"after exinq: EXSSNL =",i3,", error = ",i3)')
     1		node_list_len,ierr
 
         call exinq (exoid, EXSSDF, df_list_len, fdum, cdum, ierr)
         write (iout, '(/"after exinq: EXSSDF =",i3,", error = ",i3)')
     1		df_list_len,ierr
c
c read concatenated side sets; this produces the same information as
c the above code which reads individual side sets
c
         call exgcss (exoid, ids, num_elem_per_set, num_df_per_set, 
     1             elem_ind, df_ind, elem_list, side_list, dist_fact,
     2             ierr)
         write (iout, '(/"after exgcss, error = ", i3)' ) ierr
 
         write (iout, '("concatenated side set info")')
 
         write (iout, '("ids = ")')

         do 200 i = 1, num_side_sets
            write (iout, '(i3)') ids(i)
200      continue

         write (iout, '("num_elem_per_set = ")')

         do 210 i = 1, num_side_sets
            write (iout, '(i3)') num_elem_per_set(i)
210      continue

         write (iout, '("num_df_per_set = ")')

         do 220 i = 1, num_side_sets
            write (iout, '(i3)') num_df_per_set(i)
220      continue

         write (iout, '("elem_ind = ")')

         do 230 i = 1, num_side_sets
            write (iout, '(i3)') elem_ind(i)
230      continue

         write (iout, '("df_ind = ")')

         do 240 i = 1, num_side_sets
            write (iout, '(i3)') df_ind(i)
240      continue

         write (iout, '("elem_list = ")')

         do 250 i = 1, elem_list_len
            write (iout, '(i3)') elem_list(i)
250      continue

         write (iout, '("side_list = ")')

         do 260 i = 1, elem_list_len
            write (iout, '(i3)') side_list(i)
260      continue

         write (iout, '("dist_fact = ")')

         do 270 i = 1, df_list_len
            write (iout, '(f6.3)') dist_fact(i)
270      continue
      endif

275   continue

c
c read QA records
c
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr

      call exgqa (exoid, qa_record, ierr) 
      write (iout, '(/"after exgqa, error = ", i3)' ) ierr
 
      write (iout, '("QA records = ")')

      do 290 i = 1, num_qa_rec
         do 280 j = 1, 4
            write (iout, '(a)') qa_record(j,i)
280      continue
290   continue

c
c read information records
c
      call exinq (exoid, EXINFO, num_info, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr
 
      call exginf (exoid, inform, ierr) 
      write (iout, '(/"after exginf, error = ", i3)' ) ierr

      write (iout, '("info records = ")')

      do 300 i = 1, num_info
         write (iout, '(a81)') inform(i)
300   continue

c
c read global variables parameters and names 
c
      call exgvp (exoid, "g", num_glo_vars, ierr)
      write (iout, '(/"after exgvp, error = ", i3)' ) ierr
 
      call exgvnm (exoid, "g", 1, var_names(1), ierr)
      write (iout, '(/"after exgvnm, error = ", i3)' ) ierr
 
      write (iout, '("There are ",i2," global variables; their names ",
     1                "are :")')  num_glo_vars

      do 320 i = 1, num_glo_vars
         write (iout, '(a9)') var_names(i)
320   continue

c
c read nodal variables parameters and names
c
      call exgvp (exoid, "n", num_nod_vars, ierr)
      write (iout, '(/"after exgvp, error = ", i3)' ) ierr
 
      call exgvan (exoid, "n", num_nod_vars, var_names, ierr)
      write (iout, '(/"after exgvan, error = ", i3)' ) ierr
 
      write (iout, '("There are ",i2," nodal variables; their names ",
     1                "are :")')  num_nod_vars

      do 330 i = 1, num_nod_vars
         write (iout, '(a9)') var_names(i)
330   continue

c
c read element variables parameters and names
c
      call exgvp (exoid, "e", num_ele_vars, ierr)
      write (iout, '(/"after exgvp, error = ", i3)' ) ierr
 
      call exgvan (exoid, "e", num_ele_vars, var_names, ierr)
      write (iout, '(/"after exgvan, error = ", i3)' ) ierr
 
      write (iout, '("There are ",i2," element variables; their names ",
     1                "are :")')  num_ele_vars

      do 340 i = 1, num_ele_vars
         write (iout, '(a9)') var_names(i)
340   continue

c
c read element variable truth table
c
      call exgvtt (exoid, num_elem_blk, num_ele_vars, truth_tab, ierr)
      write (iout, '(/"after exgvtt, error = ", i3)' ) ierr

      write (iout, '("This is the element variable truth table:")')

      do 360 i = 1, num_elem_blk
         do 350 j = 1, num_ele_vars
            write (iout, '(i2)') truth_tab(j,i)
350      continue
360   continue

c
c determine how many time steps are stored
c
      call exinq (exoid, EXTIMS, num_time_steps, fdum, cdum, ierr)
      write (iout, '(/"after exinq, error = ", i3)' ) ierr

      write (iout, '("There are ",i2," time steps in the database.")') 
     1       num_time_steps
c
c read time value at one time step
c
      time_step = 3
      call exgtim (exoid, time_step, time_value, ierr)
      write (iout, '(/"after exgtim, error = ", i3)' ) ierr
 
      write (iout, '("time value at time step ",i2," = ", f5.3)') 
     1       time_step, time_value
c
c read time values at all time steps
c
      call exgatm (exoid, time_values, ierr)
      write (iout, '(/"after exgatm, error = ", i3)' ) ierr
 
      write (iout, '("time values at all time steps are:")')

      do 370 i = 1, num_time_steps
         write (iout, '(f5.3)') time_values(i)
370   continue

      var_index = 1
      beg_time = 1
      end_time = -1
c
c read all global variables at one time step
c
      call exggv (exoid, time_step, num_glo_vars, var_values, ierr)
      write (iout, '(/"after exggv, error = ", i3)' ) ierr
 
      write (iout, '("global variable values at time step ",i2)') 
     1       time_step

      do 400 i = 1, num_glo_vars
         write (iout, '(f5.3)') var_values(i)
400   continue

c
c read a single global variable through time
c
      call exggvt (exoid, var_index, beg_time, end_time, var_values,
     1             ierr)
      write (iout, '(/"after exggvt, error = ", i3)' ) ierr
 
      write (iout, '("global variable ",i2," values through time:")')
     1       var_index

      do 410 i = 1, num_time_steps
         write (iout, '(f5.3)') var_values(i)
410   continue

c
c read a nodal variable at one time step
c
      call exgnv (exoid, time_step, var_index, num_nodes, var_values,
     1            ierr)
      write (iout, '(/"after exgnv, error = ", i3)' ) ierr

      write (iout, '("nodal variable ",i2," values at time step ",i2)') 
     1       var_index, time_step

      do 420 i = 1, num_nodes
         write (iout, '(f5.3)') var_values(i)
420   continue

c
c read a nodal variable through time
c
      node_num = 1

      call exgnvt (exoid, var_index, node_num, beg_time, end_time, 
     1             var_values, ierr)
      write (iout, '(/"after exgnvt, error = ", i3)' ) ierr
 
      write (iout, '("nodal variable ",i2," values for node ",i2,
     1               " through time:")') var_index, node_num

      do 430 i = 1, num_time_steps
         write (iout, '(f5.3)') var_values(i)
430   continue

c
c read an element variable at one time step
c
      call exgebi (exoid, ids, ierr)
      write (iout, '(/"after exgebi, error = ", i3)' ) ierr

      do 450 i = 1, num_elem_blk

         call exgev (exoid, time_step, var_index, ids(i), 
     1               num_elem_in_block(i), var_values, ierr)
         write (iout, '(/"after exgev, error = ", i3)' ) ierr
 
         if (ierr .eq. 0) then
            write (iout, '("element variable ",i2," values of element ",
     1                     "block ",i2," at time step ",i2)')
     2                     var_index, ids(i), time_step
         endif

         do 440 j = 1, num_elem_in_block(i)
            write (iout, '(f6.3)') var_values(j)
440      continue

450   continue

c
c read an element variable through time
c
      var_index = 2
      elem_num = 2

      call exgevt (exoid, var_index, elem_num, beg_time, end_time, 
     1             var_values, ierr)
      write (iout, '(/"after exgevt, error = ", i3)' ) ierr
 
      write (iout, '("element variable ",i2," values for element ",i2,
     1               " through time:")') var_index, elem_num

      do 460 i = 1, num_time_steps
         write (iout, '(f5.3)') var_values(i)
460   continue

      call exclos (exoid, ierr)
      write (iout, '(/"after exclos, error = ", i3)' ) ierr

      stop
      end

