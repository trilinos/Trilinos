!! 
!! @HEADER
!! *****************************************************************************
!!  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
!!
!! Copyright 2012 NTESS and the Zoltan contributors.
!! SPDX-License-Identifier: BSD-3-Clause
!! *****************************************************************************
!! @HEADER
!!

module dr_chaco_io
use zoltan
use zoltan_user_data
use mpi_h
use dr_const
use dr_input
implicit none
private

public :: read_chaco_mesh, free_element_arrays, in_list, build_elem_comm_maps, initialize_element

!                                                                          
!--------------------------------------------------------------------------
! Author(s):  Matthew M. St.John (9226)                                    
!   Translated to Fortran by William F. Mitchell
!--------------------------------------------------------------------------
! Revision History:                                                        
!                                                                          
!    24 May 1999:      Date of creation                                    
!       1 September 1999: Translation to Fortran
!--------------------------------------------------------------------------


logical, parameter :: CHECK_INPUT = .false.
integer(Zoltan_INT), parameter :: MAP_ALLOC = 10


contains

!**************************************************************************
!**************************************************************************
!**************************************************************************

logical function read_chaco_mesh(Proc, Num_Proc, prob, pio_info, elements)
integer(Zoltan_INT) :: Proc, Num_Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info
type(ELEM_INFO), pointer :: elements(:)

!   Local declarations. 
  character(len=FILENAME_MAX+8) :: chaco_fname

  integer(Zoltan_INT) :: i, nvtxs, ios, allocstat
  integer(Zoltan_INT) :: ndim = 0
  integer(Zoltan_INT), pointer, dimension(:) :: start, adj, vtxdist

  real(Zoltan_FLOAT), pointer, dimension(:) :: vwgts, ewgts, x, y, z

  integer(Zoltan_INT) :: fp
!**************************** BEGIN EXECUTION *****************************

! Set appropriate callbacks for this file type.
  Test_Graph_Callbacks = 1
  Test_Hypergraph_Callbacks = 0

  nullify(start, adj, vwgts, vtxdist, ewgts, x, y, z)

  if (Proc == 0) then

!     Open and read the Chaco graph file. 
    fp = 12
    chaco_fname = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".graph"
    open(unit=fp,file=chaco_fname,action='read',status='old',iostat=ios)
    if (ios /= 0) then
      print *, "fatal:  Could not open Chaco graph file ", chaco_fname
      read_chaco_mesh = .false.
      return
    endif

!     read the array in on processor 0 
    if (.not.chaco_input_graph(fp, chaco_fname, start, adj, nvtxs, &
                               vwgts, ewgts)) then
      print *, "fatal: Error returned from chaco_input_graph"
      read_chaco_mesh = .false.
      return
    endif

!     Read Chaco geometry file, if provided. 
    fp = 12
    chaco_fname = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".coords"
    open(unit=fp,file=chaco_fname,action='read',status='old',iostat=ios)
    if (ios /= 0) then
      print *, "warning:  Could not open Chaco geometry file ",chaco_fname, &
              "; no geometry data will be read"
    else
!       read the coordinates in on processor 0 
      if (.not.chaco_input_geom(fp, chaco_fname, nvtxs, ndim, x, y, z)) then
        print *, "fatal: Error returned from chaco_input_geom"
        read_chaco_mesh = .false.
        return
      endif
    endif
  endif ! Proc == 0

!   Distribute graph 
  if (.not.chaco_dist_graph(zoltan_get_global_comm(), 0, nvtxs, vtxdist, start, adj, &
                            vwgts, ewgts, ndim, x, y, z)) then
      print *, "fatal: Error returned from chaco_dist_graph"
      read_chaco_mesh = .false.
      return
  endif

!   Initialize Mesh structure for Chaco mesh. 
  Mesh%num_elems = nvtxs
  Mesh%elem_array_len = Mesh%num_elems + 5
  Mesh%num_dims = ndim
  Mesh%num_el_blks = 1

  allocate(Mesh%eb_ids(0:Mesh%num_el_blks-1), &
           Mesh%eb_cnts(0:Mesh%num_el_blks-1), &
           Mesh%eb_nnodes(0:Mesh%num_el_blks-1), &
           Mesh%eb_nattrs(0:Mesh%num_el_blks-1), stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_chaco_mesh = .false.
    return
  endif

  allocate(Mesh%eb_names(0:Mesh%num_el_blks-1),stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_chaco_mesh = .false.
    return
  endif

  Mesh%eb_ids(0) = 1
  Mesh%eb_cnts(0) = nvtxs
!  
!   * Each element has one set of coordinates (i.e., node) if a coords file
!   * was provided; zero otherwise. 
!   
  if (associated(x)) then
    Mesh%eb_nnodes(0) = 1
  else
    Mesh%eb_nnodes(0) = 0
  endif
  Mesh%eb_nattrs(0) = 0
  Mesh%eb_names(0) = "chaco"

!   allocate the element structure array 
  allocate(elements(0:Mesh%elem_array_len-1), stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_chaco_mesh = .false.
    return
  endif

!  
!   * intialize all of the element structs as unused by
!   * setting the globalID to -1
!   
  do i = 0, Mesh%elem_array_len-1
    call initialize_element(elements(i))
  end do

!  
!   * now fill the element structure array with the
!   * information from the Chaco file
!   
  if (.not.fill_elements(Proc, Num_Proc, prob, elements, nvtxs, vtxdist, &
                     start, adj, vwgts, ewgts, ndim, x, y, z)) then
    print *, "fatal: Error returned from fill_elements"
    read_chaco_mesh = .false.
    return
  endif

  if (associated(adj)) deallocate(adj)
  if (associated(vwgts)) deallocate(vwgts)
  if (associated(ewgts)) deallocate(ewgts)
  if (associated(start)) deallocate(start)
  if (associated(vtxdist)) deallocate(vtxdist)
  if (associated(x)) deallocate(x)
  if (associated(y)) deallocate(y)
  if (associated(z)) deallocate(z)

  read_chaco_mesh = .true.
end function read_chaco_mesh

!***************************************************************************
!***************************************************************************
!***************************************************************************

logical function fill_elements(Proc, Num_Proc, prob, elem, nvtxs, vtxdist, &
                               start, adj, vwgts, ewgts, ndim, x, y, z)
  integer(Zoltan_INT) ::  Proc
  integer(Zoltan_INT) ::  Num_Proc
  type(PROB_INFO) :: prob                ! problem description
  type(ELEM_INFO), pointer :: elem(:)    ! array of element information
  integer(Zoltan_INT) ::  nvtxs               ! number of vertices in graph
  integer(Zoltan_INT), pointer ::  vtxdist(:) ! vertex distribution data
  integer(Zoltan_INT), pointer ::  start(:)   ! start of edge list for each vertex
  integer(Zoltan_INT), pointer ::  adj(:)     ! edge list data
  real(Zoltan_FLOAT), pointer  ::  vwgts(:)   ! vertex weight list data
  real(Zoltan_FLOAT), pointer  ::  ewgts(:)   ! edge weight list data
  integer(Zoltan_INT) ::  ndim                ! dimension of the geometry
  real(Zoltan_FLOAT), pointer  ::  x(:)       ! x-coordinates of the vertices
  real(Zoltan_FLOAT), pointer  ::  y(:)       ! y-coordinates of the vertices
  real(Zoltan_FLOAT), pointer  ::  z(:)       ! z-coordinates of the vertices

!   Local declarations. 
  integer(Zoltan_INT) :: i, j, k, start_id, elem_id, local_id, allocstat
!**************************** BEGIN EXECUTION *****************************

  start_id = vtxdist(Proc)+1  ! global ids start at 1 

  do i = 0, Mesh%num_elems-1
    elem(i)%globalID = start_id + i
    if (associated(vwgts)) then
      elem(i)%cpu_wgt = vwgts(i)
    else
      elem(i)%cpu_wgt = 1.0
    endif
    elem(i)%elem_blk = 0        ! only one element block for all vertices 
    elem(i)%my_part = Proc
    elem(i)%perm_value = -1
    elem(i)%invperm_value = -1
    if (Mesh%num_dims > 0) then
!       One set of coords per element. 
      allocate(elem(i)%connect(0:0))
      elem(i)%connect(0) = elem(i)%globalID
      allocate(elem(i)%coord(0:Mesh%num_dims-1,0:0))
      elem(i)%coord(0,0) = x(i)
      if (Mesh%num_dims > 1) then
        elem(i)%coord(1,0) = y(i)
        if (Mesh%num_dims > 2) then
          elem(i)%coord(2,0) = z(i)
        endif
      endif
    endif

!     now start with the adjacencies 
    if (associated(start)) then
      elem(i)%nadj = start(i+1) - start(i)
    else
      elem(i)%nadj = 0
    endif
    if (elem(i)%nadj > 0) then
      elem(i)%adj_len = elem(i)%nadj
      allocate(elem(i)%adj(0:elem(i)%nadj-1), &
               elem(i)%adj_proc(0:elem(i)%nadj-1),stat=allocstat)
      if (allocstat /= 0) then
        print *, "fatal: insufficient memory"
        fill_elements= .false.
        return
      endif
      if (associated(ewgts)) then
        allocate(elem(i)%edge_wgt(0:elem(i)%nadj-1),stat=allocstat)
        if (allocstat /= 0) then
           print *, "fatal: insufficient memory"
           fill_elements= .false.
           return
        endif
      else
        nullify(elem(i)%edge_wgt)
      endif

      do j = 0, elem(i)%nadj-1
        elem_id = adj(start(i) + j)

!         determine which processor the adjacent vertex is on 
        do k = 0, Num_Proc-1
!          Compare with <= since elem_id is 1-based and vtxdist is 0-based. 
          if (elem_id <= vtxdist(k+1)) exit
        end do

!         sanity check 
        if (k == Num_Proc) then
          print *, "fatal:  adjacent element ",elem_id, &
                   " not in vtxdist array ",i,j
          fill_elements = .false.
          return
        endif 

!        
!         * if the adjacent element is on this processor
!         * then find the local id for that element
!         
        if (k == Proc) then
          local_id = elem_id - start_id
          elem(i)%adj(j) = local_id
        else !  use the global id 
          elem(i)%adj(j) = elem_id
        endif

        elem(i)%adj_proc(j) = k

        if (associated(ewgts)) then
          elem(i)%edge_wgt(j) = ewgts(start(i) + j)
        endif
      end do
    endif !  End: "if (elem(i)%nadj > 0)" 
  end do !  End: "do i = 0, Mesh.num_elems-1" 

  if (.not.build_elem_comm_maps(Proc, elem)) then
    print *, "Fatal: error building initial elem comm maps"
    fill_elements = .false.
    return
  endif

  fill_elements = .true.
end function fill_elements

!********************************************************

logical function chaco_input_graph(fin, inname, start, adjacency, nvtxs, vweights, eweights)
integer(Zoltan_INT) :: fin                         ! input file 
character(len=*) :: inname                 ! name of input file 
integer(Zoltan_INT), pointer :: start(:)     ! start of edge list for each vertex
integer(Zoltan_INT), pointer :: adjacency(:) ! edge list data 
integer(Zoltan_INT) :: nvtxs                 ! number of vertices in graph 
real(Zoltan_FLOAT), pointer :: vweights(:)         ! vertex weight list data 
real(Zoltan_FLOAT), pointer :: eweights(:)   ! edge weight list data 
 
integer(Zoltan_INT) :: adjptr                ! loops through adjacency data 
integer(Zoltan_INT) :: ewptr                 ! loops through edge weight data 
integer(Zoltan_INT) :: narcs                ! number of edges expected in graph 
integer(Zoltan_INT) :: nedges        ! twice number of edges really in graph 
integer(Zoltan_INT) :: nedge        ! loops through edges for each vertex 
logical :: found_flag        ! is vertex found in adjacency list? 
logical :: skip_flag        ! should this edge be ignored? 
integer(Zoltan_INT) :: vtx                ! vertex in graph 
integer(Zoltan_INT) :: sum_edges        ! total number of edges read so far 
integer(Zoltan_INT) :: option                ! input option 
logical :: using_ewgts        ! are edge weights in input file? 
logical :: using_vwgts        ! are vertex weights in input file? 
logical :: vtxnums                ! are vertex numbers in input file? 
integer(Zoltan_INT) :: vertex                ! current vertex being read 
logical :: new_vertex        ! new vertex being read 
real(Zoltan_FLOAT) :: weight        ! weight being read 
real(Zoltan_FLOAT) :: eweight        ! edge weight being read 
integer(Zoltan_INT) :: neighbor                ! neighbor of current vertex 
integer(Zoltan_INT) :: self_edge        ! is a self edge encountered? 
logical :: ignore_me        ! is this edge being ignored? 
integer(Zoltan_INT) :: ignored                ! how many edges are ignored? 
logical :: error_flag        ! error reading input? 
integer(Zoltan_INT) :: j                ! loop counters 
integer(Zoltan_INT) :: i ! current data index on input line
integer(Zoltan_INT) :: ints_read(32) ! array of integers from one input line
integer(Zoltan_INT) :: nints_read    ! number of integers on last input line read
real(Zoltan_FLOAT)  :: vals_read(32) ! array of values from one input line
integer(Zoltan_INT) :: nvals_read    ! number of values on last input line read

    nullify(start, adjacency, vweights, eweights)
    error_flag = .false.

!     Read first line  of input (= nvtxs, narcs, option). 
!     The (decimal) digits of the option variable mean:
!           1's digit not zero => input edge weights
!          10's digit not zero => input vertex weights
!         100's digit not zero => include vertex numbers 

    call read_graph_line(fin,ints_read,nints_read)
    if (nints_read < 2) then
        print *,"ERROR in graph file ",inname, &
                ": Less than two integers on the first noncomment line"
        close(fin)
        chaco_input_graph = .false.
        return
    endif

    nvtxs = ints_read(1)
    narcs = ints_read(2)
    if (nints_read > 2) then
        option = ints_read(3)
    else
        option = 0
    endif

    if (nvtxs <= 0) then
        print *,"ERROR in graph file ", inname, &
                ": Invalid number of vertices ", nvtxs
        close(fin)
        chaco_input_graph = .false.
        return
    endif

    if (narcs < 0) then
        print *,"ERROR in graph file ", inname, &
                ": Invalid number of expected edges ", narcs
        close(fin)
        chaco_input_graph = .false.
        return
    endif

    using_ewgts = (option - 10 * (option / 10)) /= 0
    option = option/10
    using_vwgts = (option - 10 * (option / 10)) /= 0
    option = option/10
    vtxnums = (option - 10 * (option / 10)) /= 0

!     Allocate space for rows and columns. 
    allocate(start(0:nvtxs))
    if (narcs /= 0) then
        allocate(adjacency(0:2*narcs))
    endif

    if (using_vwgts) then
        allocate(vweights(0:nvtxs-1))
    endif

    if (using_ewgts) then
        allocate(eweights(0:2*narcs))
    endif

    adjptr = 0
    ewptr = 0
    self_edge = 0
    ignored = 0

    sum_edges = 0
    nedges = 0
    start(0) = 0
    vertex = 0
    vtx = 0
    new_vertex = .true.

! for each input line

  if (narcs > 0) then
    do
        call read_real_line(fin,vals_read,nvals_read)
        i = 1
        if (nvals_read == 0) then ! end of data
           if (vertex == nvtxs) exit
           print *,"ERROR in graph file ", inname, &
                   ": end of data before assigning all vertices"
           close(fin)
           chaco_input_graph = .false.
           return
        endif

! If multiple input lines per vertex, read vertex number. 
        if (vtxnums) then
            j = NINT(vals_read(i))
            i = i+1
            if (j /= vertex .and. j /= vertex + 1) then
                print *,"ERROR in graph file ", inname, &
                        ": out-of-order vertex number ",j
                close(fin)
                chaco_input_graph = .false.
                return
            endif
            if (j /= vertex) then
                new_vertex = .true.
                vertex = j
            else
                new_vertex = .false.
            endif
        else
            vtx = vtx + 1
            vertex = vtx
        endif

        if (vertex > nvtxs) exit

! If vertices are weighted, read vertex weight. 
        if (using_vwgts .and. new_vertex) then
            if (nvals_read < i) then
                print *,"ERROR in graph file ", inname, &
                        ": no weight for vertex ", vertex
                close(fin)
                chaco_input_graph = .false.
                return
            endif
            weight = vals_read(i)
            i = i+1
            if (weight < 0) then
                print *,"ERROR in graph file ", inname, &
                        ": zero or negative weight entered for vertex ", vertex
                close(fin)
                chaco_input_graph = .false.
                return
            endif
            vweights(vertex-1) = weight
        endif

        nedge = 0;

        do
            if (i > nvals_read) exit

! Read number of adjacent vertex. 
            neighbor = NINT(vals_read(i))
            i = i+1

            skip_flag = .false.
            ignore_me = .false.

          if (CHECK_INPUT) then
            if (neighbor > nvtxs) then
                print *,"ERROR in graph file ", inname, &
                        ": nvtxs=",nvtxs,", but edge (",vertex,",",neighbor, &
                        ") was input."
                close(fin)
                chaco_input_graph = .false.
                return
            endif
            if (neighbor < 0) then
                print *,"ERROR in graph file ", inname, &
                        ": zero or negative vertex in edge (", &
                        vertex,",",neighbor,")"
                close(fin)
                chaco_input_graph = .false.
                return
            endif

            if (neighbor == vertex) then
                if ((self_edge==0) .and. CHECK_INPUT) then
                    print *,"WARNING: Self edge (",vertex,",",vertex, &
                            ") being ignored."
                endif
                skip_flag = .true.
                self_edge = self_edge + 1
            endif

!  Check if adjacency is repeated. 
            if (.not. skip_flag) then
                found_flag = .false.
                do j = start(vertex-1), sum_edges+nedge-1
                    if (adjacency(j) == neighbor) then
                        found_flag = .true.
                        exit
                    endif
                end do
                if (found_flag) then
                    print *,"WARNING: Multiple occurences of edge (",vertex, &
                            ",",neighbor,") ignored"
                    skip_flag = .true.
                    if (.not.ignore_me) then
                        ignore_me = .true.
                        ignored = ignored + 1
                    endif
                endif
            endif
          endif !CHECK_INPUT

! Read edge weight if it's being input. 
            if (using_ewgts) then
                if (nvals_read < i) then
                    print *,"ERROR in graph file ", inname, &
                            ": no weight for edge ",vertex,",",neighbor,")."
                    close(fin)
                    chaco_input_graph = .false.
                    return
                endif

                eweight = vals_read(i)
                i = i+1

                if (eweight <= 0 .and. CHECK_INPUT) then
                    print *,"WARNING: Bad weight entered for edge ",vertex, &
                            ",",neighbor,").  Edge ignored."
                    skip_flag = .true.
                    if (.not. ignore_me) then
                        ignore_me = .true.
                        ignored = ignored + 1
                    endif
                else
                    eweights(ewptr) = eweight
                    ewptr = ewptr + 1
                endif
            endif

! Check for edge only entered once. 
            if (neighbor < vertex .and. .not.skip_flag) then
                found_flag = .false.
                do j = start(neighbor-1), start(neighbor)-1
                    if (adjacency(j) == vertex) then
                        found_flag = .true.
                        exit
                    endif
                end do
                if (.not. found_flag) then
                    print *,"ERROR in graph file ", inname, &
                        ": Edge (",vertex,",",neighbor,") entered but not (", &
                            neighbor,",",vertex,")."
                    error_flag = .true.
                endif
            endif

! Add edge to data structure. 
            if (.not. skip_flag) then
                nedges = nedges + 1
                if (nedges > 2*narcs) then
                    print *,"ERROR in graph file ", inname, &
                            ": at least ",nedges," adjacencies entered,", &
                            " but nedges = ",narcs
                    close(fin)
                    chaco_input_graph = .false.
                    return
                endif
                adjacency(adjptr) = neighbor
                adjptr = adjptr + 1
                nedge = nedge + 1
            endif

        end do

        sum_edges = sum_edges + nedge
        start(vertex) = sum_edges
    end do
  endif ! narcs > 0

! Make sure there's nothing else in file. 
    call read_real_line(fin,vals_read,nvals_read)
    if (nvals_read /= 0 .and. CHECK_INPUT) then
        print *,"WARNING: Possible error in graph file ", inname
        print *,"         Data found after expected end of file"
    endif

    start(nvtxs) = sum_edges

    if (self_edge > 1 .and. CHECK_INPUT) then
        print *,"WARNING: ",self_edge," self edges were read and ignored."
    endif

    if (vertex /= 0) then ! Normal file was read. 
!         Make sure narcs was reasonable. 
        if (nedges + 2 * self_edge /= 2 * narcs .and. &
            nedges + 2 * self_edge + ignored /= 2 * narcs .and. &
                nedges + self_edge /= 2 * narcs .and. & 
                nedges + self_edge + ignored /= 2 * narcs .and. & 
                nedges /= 2 * narcs .and. &
                nedges + ignored /= 2 * narcs .and. &
                CHECK_INPUT) then
            print *,"WARNING: I expected ",narcs, &
                    " edges entered twice, but I only count ",nedges
        endif

    else
! Graph was empty => must be using inertial method. 
        deallocate(start)
        if (associated(adjacency)) deallocate(adjacency)
        if (associated(vweights)) deallocate(vweights)
        if (associated(eweights)) deallocate(eweights)
    endif

    close(fin)

    chaco_input_graph = .not. error_flag
end function chaco_input_graph

! This software was developed by Bruce Hendrickson and Robert Leland   *
! * at Sandia National Laboratories under US Department of Energy        *
! * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. 

logical function chaco_input_geom(fingeom, geomname, nvtxs, igeom, x, y, z)
integer(Zoltan_INT) ::     fingeom        ! geometry input file 
character(len=FILENAME_MAX) :: geomname ! name of geometry file 
integer(Zoltan_INT) :: nvtxs        ! number of coordinates to read 
integer(Zoltan_INT) :: igeom        ! dimensionality of geometry 
real(Zoltan_FLOAT), pointer :: x(:), y(:), z(:) ! coordiates of vertices 

    real(Zoltan_FLOAT) :: xc, yc, zc ! first x, y, z coordinate 
    integer(Zoltan_INT) :: nread     ! number of lines of coordinates read 
    integer(Zoltan_INT) :: ndims     ! number of values in an input line 
    integer(Zoltan_INT) :: i         ! loop counter 
    real(Zoltan_FLOAT) :: floats_read(32) ! array of floats from one input line
    integer(Zoltan_INT) :: num_read  ! number of entries in floats_read
    integer :: iostat            ! read status

    nullify(x,y,z)
    call read_real_line(fingeom,floats_read,num_read)

    if (num_read == 0) then
        print *,"No values found in geometry file ", geomname
        close(fingeom)
        chaco_input_geom = .false.
        return
    endif

    xc = floats_read(1)
    ndims = 1
    if (num_read > 1) then
        ndims = 2
        yc = floats_read(2)
        if (num_read > 2) then
            ndims = 3
            zc = floats_read(3)
            if (num_read > 3) then
                print *,"Too many values on input line of geometry file ", &
                       geomname

                print *," Maximum dimensionality is 3"
                close(fingeom)
                chaco_input_geom = .false.
                return
            endif
        endif
    endif

    igeom = ndims

    allocate(x(0:nvtxs))
    x(0) = xc
    if (ndims > 1) then
        allocate(y(0:nvtxs))
        y(0) = yc
    endif
    if (ndims > 2) then
        allocate(z(0:nvtxs))
        z(0) = zc
    endif

    do nread = 1, nvtxs-1
        if (ndims == 1) then
            read(fingeom,*,iostat=iostat) x(nread)
        else if (ndims == 2) then
            read(fingeom,*,iostat=iostat) x(nread),y(nread)
        else if (ndims == 3) then
            read(fingeom,*,iostat=iostat) x(nread),y(nread),z(nread)
        endif

        if (iostat /= 0) then
            print *,"Too few lines or wrong number of values in geometry file ",geomname
            print *,"nvtxs = ",nvtxs,"; have read ",nread+1
            close(fingeom)
            chaco_input_geom = .false.
            return
        endif
    end do

!     Check for spurious extra stuff in file. 
    call read_real_line(fingeom,floats_read,num_read)
    if (num_read > 0 .and. CHECK_INPUT) then
        print *,"Warning: possible error in geometry file ", geomname
        print *," Numerical data found after expected end of file"
    endif

    close(fingeom)

    chaco_input_geom = .true.
end function chaco_input_geom

logical function chaco_dist_graph(comm, host_proc, nvtxs, vtxdist, xadj, &
                                  adjncy, vwgts, ewgts, ndim, x, y, z)

  integer :: comm                        ! MPI Communicator 
  integer(Zoltan_INT) :: host_proc        ! processor where all the data is initially 
  integer(Zoltan_INT) :: nvtxs               ! number of vertices in graph 
  integer(Zoltan_INT), pointer :: vtxdist(:) ! vertex distribution data 
  integer(Zoltan_INT), pointer :: xadj(:)    ! start of edge list for each vertex
  integer(Zoltan_INT), pointer :: adjncy(:)  ! edge list data 
  real(Zoltan_FLOAT), pointer :: vwgts(:)   ! vertex weight list data 
  real(Zoltan_FLOAT), pointer :: ewgts(:)    ! edge weight list data 
  integer(Zoltan_INT) :: ndim                ! dimension of the geometry 
  real(Zoltan_FLOAT), pointer :: x(:)        ! x-coordinates of the vertices 
  real(Zoltan_FLOAT), pointer :: y(:)        ! y-coordinates of the vertices 
  real(Zoltan_FLOAT), pointer :: z(:)        ! z-coordinates of the vertices 

!
! * Distribute a graph from one processor to all processors.
! * The ParMetis format is used for the distributed graph.
! * The memory for the graph on the host node is freed
! * and fresh memory is allocated for the distr. graph.
! 

  integer(Zoltan_INT) :: nprocs, myproc, i, n, p, nedges, nsend, rest
  integer(Zoltan_INT) :: ierr, allocstat
  integer(Zoltan_INT) :: offset, use_vwgts, use_ewgts, use_graph
  integer(Zoltan_INT), pointer, dimension(:) :: old_xadj, old_adjncy,  &
                                            size
  real(Zoltan_FLOAT), pointer, dimension(:) :: old_x, old_y, old_z, old_vwgts, &
                                           old_ewgts
  integer :: status(MPI_STATUS_SIZE)

  nullify(old_xadj, old_adjncy, old_vwgts, size, old_x, old_y, old_z, old_ewgts)

!   Determine number of processors and my rank. 
  call MPI_Comm_size (comm, nprocs, ierr )
  call MPI_Comm_rank (comm, myproc, ierr )

!   Initialize 
  if (associated(ewgts)) then
     use_ewgts = 1
  else
     use_ewgts = 0
  endif
  if (associated(vwgts)) then
     use_vwgts = 1
  else
     use_vwgts = 0
  endif
  if (associated(xadj)) then
     use_graph = 1
  else
     use_graph = 0
  endif
 
!   Broadcast to all procs 
  call MPI_Bcast( use_vwgts, 1, MPI_INTEGER, host_proc, comm, ierr)
  call MPI_Bcast( use_ewgts, 1, MPI_INTEGER, host_proc, comm, ierr)
  call MPI_Bcast( use_graph, 1, MPI_INTEGER, host_proc, comm, ierr)
  call MPI_Bcast( ndim, 1, MPI_INTEGER, host_proc, comm, ierr)
  call MPI_Bcast( nvtxs, 1, MPI_INTEGER, host_proc, comm, ierr)
  
!   Set up vtxdist data 
  if (.not. associated(vtxdist)) then
    allocate(vtxdist(0:nprocs), stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      chaco_dist_graph = .false.
      return
    endif
  endif
!   Calculate uniform vertex distribution 
  vtxdist(0) = 0
  rest = nvtxs
  do i=0, nprocs-1
    n = rest/(nprocs-i)
    vtxdist(i+1) = vtxdist(i) + n
    rest = rest - n
  end do

!   Store pointers to original data 
  if (myproc == host_proc) then
    old_xadj   => xadj
    old_adjncy => adjncy
    old_x      => x
    old_y      => y
    old_z      => z
  endif

!   Allocate space for new distributed graph data 
  n = vtxdist(myproc+1)- vtxdist(myproc) ! local # of nodes 
  nvtxs = n
  if (use_graph /= 0) then
    allocate(xadj(0:n),stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      chaco_dist_graph = .false.
      return
    endif
  endif
  if (use_vwgts /= 0) then
    old_vwgts => vwgts
    allocate(vwgts(0:n-1),stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      chaco_dist_graph = .false.
      return
    endif
  endif
  if (ndim > 0) then
    allocate(x(0:n-1))
    if (ndim > 1) then
      allocate(y(0:n-1))
      if (ndim > 2) then
        allocate(z(0:n-1))
      endif
    endif
  endif


!   Distribute graph data to all procs 
!   Send xadj and coordinates, if appropriate 

  if (myproc == host_proc) then
    if (use_graph /= 0) then
      allocate(size(0:nprocs-1),stat=allocstat)
      if (allocstat /= 0) then
        print *, "fatal: insufficient memory"
        chaco_dist_graph = .false.
        return
      endif
    endif

    do p = 0, nprocs-1
      if (use_graph /= 0) then
        size(p) = old_xadj(vtxdist(p+1))-old_xadj(vtxdist(p))
      endif
      offset = vtxdist(p)
      nsend = vtxdist(p+1) - offset
      if (p == myproc) then
!         do a local copy 
        if (use_graph /= 0) then
          do i=0, nsend
            xadj(i) = old_xadj(offset+i)
          end do
        endif
        do i=0, nsend-1
          if (use_vwgts /= 0) then
            vwgts(i) = old_vwgts(offset+i)
          endif
          if (ndim > 0) then
            x(i) = old_x(offset+i)
            if (ndim > 1) then
              y(i) = old_y(offset+i)
              if (ndim > 2) then
                z(i) = old_z(offset+i)
              endif
            endif
          endif
        end do
      else
        if (use_graph /= 0) then
          call MPI_Send( old_xadj(offset), nsend+1, MPI_INTEGER, p, 1, comm, ierr)
        endif
        if (use_vwgts /= 0) then
          call MPI_Send( old_vwgts(offset), nsend, MPI_REAL, p, 2, comm, ierr)
        endif
        if (ndim > 0) then
          call MPI_Send(old_x(offset), nsend, MPI_REAL, p, 3, comm, ierr)
          if (ndim > 1) then
            call MPI_Send(old_y(offset), nsend, MPI_REAL, p, 4, comm, ierr)
            if (ndim > 2) then
              call MPI_Send(old_z(offset), nsend, MPI_REAL, p, 5, comm, ierr)
            endif
          endif
        endif
      endif
    end do
  else
    if (use_graph /= 0) then
      call MPI_Recv (xadj, nvtxs+1, MPI_INTEGER, host_proc, 1, comm, status, ierr)
    endif
    if (use_vwgts /= 0) then
      call MPI_Recv (vwgts, nvtxs, MPI_REAL, host_proc, 2, comm, status, ierr)
    endif
    if (ndim > 0) then
      call MPI_Recv(x, nvtxs,  MPI_REAL, host_proc, 3, comm, status, ierr)
      if (ndim > 1) then
        call MPI_Recv(y, nvtxs, MPI_REAL, host_proc, 4, comm, status, ierr)
        if (ndim > 2) then
          call MPI_Recv(z, nvtxs, MPI_REAL, host_proc, 5, comm, status, ierr)
        endif
      endif
    endif
  endif

  if (use_graph /= 0) then
!     Adjust xadj values 
    offset = xadj(0)
    do i=0, nvtxs
      xadj(i) = xadj(i) - offset
    end do

!     Distribute adjacency info 
    nedges = xadj(nvtxs)
    if (nedges > 0) then
      allocate(adjncy(0:nedges-1),stat=allocstat)
      if (allocstat /= 0) then
        print *, "fatal: insufficient memory"
        chaco_dist_graph = .false.
        return
      endif
      if (use_ewgts /= 0) then
        old_ewgts => ewgts
        allocate(ewgts(0:nedges-1),stat=allocstat)
        if (allocstat /= 0) then
          print *, "fatal: insufficient memory"
          chaco_dist_graph = .false.
          return
        endif
      endif
    endif

!     Next send adjncy data 
    if (myproc== host_proc) then
      do p=0, nprocs-1
        if (size(p) == 0) cycle
        offset = old_xadj(vtxdist(p))
        if (p == myproc) then
!           do a local copy 
          do i=0, size(p)-1
            adjncy(i) = old_adjncy(offset+i)
          end do
          if (use_ewgts /= 0) then
            do i=0, size(p)-1
              ewgts(i) = old_ewgts(offset+i)
            end do
          endif
        else
          call MPI_Send(old_adjncy(offset), size(p), MPI_INTEGER, p, 6, comm, ierr)
          if (use_ewgts /= 0) then
            call MPI_Send(old_ewgts(offset), size(p), MPI_REAL, p, 7, comm, ierr)
          endif
        endif
      end do
    else
      if (nedges > 0) then
        call MPI_Recv(adjncy, nedges, MPI_INTEGER, host_proc, 6, comm, status, ierr)
        if (use_ewgts /= 0) then
          call MPI_Recv(ewgts, nedges, MPI_REAL, host_proc, 7, comm, status, ierr)
        endif
      endif
    endif
  endif

!   Free space on host proc 
  if (myproc == host_proc) then
    if (associated(old_xadj)) deallocate(old_xadj)
    if (associated(old_adjncy)) deallocate(old_adjncy)
    if (use_vwgts /= 0) deallocate(old_vwgts)
    if (use_ewgts /= 0) deallocate(old_ewgts)

    if (ndim > 0) then
      deallocate(old_x)
      if (ndim > 1) then
        deallocate(old_y)
        if (ndim > 2) then
          deallocate(old_z)
        endif
      endif
    endif
  end if
  if (associated(size)) deallocate(size)
   
  chaco_dist_graph = .true.
end function chaco_dist_graph

!************************************************************

subroutine read_graph_line(fp,ints,num)
integer(Zoltan_INT) :: fp,ints(:),num

! Finds and reads the next noncomment line of a chaco graph file and
! returns all the integers in array ints and the number of ints in num.
! num==0 indicates the end of the file.
! fp is the unit to read from

character(len=256) :: inp_line
integer :: iostat

! read the next noncomment line

read(fp,"(a)",iostat=iostat) inp_line
do while ((index(inp_line,"%") /= 0 .or. index(inp_line,"#") /= 0) &
           .and. iostat==0)
   read(fp,"(a)",iostat=iostat) inp_line
end do

! check for end of file

if (iostat /= 0) then
   num = 0
   return
endif

! read integers from the line

! get rid of leading blanks
inp_line = adjustl(inp_line)
num = 1
do
! read the next integer
   read(inp_line,*,iostat=iostat) ints(num)
   if (iostat /= 0) exit
   num = num + 1
! shift the line to the first nonblank character after the first blank
   inp_line = adjustl(inp_line(index(inp_line," "):))
end do
num = num - 1

end subroutine read_graph_line

!************************************************************

subroutine read_real_line(fp,floats,num)
integer(Zoltan_INT) :: fp,num
real(Zoltan_FLOAT) :: floats(:)

! Finds and reads the next noncomment line of a chaco geometry file and
! returns all the floats in array floats and the number of floats in num.
! num==0 indicates the end of the file.
! fp is the unit to read from

character(len=256) :: inp_line
integer :: iostat

! read the next noncomment line

read(fp,"(a)",iostat=iostat) inp_line
do while (index(inp_line,"%") /= 0 .and. iostat==0)
   read(fp,"(a)",iostat=iostat) inp_line
end do

! check for end of file

if (iostat /= 0) then
   num = 0
   return
endif

! read floats from the line

! get rid of leading blanks
inp_line = adjustl(inp_line)
num = 1
do
! read the next float
   read(inp_line,*,iostat=iostat) floats(num)
   if (iostat /= 0) exit
   num = num + 1
! shift the line to the first nonblank character after the first blank
   inp_line = adjustl(inp_line(index(inp_line," "):))
end do
num = num - 1

end subroutine read_real_line

!************************************************************

subroutine initialize_element(elem)
type(ELEM_INFO) :: elem

!
! * Initializes all fields of an element.
! 
  elem%globalID = -1
  elem%border = 0
  elem%my_part = -1
  elem%elem_blk = -1
  elem%cpu_wgt = 0
  elem%perm_value = -1
  elem%invperm_value = -1
  elem%mem_wgt = 0
  elem%nadj = 0
  elem%adj_len = 0
  nullify(elem%coord, elem%connect, elem%adj, elem%adj_proc, elem%edge_wgt)
end subroutine initialize_element

!************************************************************

subroutine free_element_arrays(elem)
type(ELEM_INFO) :: elem

!
! * Frees all memory malloc'ed for an individual element.
! 

  if (associated(elem%coord)) deallocate(elem%coord)
  if (associated(elem%connect)) deallocate(elem%connect)
  if (associated(elem%adj)) deallocate(elem%adj)
  if (associated(elem%adj_proc)) deallocate(elem%adj_proc)
  if (associated(elem%edge_wgt)) deallocate(elem%edge_wgt)
  elem%globalID = -1
  elem%border = 0
  elem%my_part = -1
  elem%nadj = 0
  elem%adj_len = 0
  elem%elem_blk = -1
  elem%cpu_wgt = 0
  elem%perm_value = -1
  elem%invperm_value = -1
  elem%mem_wgt = 0
end subroutine free_element_arrays

!***************************************************************************
!***************************************************************************
!***************************************************************************

logical function build_elem_comm_maps(proc, elements)
use dr_sort
integer(Zoltan_INT) :: proc
type(ELEM_INFO), target :: elements(0:)

!
! * Build element communication maps, given a distributed mesh.
! * This routine builds initial communication maps for Chaco input
! * (for Nemesis, initial communication maps are read from the Nemesis file)
! * and rebuilds communication maps after data migration.
! *
! * One communication map per neighboring processor is built.
! * The corresponding maps on neighboring processors
! * must be sorted in the same order, so that neighboring processors do not
! * have to use ghost elements.   For each communication map's pair of
! * processors, the lower-numbered processor determines the order of the
! * elements in the communication map.  The sort key is the elements' global
! * IDs on the lower-number processor; the secondary key is the neighboring
! * elements global IDs.  The secondary key is used when a single element
! * must communicate with more than one neighbor.
! 

integer(Zoltan_INT) :: i, j
type(ELEM_INFO), pointer :: elem
integer(Zoltan_INT) :: iadj_elem
integer(Zoltan_INT) :: iadj_proc
integer(Zoltan_INT) :: indx
integer(Zoltan_INT) :: num_alloc_maps
integer(Zoltan_INT) :: max_adj = 0
integer(Zoltan_INT) :: max_adj_per_map
integer(Zoltan_INT) :: cnt, offset
integer(Zoltan_INT), allocatable :: sindex(:)
integer(Zoltan_INT) :: tmp
integer :: astat, astat1, astat2, astat3, astat4
type map_list_head
  integer(Zoltan_INT) :: map_alloc_size
  integer(Zoltan_INT), pointer :: glob_id(:)
  integer(Zoltan_INT), pointer :: elem_id(:)
  integer(Zoltan_INT), pointer :: side_id(:)
  integer(Zoltan_INT), pointer :: neigh_id(:)
end type map_list_head

type(map_list_head), pointer :: tmp_maps(:), map, tmp_map_ptr(:)

!  
!   *  Free the old maps, if they exist.
!   

  if (associated(Mesh%ecmap_id)) then
    deallocate(Mesh%ecmap_id)
    if (associated(Mesh%ecmap_cnt)) deallocate(Mesh%ecmap_cnt)
    if (associated(Mesh%ecmap_elemids)) deallocate(Mesh%ecmap_elemids)
    if (associated(Mesh%ecmap_sideids)) deallocate(Mesh%ecmap_sideids)
    if (associated(Mesh%ecmap_neighids)) deallocate(Mesh%ecmap_neighids)
    Mesh%necmap = 0
  endif

!  
!   *  Look for off-processor adjacencies.
!   *  Loop over all elements 
!   

  num_alloc_maps = MAP_ALLOC
  allocate(Mesh%ecmap_id(0:num_alloc_maps-1), &
           Mesh%ecmap_cnt(0:num_alloc_maps-1), tmp_maps(0:num_alloc_maps-1), &
           stat=astat)

  if (astat /= 0) then
    print *, "Fatal:  insufficient memory"
    build_elem_comm_maps = .false.
    return
  endif

  do i = 0, Mesh%num_elems-1
    elem => elements(i)
    do j = 0, elem%adj_len-1

!      Skip NULL adjacencies (sides that are not adjacent to another elem). 
      if (elem%adj(j) == -1) cycle

      iadj_elem = elem%adj(j)
      iadj_proc = elem%adj_proc(j)

      if (iadj_proc /= proc) then
!         
!         * Adjacent element is off-processor.
!         * Add this element to the temporary data structure for 
!         * the appropriate neighboring processor.
!         
        indx = in_list(iadj_proc, Mesh%necmap, Mesh%ecmap_id)
        if (indx == -1) then
!          
!           * Start a new communication map.
!           

          if (Mesh%necmap >= num_alloc_maps) then
            num_alloc_maps = num_alloc_maps + MAP_ALLOC
            call realloc(Mesh%ecmap_id,num_alloc_maps,astat1)
            call realloc(Mesh%ecmap_cnt,num_alloc_maps,astat2)
            allocate(tmp_map_ptr(0:num_alloc_maps-1),stat=astat)
            tmp_map_ptr(0:size(tmp_maps)-1) = tmp_maps
            deallocate(tmp_maps)
            tmp_maps => tmp_map_ptr
            if (astat /= 0 .or. astat1 /= 0 .or. astat2 /= 0) then
              print *, "Fatal:  insufficient memory"
              build_elem_comm_maps = .false.
              return
            endif
          endif
          Mesh%ecmap_id(Mesh%necmap) = iadj_proc
          Mesh%ecmap_cnt(Mesh%necmap) = 0
          map => tmp_maps(Mesh%necmap)
          allocate(map%glob_id(0:MAP_ALLOC-1),map%elem_id(0:MAP_ALLOC-1), &
             map%side_id(0:MAP_ALLOC-1),map%neigh_id(0:MAP_ALLOC-1), stat=astat)
          if (astat /= 0) then
            print *, "Fatal:  insufficient memory"
            build_elem_comm_maps = .false.
            return
          endif
          map%map_alloc_size = MAP_ALLOC
          indx = Mesh%necmap
          Mesh%necmap = Mesh%necmap + 1
        endif
!         Add to map for indx. 
        map => tmp_maps(indx)
        if (Mesh%ecmap_cnt(indx) >= map%map_alloc_size) then
          map%map_alloc_size = map%map_alloc_size + MAP_ALLOC
          call realloc(map%glob_id,map%map_alloc_size,astat1)
          call realloc(map%elem_id,map%map_alloc_size,astat2)
          call realloc(map%side_id,map%map_alloc_size,astat3)
          call realloc(map%neigh_id,map%map_alloc_size,astat4)
          if (astat1 /= 0 .or. astat2 /= 0 .or. astat3 /= 0 .or. astat4 /= 0) then
            print *, "Fatal:  insufficient memory"
            build_elem_comm_maps = .false.
            return
          endif
        endif       
        tmp = Mesh%ecmap_cnt(indx)
        map%glob_id(tmp) = elem%globalID
        map%elem_id(tmp) = i
        map%side_id(tmp) = j+1  ! side is determined by position in
                                !   adj array (+1 since not 0-based). 
        map%neigh_id(tmp) = iadj_elem
        Mesh%ecmap_cnt(indx) = Mesh%ecmap_cnt(indx) + 1
        max_adj = max_adj + 1
      endif
    end do
  end do

!  
!   * Allocate data structure for element communication map arrays.
!   

  allocate(Mesh%ecmap_elemids(0:max_adj-1), Mesh%ecmap_sideids(0:max_adj-1), &
           Mesh%ecmap_neighids(0:max_adj-1))

!  
!   * Allocate temporary memory for sort index.
!   
  max_adj_per_map = 0
  do i = 0, Mesh%necmap-1
    if (Mesh%ecmap_cnt(i) > max_adj_per_map) then
      max_adj_per_map = Mesh%ecmap_cnt(i)
    endif
  end do
  allocate(sindex(0:max_adj_per_map))

  cnt = 0
  do i = 0, Mesh%necmap-1

    map => tmp_maps(i)
    do j = 0, Mesh%ecmap_cnt(i)-1
      sindex(j) = j
    end do

!    
!     * Sort the map so that adjacent processors have the same ordering
!     * for the communication.  
!     * Assume the ordering of the lower-numbered processor in the pair
!     * of communicating processors.
!     

    if (proc < Mesh%ecmap_id(i)) then
      call dr_sort2_index(0, Mesh%ecmap_cnt(i)-1, map%glob_id, map%neigh_id, sindex)
    else
      call dr_sort2_index(0, Mesh%ecmap_cnt(i)-1, map%neigh_id, map%glob_id, sindex)
    endif

!    
!     * Copy sorted data into elem map arrays. 
!     

    offset = cnt
    do j = 0, Mesh%ecmap_cnt(i)-1
      Mesh%ecmap_elemids(offset)  = map%elem_id(sindex(j))
      Mesh%ecmap_sideids(offset)  = map%side_id(sindex(j))
      Mesh%ecmap_neighids(offset) = map%neigh_id(sindex(j))
      offset = offset + 1
    end do

    cnt = cnt + Mesh%ecmap_cnt(i)
  end do

!   Free temporary data structure. 
  do i = 0, Mesh%necmap-1
    deallocate(tmp_maps(i)%glob_id, tmp_maps(i)%elem_id, &
               tmp_maps(i)%side_id, tmp_maps(i)%neigh_id)
  end do
  deallocate(tmp_maps)
  deallocate(sindex)

  build_elem_comm_maps = .true.
end function build_elem_comm_maps

!***************************************************************************
!***************************************************************************
!***************************************************************************
! Function in_list() begins:
! *----------------------------------------------------------------------------
! * This function searches a vector for the input value. If the value is
! * found in the vector then it's index in that vector is returned, otherwise
! * the function returns -1;
! ****************************************************************************
integer(Zoltan_INT) function in_list(value, count, vector)
integer(Zoltan_INT), intent(in) :: value, count
integer(Zoltan_INT) :: vector(0:)

  integer i


  do i=0, count-1
  
    if(vector(i) == value) then
      in_list = i
      return
    endif
  end do

  in_list = -1
end function in_list

!***************************************************************************
!***************************************************************************
!***************************************************************************

subroutine realloc(array,n,stat)
integer(Zoltan_INT), pointer :: array(:)
integer(Zoltan_INT) :: n, stat

integer(Zoltan_INT), pointer :: tmp(:)
integer(Zoltan_INT) :: lb,ub

lb = lbound(array,dim=1)
ub = ubound(array,dim=1)

allocate(tmp(lb:n+lb-1), stat=stat)
if (stat==0) then
   tmp(lb:ub) = array
   deallocate(array)
   array => tmp
endif

end subroutine realloc

end module dr_chaco_io
