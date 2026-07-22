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

module dr_mm_io
use zoltan
use zoltan_user_data
use mpi_h
use dr_const
use dr_input
use dr_chaco_io
use dr_sort
implicit none
private

public :: read_mm_file

! Pin distribution is assumed to be linear always.

contains

!**************************************************************************
!**************************************************************************
!**************************************************************************

! Function to read MatrixMarket input; for now, reads only standard 
! MatrixMarket, not MatrixMarket+.

logical function read_mm_file(Proc, Num_Proc, prob, pio_info)
integer(Zoltan_INT) :: Proc, Num_Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

!   Local declarations. 
  character(len=FILENAME_MAX+8) :: mm_fname
  character(len=10) :: mm_rep
  character(len=7) :: mm_field
  character(len=19) :: mm_symm
  integer :: i, rest, cnt, sum, n, share, pin, p, itmp, mynext
  integer :: prev_edge, pincnt, edgecnt

! Values read from matrix market
  integer :: mm_nrow, mm_ncol, mm_nnz, mm_max
  integer, pointer :: mm_iidx(:), mm_jidx(:)
  integer, pointer :: mm_ival(:)
  double precision, pointer :: mm_rval(:)
  complex, pointer :: mm_cval(:)

  integer(Zoltan_INT) :: fp, iostat, allocstat, ierr
  integer ::  status(MPI_STATUS_SIZE)
  integer(Zoltan_INT), pointer ::  vtxdist(:) ! vertex distribution data
  integer(Zoltan_INT), pointer ::  pindist(:) ! pin distribution data
  integer :: sendsize
! Local values
  integer(Zoltan_INT) :: npins, nedges, nvtxs
  integer(Zoltan_INT), allocatable ::  iidx(:) ! pin data
  integer(Zoltan_INT), allocatable ::  jidx(:) ! pin data
  integer(Zoltan_INT), allocatable ::  idx(:)  ! temp index 
  integer(Zoltan_INT), allocatable ::  tmp(:)  ! temp values 
  integer :: prev_i, prev_j, temp
  logical :: sorted

!**************************** BEGIN EXECUTION *****************************

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Set appropriate callbacks for this file type.
  Test_Hypergraph_Callbacks = 1
  Test_Graph_Callbacks = 0

! Read the MatrixMarket file.

  if (Proc == 0) then

!   Open and read the MatrixMarket file. 
!   Use the MatrixMarket reader from NIST.
    fp = 12
    mm_fname = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".mtx"
    open(unit=fp,file=mm_fname,action='read',iostat=iostat)
    if (iostat /= 0) then
      print *, "fatal:  Could not open MatrixMarket file ", mm_fname
      read_mm_file = .false.
      return
    endif

    call mminfo(fp, mm_rep, mm_field, mm_symm, mm_nrow, mm_ncol, mm_nnz)

!   read the matrix in on processor 0.
    nullify(mm_ival, mm_cval)
!KDD  Valgrind reports some errors if mm_ival and mm_cval are not allocated,
!KDD  but we don't need them.   32-bit runs fail on the Mac if we don't
!KDD  allocate them.  The error seems to occur in gfortran as it prepares
!KDD  to call mmread.  So we'll allocate them and then deallocate them below.
!KDD  It may be possible to allocate them smaller if needed, but these sizes
!KDD  are OK for our nightly tests.
    allocate(mm_ival(0:mm_nnz-1), stat=allocstat) !KDD
    allocate(mm_cval(0:mm_nnz-1), stat=allocstat) !KDD
    allocate(mm_iidx(0:mm_nnz-1), stat=allocstat)
    allocate(mm_jidx(0:mm_nnz-1), stat=allocstat)
    allocate(mm_rval(0:mm_nnz-1), stat=allocstat)
    allocate(idx(0:mm_nnz-1), stat=allocstat)
    allocate(tmp(0:mm_nnz-1), stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      read_mm_file = .false.
      return
    endif

    mm_max = mm_nnz
    call mmread(fp, mm_rep, mm_field, mm_symm, mm_nrow, mm_ncol, mm_nnz, &
                mm_max, mm_iidx, mm_jidx, mm_ival, mm_rval, mm_cval)

    if (associated(mm_ival)) deallocate(mm_ival) !KDD
    if (associated(mm_cval)) deallocate(mm_cval) !KDD
!   Don't need the numerical values.
    if (associated(mm_rval)) deallocate(mm_rval)

!   Check if pins are sorted by (i,j) values, with row (i) the major index.
!   We could alternatively skip this test and always sort. 
    sorted = .true.
    prev_i = 0
    prev_j = 0
    do i = 0, mm_nnz-1
      if ((mm_iidx(i) < prev_i) .or. ((mm_iidx(i) ==  prev_i) .and. &
          mm_jidx(i) < prev_j)) then
        sorted = .false.
        exit
      endif
      prev_i = mm_iidx(i)
      prev_j = mm_jidx(i)
    enddo

!   If not sorted by (i,j), then sort and permute arrays.
    if (.not. sorted) then
      do i = 0, mm_nnz-1
        idx(i) = i
        ! EBEB For large matrices, the formula below may cause overflow!
        tmp(i) = mm_ncol*mm_iidx(i)+mm_jidx(i) ! Row major, column minor
      enddo
      !print *, 'Before sort (i):', mm_iidx(0), mm_iidx(1), mm_iidx(2)
      !print *, 'Before sort (j):', mm_jidx(0), mm_jidx(1), mm_jidx(2)
      call dr_sort_index(0, mm_nnz-1, tmp, idx) ! TEST
      ! Permute mm_iidx and mm_jidx
      do i = 0, mm_nnz-1
        tmp(i) = mm_iidx(idx(i))
      enddo
      do i = 0, mm_nnz-1
        mm_iidx(i) = tmp(i)
      enddo
      do i = 0, mm_nnz-1
        tmp(i) = mm_jidx(idx(i))
      enddo
      do i = 0, mm_nnz-1
        mm_jidx(i) = tmp(i)
      enddo
      !print *, 'After sort (i):', mm_iidx(0), mm_iidx(1), mm_iidx(2)
      !print *, 'After sort (j):', mm_jidx(0), mm_jidx(1), mm_jidx(2)

    endif

    do i = 0, mm_nnz-1    !  Decrement edge IDs to match C version
      mm_iidx(i) = mm_iidx(i) - 1
    enddo

    deallocate(idx)
    deallocate(tmp)

  endif ! Proc == 0

! BCast pertinent info to all procs.
  call MPI_Bcast(mm_ncol, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(mm_nrow, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)
  call MPI_Bcast(mm_nnz, 1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Assume linear distribution of vertices.
!  Calculate uniform vertex distribution.
  if (.not. associated(vtxdist)) then
    allocate(vtxdist(0:Num_Proc), stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      read_mm_file = .false.
      return
    endif
  endif
  vtxdist(0) = 0
  rest = mm_ncol
  do i=0, Num_Proc-1
    n = rest/(Num_Proc-i)
    vtxdist(i+1) = vtxdist(i) + n
    rest = rest - n
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create elements associated with owned vertices.
! Initialize Mesh structure for MM mesh. 
  nvtxs = vtxdist(Proc+1) - vtxdist(Proc)
  Mesh%num_elems = nvtxs
  Mesh%elem_array_len = Mesh%num_elems + 5
  Mesh%num_dims = 0
  Mesh%num_el_blks = 1

  allocate(Mesh%eb_ids(0:Mesh%num_el_blks-1), &
           Mesh%eb_cnts(0:Mesh%num_el_blks-1), &
           Mesh%eb_nnodes(0:Mesh%num_el_blks-1), &
           Mesh%eb_nattrs(0:Mesh%num_el_blks-1), stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_mm_file = .false.
    return
  endif

  allocate(Mesh%eb_names(0:Mesh%num_el_blks-1),stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_mm_file = .false.
    return
  endif

  Mesh%eb_ids(0) = 1
  Mesh%eb_cnts(0) = nvtxs
! Assume no coordinates for MatrixMarket vertices.
  Mesh%eb_nnodes(0) = 0
  Mesh%eb_nattrs(0) = 0
  Mesh%eb_names(0) = "mm"

! allocate the element structure array.
  allocate(Mesh%elements(0:Mesh%elem_array_len-1), stat=allocstat)
  if (allocstat /= 0) then
    print *, "fatal: insufficient memory"
    read_mm_file = .false.
    return
  endif

! intialize all of the element structs as unused by
! setting the globalID to -1
  do i = 0, Mesh%elem_array_len-1
    call initialize_element(Mesh%elements(i))
  end do

  do i = 0,nvtxs-1
    Mesh%elements(i)%globalID = 1 + vtxdist(Proc) + i
    Mesh%elements(i)%elem_blk = 0
    Mesh%elements(i)%my_part = Proc
    Mesh%elements(i)%perm_value = -1
    Mesh%elements(i)%invperm_value = -1
    Mesh%elements(i)%cpu_wgt = 1
    Mesh%elements(i)%mem_wgt = 1
  enddo
  if (associated(vtxdist)) deallocate(vtxdist)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculate edge and pin distribution 
!  Send pins for those edges to owning processor

  allocate(pindist(0:Num_Proc))
! ONLY INITIAL_LINEAR edge distribution is supported.
  if (Proc == 0) then
!   Assuming pins are sorted by edge number.
    do i = 0, Num_Proc
      pindist(i) = 0
    enddo
    do i = 0, mm_nnz-1
!     Compute the processor to which the edge goes.
      p = int(float(mm_iidx(i) * Num_Proc) / float(mm_nrow));
      pindist(p) = pindist(p)+1
    enddo
!   Compute prefix sum.
    sum = 0
    do i = 0, Num_Proc-1
      itmp = pindist(i)
      pindist(i) = sum
      sum = sum + itmp
    enddo
    pindist(Num_Proc) = sum   
  endif

! Allocate arrays to receive pins.
  call MPI_Bcast(pindist, Num_Proc+1, MPI_INTEGER, 0, zoltan_get_global_comm(), ierr);
  npins = pindist(Proc+1) - pindist(Proc)
  allocate(iidx(0:npins-1),stat=allocstat)
  allocate(jidx(0:npins-1),stat=allocstat)

  if (Proc == 0) then

!   Fill communication buffer with pins to be sent.
!   Assume INITIAL_LINEAR edge distribution.
    do i = 1, Num_Proc-1
      sendsize = pindist(i+1)-pindist(i)
      call MPI_Send(mm_iidx(pindist(i)), sendsize, MPI_INTEGER, &
                    i, 1, zoltan_get_global_comm(), ierr)
      call MPI_Send(mm_jidx(pindist(i)), sendsize, MPI_INTEGER, &
                    i, 2, zoltan_get_global_comm(), ierr)
    enddo
!   Copy Proc zero's pins.
    do i = 0, pindist(1)-1
      iidx(i) = mm_iidx(i)
      jidx(i) = mm_jidx(i)
    enddo
  else
    call MPI_Recv(iidx, npins, MPI_INTEGER, 0, 1, zoltan_get_global_comm(), &
                  status, ierr)
    call MPI_Recv(jidx, npins, MPI_INTEGER, 0, 2, zoltan_get_global_comm(), &
                  status, ierr)
  endif
     
  if (associated(pindist)) deallocate(pindist)
  if (Proc == 0) then
    if (associated(mm_iidx)) deallocate(mm_iidx)
    if (associated(mm_jidx)) deallocate(mm_jidx)
  endif

! KDDKDD We assume the MatrixMarket file is sorted by row numbers.
! KDDKDD This sort was done on a single processor.

! Count number of unique edge IDs on this processor.
  prev_edge = -1
  nedges = 0
  do i = 0, npins-1
    if (iidx(i) .ne. prev_edge) nedges = nedges + 1
    if (iidx(i) < prev_edge) then
!     KDDKDD see note above.
      print *, "Error in MatrixMarket file.  Entries are not sorted by I index."
      read_mm_file = .false.
      return
    endif
    prev_edge = iidx(i)
  enddo
  Mesh%nhedges = nedges

! Allocate the index and pin arrays.
  allocate(Mesh%hgid(0:nedges-1),Mesh%hindex(0:nedges), &
           Mesh%hvertex(0:npins-1),stat=allocstat)

! Fill the index and pin arrays.
  pincnt = 0
  edgecnt = 0
  prev_edge = -1
  do i = 0, npins-1
    if (iidx(i) .ne. prev_edge) then
      Mesh%hindex(edgecnt) = pincnt
      Mesh%hgid(edgecnt) = iidx(i)
      edgecnt = edgecnt + 1
      prev_edge = iidx(i)
    endif
    Mesh%hvertex(pincnt) = jidx(i)
    pincnt = pincnt + 1
  enddo
  Mesh%hindex(nedges) = npins

! Almost done.
!  if (associated(iidx)) deallocate(iidx)
!  if (associated(jidx)) deallocate(jidx)
  deallocate(iidx)
  deallocate(jidx)
  read_mm_file = .true.
end function read_mm_file

end module dr_mm_io
