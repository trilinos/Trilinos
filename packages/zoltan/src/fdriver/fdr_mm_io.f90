!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zoltan Library for Parallel Applications                                   !
! For more info, see the README file in the top-level Zoltan directory.      ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CVS File Information :
!     $RCSfile$
!     $Author$
!     $Date$
!     $Revision$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dr_mm_io
use zoltan
use zoltan_user_data
use mpi_h
use dr_const
use dr_input
use dr_chaco_io
implicit none
private

public :: read_mm_file

integer(Zoltan_INT), parameter :: INITIAL_ZERO = 0
integer(Zoltan_INT), parameter :: INITIAL_LINEAR = 1
integer(Zoltan_INT), parameter :: INITIAL_CYCLIC = 2
integer(Zoltan_INT), parameter :: INITIAL_ROW = 3
integer(Zoltan_INT), parameter :: INITIAL_COL = 4

contains

!/****************************************************************************/
!/****************************************************************************/
!/****************************************************************************/

! Function to read MatrixMarket input; for now, reads only standard 
! MatrixMarket, not MatrixMarket+.

logical function read_mm_file(Proc, Num_Proc, prob, pio_info)
integer(Zoltan_INT) :: Proc, Num_Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

!  /* Local declarations. */
  character(len=FILENAME_MAX+8) :: mm_fname
  character(len=10) :: mm_rep
  character(len=7) :: mm_field
  character(len=19) :: mm_symm
  integer :: i, rest, cnt, sum, n, share, pin, pinproc, tmp, mynext
  integer :: prev_edge, pincnt, edgecnt
  integer, pointer :: sizes(:)

! Values read from matrix market
  integer :: mm_nrow, mm_ncol, mm_nnz, mm_max
  integer, pointer :: mm_iidx(:), mm_jidx(:)
  integer, pointer :: mm_ival(:)
  double precision, pointer :: mm_rval(:)
  complex, pointer :: mm_cval(:)

  integer(Zoltan_INT) :: fp, iostat, allocstat, ierr, status
  integer(Zoltan_INT), pointer ::  vtxdist(:) ! vertex distribution data
  integer(Zoltan_INT), pointer ::  pindist(:) ! pin distribution data
  integer(Zoltan_INT), pointer ::  ibuf(:)    ! pin comm buffer
  integer(Zoltan_INT), pointer ::  jbuf(:)    ! pin comm buffer
  integer(Zoltan_INT), pointer ::  next(:)    ! counters
! Local values
  integer(Zoltan_INT) :: npins, nedges, nvtxs
  integer(Zoltan_INT), pointer ::  iidx(:) ! pin data
  integer(Zoltan_INT), pointer ::  jidx(:) ! pin data
!/***************************** BEGIN EXECUTION ******************************/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

! KDDKDD  Switch fields for testing without sort.
! KDDKDD   call mminfo(fp, mm_rep, mm_field, mm_symm, mm_nrow, mm_ncol, mm_nnz)
    call mminfo(fp, mm_rep, mm_field, mm_symm, mm_ncol, mm_nrow, mm_nnz)
! KDDKDD

!   read the matrix in on processor 0.
    allocate(mm_iidx(0:mm_nnz-1), stat=allocstat)
    allocate(mm_jidx(0:mm_nnz-1), stat=allocstat)
    allocate(mm_ival(0:mm_nnz-1), stat=allocstat)
    allocate(mm_rval(0:mm_nnz-1), stat=allocstat)
    allocate(mm_cval(0:mm_nnz-1), stat=allocstat)
!    allocate(mm_iidx(mm_nnz), stat=allocstat)
!    allocate(mm_jidx(mm_nnz), stat=allocstat)
!    allocate(mm_ival(mm_nnz), stat=allocstat)
!    allocate(mm_rval(mm_nnz), stat=allocstat)
!    allocate(mm_cval(mm_nnz), stat=allocstat)
    if (allocstat /= 0) then
      print *, "fatal: insufficient memory"
      read_mm_file = .false.
      return
    endif

    mm_max = mm_nnz
! KDDKDD  Switch fields for testing without sort.
! KDDKDD    call mmread(fp, mm_rep, mm_field, mm_symm, mm_nrow, mm_ncol, mm_nnz, &
! KDDKDD                mm_max, mm_iidx, mm_jidx, mm_ival, mm_rval, mm_cval)
    call mmread(fp, mm_rep, mm_field, mm_symm, mm_ncol, mm_nrow, mm_nnz, &
                mm_max, mm_jidx, mm_iidx, mm_ival, mm_rval, mm_cval)
! KDDKDD

!   Don't need the numerical values.
    if (associated(mm_ival)) deallocate(mm_ival)
    if (associated(mm_rval)) deallocate(mm_rval)
    if (associated(mm_cval)) deallocate(mm_cval)

  endif ! Proc == 0

! BCast pertinent info to all procs.
  call MPI_Bcast(mm_ncol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mm_nrow, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(mm_nnz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

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
    Mesh%elements(i)%globalID = vtxdist(Proc) + i
    Mesh%elements(i)%elem_blk = 0
    Mesh%elements(i)%my_part = Proc
    Mesh%elements(i)%perm_value = -1
    Mesh%elements(i)%invperm_value = -1
    Mesh%elements(i)%cpu_wgt = 1
    Mesh%elements(i)%mem_wgt = 1
  enddo
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculate pin distribution 
!  Send pins to owning processor

  allocate(sizes(0:Num_Proc-1))
  if (Proc == 0) then
    if (pio_info%init_dist_pins == INITIAL_LINEAR) then
      allocate(pindist(0:Num_Proc))
      share = mm_nnz / Num_Proc;
      rest = mm_nnz - (Num_Proc * share);
      sum = 0
      do i = 0, Num_Proc
        pindist(i) = sum
        sum = sum + share
        if (i < rest) sum = sum + 1
      enddo
    else
      print *, "INITIAL_LINEAR IS ONLY DISTRIBUTION SUPPORTED FOR PINS"
      read_mm_file = .false.
      return
    endif

!   Compute number of pins to send to each proc.
    do i = 0, Num_Proc-1
      sizes(i) = 0
    enddo
    do pin = 0, mm_nnz-1
      pinproc = getpinproc(pio_info, pin, mm_jidx(pin), mm_iidx(pin), &
                           Proc, Num_Proc, pindist)
      sizes(pinproc) = sizes(pinproc) + 1
    enddo
  endif

  call MPI_Bcast(sizes, Num_Proc, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

! Allocate arrays to receive pins.
  npins = sizes(Proc)
  allocate(iidx(0:npins-1),jidx(0:npins-1),stat=allocstat)

! Compute prefix sum of sizes
  sum = 0
  do i = 0, Num_Proc-1
    tmp = sizes(i)
    sizes(i) = sum
    sum = sum + tmp
  enddo
  sizes(Num_Proc) = sum

  if (Proc == 0) then
!   Fill communication buffer with pins to be sent.
    allocate(ibuf(0:mm_nnz-1))
    allocate(jbuf(0:mm_nnz-1))
    allocate(next(0:Num_Proc-1))
    do i = 0, Num_Proc-1
      next(i) = sizes(i)
    enddo
    mynext = 0
    do pin = 0, mm_nnz-1
      pinproc = getpinproc(pio_info, pin, mm_jidx(pin), mm_iidx(pin), &
                           Proc, Num_Proc, pindist)
      if (pinproc == 0) then
!       Keep this pin here; local copy.
        iidx(mynext) = mm_iidx(pin)
        jidx(mynext) = mm_jidx(pin)
        mynext = mynext + 1
      else
!       Put this pin in buffer for communication.
        ibuf(next(pinproc)) = mm_iidx(pin)
        jbuf(next(pinproc)) = mm_jidx(pin)
        next(pinproc) = next(pinproc) + 1
      endif
    enddo
    do i = 1, Num_Proc-1
      call MPI_Send(ibuf(i), (sizes(i+1)-sizes(i)), MPI_INTEGER, i, 1, &
                    MPI_COMM_WORLD, ierr)
      call MPI_Send(jbuf(i), (sizes(i+1)-sizes(i)), MPI_INTEGER, i, 2, &
                    MPI_COMM_WORLD, ierr)
    enddo
    if (associated(pindist)) deallocate(pindist)
    if (associated(ibuf)) deallocate(ibuf)
    if (associated(jbuf)) deallocate(jbuf)
    if (associated(next)) deallocate(next)
  else
    call MPI_Recv(iidx, npins, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &
                  status, ierr)
    call MPI_Recv(jidx, npins, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, &
                  status, ierr)
  endif
      
  if (associated(vtxdist)) deallocate(vtxdist)
  if (associated(sizes)) deallocate(sizes)
  if (associated(mm_iidx)) deallocate(mm_iidx)
  if (associated(mm_jidx)) deallocate(mm_jidx)

! Sort the received pins here by edge number.
! KDDKDD TODO
! KDDKDD TODO

! Count number of unique edge IDs.
  prev_edge = -1
  nedges = 0
  do i = 0, npins-1
    if (iidx(i) .ne. prev_edge) nedges = nedges + 1
    prev_edge = iidx(i)
  enddo
  Mesh%nhedges = nedges

! Allocate the index and pin arrays.
  allocate(Mesh%hindex(0:nedges),Mesh%hvertex(0:npins-1),stat=allocstat)

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
  if (associated(iidx)) deallocate(iidx)
  if (associated(jidx)) deallocate(jidx)
  read_mm_file = .true.
end function read_mm_file

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/

integer function getpinproc(pio_info, pin, vid, eid, myrank, nprocs, pindist)
type(PARIO_INFO) :: pio_info
integer pin, vid, eid, myrank, nprocs
integer, pointer ::  pindist(:) 
integer :: i

! Function to return the processor that will own a pin.
! Copied from the C version.  Currently, supporting only INITIAL_LINEAR.

  if (pio_info%init_dist_pins == INITIAL_ZERO) then
!   Node zero initially has all pins.
    getpinproc = 0
  else if (pio_info%init_dist_pins == INITIAL_CYCLIC) then
!   Deal out the pins in a cyclic fashion
    getpinproc = modulo(pin, nprocs)
  else if (pio_info%init_dist_pins == INITIAL_LINEAR) then
!   First process gets first npins/nprocs pins, and so on 
    i = pin / nprocs
    do while (pin < pindist(i))
      i = i - 1
    enddo
    do while (pin >= pindist(i+1))
      i = i + 1
    enddo
    getpinproc = i
  else if (pio_info%init_dist_pins == INITIAL_ROW) then
!   Each process gets entire rows of pins, no row is split across procs 
    getpinproc = modulo(eid, nprocs)
  else if (pio_info%init_dist_pins == INITIAL_COL) then
!   Each process gets entire columns of pins, no column is split across procs 
    getpinproc = modulo(vid, nprocs)
  endif

end function getpinproc

end module dr_mm_io
