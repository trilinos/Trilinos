module dr_loadbal
use zoltan
use lb_user_const
use mpi_h
use dr_const
use dr_migrate
implicit none
private

public :: run_zoltan

!/*--------------------------------------------------------------------------*/
!/* Purpose: Call Zoltan to determine a new load balance.                    */
!/*          Contains all of the callback functions that Zoltan needs        */
!/*          for the load balancing.                                         */
!/*--------------------------------------------------------------------------*/
!/* Author(s):  Matthew M. St.John (9226)                                    */
!   Translated to Fortran by William F. Mitchell
!/*--------------------------------------------------------------------------*/
!/*--------------------------------------------------------------------------*/
!/* Revision History:                                                        */
!/*    10 May 1999:       Date of creation.                                  */
!       1 September 1999: Fortran translation
!/*--------------------------------------------------------------------------*/

contains

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/

logical function run_zoltan(Proc, prob, elements)
integer(LB_INT) :: Proc
type(PROB_INFO) :: prob
type(ELEM_INFO), pointer :: elements(:)

!/* Local declarations. */
  type(LB_Struct), pointer :: lb_obj

!  /* Variables returned by the load balancer */
  integer(LB_INT), pointer :: import_gids(:)  !/* Global node nums of nodes to
                                             ! be imported
  integer(LB_INT), pointer :: import_lids(:)  !/* Pointers to nodes to be
                                             ! imported
  integer(LB_INT), pointer :: import_procs(:) !/* Proc IDs of procs owning
                                             ! nodes to be imported.
  integer(LB_INT), pointer :: export_gids(:)  !/* Global node nums of nodes to
                                             ! be exported
  integer(LB_INT), pointer :: export_lids(:)  !/* Pointers to nodes to be
                                             ! exported
  integer(LB_INT), pointer :: export_procs(:) !/* Proc IDs of destination procs
                                             ! for nodes to be exported.
  integer(LB_INT) :: num_imported !/* Number of nodes to be imported.
  integer(LB_INT) :: num_exported !/* Number of nodes to be exported.
  logical :: new_decomp           !/* Flag indicating whether the decomposition
                                  !   has changed

  integer(LB_INT) :: i            !/* Loop index
  integer(LB_INT) :: ierr         !   Return code
  type(LB_User_Data_1) :: elements_wrapper ! wrapper to pass elements to query

!/***************************** BEGIN EXECUTION ******************************/

  nullify(lb_obj, import_gids, import_lids, import_procs, &
                  export_gids, export_lids, export_procs)

! make elements passable to the callback functions

  elements_wrapper%ptr => elements

!  /*
!   *  Create a load-balancing object.
!   */
  lb_obj => LB_Create_Object(MPI_COMM_WORLD)
  if (.not.associated(lb_obj)) then
    print *, "fatal:  NULL object returned from LB_Create_Object()"
    run_zoltan = .false.
    return
  endif


!  /* Set the user-specified parameters */
  do i = 0, prob%num_params-1
    ierr = LB_Set_Param(lb_obj, trim(prob%params(i)%str(0)), &
                                trim(prob%params(i)%str(1)))
  end do


!  /* Set the method */
  if (LB_Set_Method(lb_obj, prob%method) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Method()"
    run_zoltan = .false.
    return
  endif

!  /*
!   * Set the callback functions
!   */

  if (LB_Set_Fn(lb_obj, LB_NUM_OBJ_FN_TYPE, get_num_elements) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

  if (LB_Set_Fn(lb_obj, LB_FIRST_OBJ_FN_TYPE, get_first_element, &
                elements_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

  if (LB_Set_Fn(lb_obj, LB_NEXT_OBJ_FN_TYPE, get_next_element, &
                elements_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
  if (LB_Set_Fn(lb_obj, LB_NUM_GEOM_FN_TYPE, get_num_geom) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

  if (LB_Set_Fn(lb_obj, LB_GEOM_FN_TYPE, get_geom, &
      elements_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
  if (LB_Set_Fn(lb_obj, LB_NUM_EDGES_FN_TYPE, get_num_edges, &
                elements_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

  if (LB_Set_Fn(lb_obj, LB_EDGE_LIST_FN_TYPE, get_edge_list, &
                elements_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Evaluate the old balance */
!  if (Debug_Driver > 0) {
    if (Proc == 0) then
       print *,"BEFORE load balancing"
    endif
!    driver_eval();
    call LB_Eval(lb_obj, .true., ierr = ierr)
!  }

!  /*
!   * Call the load balancer
!   */
  if (LB_Balance(lb_obj, new_decomp, num_imported, import_gids,        &
                 import_lids, import_procs, num_exported, export_gids, &
                 export_lids, export_procs) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Balance()"
    run_zoltan = .false.
    return
  endif

!  /*
!   * Call another routine to perform the migration
!   */
  if (new_decomp) then
    if (.not.migrate_elements(Proc,elements,lb_obj,num_imported,import_gids, &
                          import_lids,import_procs,num_exported,export_gids, &
                          export_lids,export_procs)) then
      print *, "fatal:  error returned from migrate_elements()"
      run_zoltan = .false.
      return
    endif
  endif

!  /* Evaluate the new balance */
!  if (Debug_Driver > 0) {
    if (Proc == 0) then
      print *,"AFTER load balancing"
    endif
!    driver_eval();
    call LB_Eval(lb_obj, .true., ierr = ierr)
!  }

!  /* Clean up */
  ierr = LB_Free_Data(import_gids, import_lids, import_procs, &
                      export_gids, export_lids, export_procs)

  call LB_Destroy_Object(lb_obj)

  run_zoltan = .true.

end function run_zoltan

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_num_elements(data, ierr)
integer(LB_INT) :: data(*)
integer(LB_INT) :: ierr

  ierr = LB_OK !/* set error code */

  get_num_elements = Mesh%num_elems
end function get_num_elements

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_first_element(data, global_id, local_id, &
                                          wdim, wgt, ierr)

  type(LB_User_Data_1) :: data
  integer(LB_INT) :: global_id
  integer(LB_INT) :: local_id
  integer(LB_INT) :: wdim
  real(LB_FLOAT) :: wgt
  integer(LB_INT) :: ierr

  type(ELEM_INFO), pointer :: elem(:)

  elem => data%ptr

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_first_element = 0
    return
  endif
  
  local_id = 0
  global_id = elem(local_id)%globalID

  if (wdim>0) then
    wgt = elem(local_id)%cpu_wgt
  endif

  if (wdim>1) then
    ierr = LB_WARN ! /* we didn't expect multidimensional weights */
  else
    ierr = LB_OK
  endif

  get_first_element = 1
end function get_first_element

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_next_element(data, global_id, local_id, &
                     next_global_id, next_local_id, wdim, next_wgt, ierr)
  type(LB_User_Data_1) :: data
  integer(LB_INT) :: global_id, next_global_id
  integer(LB_INT) :: local_id, next_local_id
  integer(LB_INT) :: wdim
  real(LB_FLOAT)  :: next_wgt
  integer(LB_INT) :: ierr

  integer(LB_INT) :: found
  type(ELEM_INFO), pointer :: elem(:)

  found = 0
  elem => data%ptr

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_next_element = 0
    return
  endif
  
  if (local_id+1 < Mesh%num_elems) then
    found = 1
    next_local_id = local_id + 1
    next_global_id = elem(next_local_id)%globalID

    if (wdim>0) then
      next_wgt = elem(next_local_id)%cpu_wgt
    endif

    if (wdim>1) then
      ierr = LB_WARN !/* we didn't expect multidimensional weights */
    else
      ierr = LB_OK
    endif
  endif

  get_next_element = found
end function get_next_element

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_num_geom(data, ierr)
integer(LB_INT) :: data(*)
integer(LB_INT) :: ierr

  ierr = LB_OK ! /* set error flag */

  get_num_geom = Mesh%num_dims
end function get_num_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_geom(data, global_id, local_id, coor, ierr)
type (LB_User_Data_1) :: data
integer(LB_INT) :: global_id
integer(LB_INT) :: local_id
real(LB_DOUBLE) :: coor(*)
integer(LB_INT) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  integer(LB_INT) :: i, j
  real(LB_DOUBLE) :: tmp

  elem => data%ptr

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    return
  endif

  if (Mesh%eb_nnodes(elem(local_id)%elem_blk) == 0) then
    !/* No geometry info was read. */
    ierr = LB_FATAL
    return
  endif
  
!  /*
!   * calculate the geometry of the element by averaging
!   * the coordinates of the nodes in its connect table
!   */
  do i = 0, Mesh%num_dims-1
    tmp = 0.0_LB_DOUBLE
    do j = 0, Mesh%eb_nnodes(elem(local_id)%elem_blk)-1
      tmp = tmp + elem(local_id)%coord(i,j)
    end do

    coor(i+1) = tmp / Mesh%eb_nnodes(elem(local_id)%elem_blk)
  end do

  ierr = LB_OK
end subroutine get_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_num_edges(data, global_id, local_id, ierr)
type(LB_User_Data_1) :: data
integer(LB_INT) :: global_id
integer(LB_INT) :: local_id
integer(LB_INT) :: ierr

type(ELEM_INFO), pointer :: elem(:)

  elem => data%ptr

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_num_edges = 0
    return
  endif

  ierr = LB_OK;

  get_num_edges = elem(local_id)%nadj
end function get_num_edges

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_edge_list (data, global_id, local_id, nbor_global_id, &
                          nbor_procs, get_ewgts, nbor_ewgts, ierr)
type(LB_User_Data_1) :: data
integer(LB_INT) :: global_id, nbor_global_id(*)
integer(LB_INT) :: local_id
integer(LB_INT) :: nbor_procs(*), get_ewgts, nbor_ewgts(*), ierr

  type(ELEM_INFO), pointer :: elem(:)
  integer(LB_INT) :: i, j, proc, local_elem, mpierr

  elem => data%ptr

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    return
  endif

!  /* get the processor number */
  call MPI_Comm_rank(MPI_COMM_WORLD, proc, mpierr)

  j = 1
  do i = 0, elem(local_id)%adj_len-1

!    /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
    if (elem(local_id)%adj(i) == -1) cycle

    if (elem(local_id)%adj_proc(i) == proc) then
      local_elem = elem(local_id)%adj(i)
      nbor_global_id(j) = elem(local_elem)%globalID
    else  ! /* adjacent element on another processor */
      nbor_global_id(j) = elem(local_id)%adj(i)
    endif
    nbor_procs(j) = elem(local_id)%adj_proc(i)

    if (get_ewgts /= 0) then
      if (.not. associated(elem(local_id)%edge_wgt)) then
        nbor_ewgts(j) = 1 !/* uniform weights is default */
      else
        nbor_ewgts(j) = elem(local_id)%edge_wgt(i)
      endif
    endif
    j = j+1
  end do

  ierr = LB_OK
end subroutine get_edge_list

end module dr_loadbal
