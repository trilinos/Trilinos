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

logical function run_zoltan(Proc, prob)
integer(LB_INT) :: Proc
type(PROB_INFO) :: prob

!/* Local declarations. */
  type(LB_Struct), pointer :: lb_obj

!  /* Variables returned by the load balancer */
  integer(LB_INT),pointer :: import_gids(:)  !/* Global node nums of nodes to
                                             ! be imported
  integer(LB_INT),pointer :: import_lids(:)  !/* Pointers to nodes to be
                                             ! imported
  integer(LB_INT), pointer :: import_procs(:) !/* Proc IDs of procs owning
                                             ! nodes to be imported.
  integer(LB_INT),pointer :: export_gids(:)  !/* Global node nums of nodes to
                                             ! be exported
  integer(LB_INT),pointer :: export_lids(:)  !/* Pointers to nodes to be
                                             ! exported
  integer(LB_INT), pointer :: export_procs(:) !/* Proc IDs of destination procs
                                             ! for nodes to be exported.
  integer(LB_INT) :: num_imported !/* Number of nodes to be imported.
  integer(LB_INT) :: num_exported !/* Number of nodes to be exported.
  logical :: new_decomp           !/* Flag indicating whether the decomposition
                                  !   has changed

  integer(LB_INT) :: i            !/* Loop index
  integer(LB_INT) :: ierr         !   Return code
  integer(LB_INT) :: num_gid_entries  ! # of array entries in global IDs
  integer(LB_INT) :: num_lid_entries  ! # of array entries in local IDs
  type(LB_User_Data_2) :: mesh_wrapper ! wrapper to pass mesh to query

!/***************************** BEGIN EXECUTION ******************************/

  nullify(lb_obj, import_gids, import_lids, import_procs, &
                  export_gids, export_lids, export_procs)

! make Mesh passable to the callback functions
  mesh_wrapper%ptr => Mesh

!  /*
!   *  Create a load-balancing object.
!   */
  lb_obj => LB_Create(MPI_COMM_WORLD)
  if (.not.associated(lb_obj)) then
    print *, "fatal:  NULL object returned from LB_Create()"
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

! if (LB_Set_Fn(lb_obj, LB_NUM_OBJ_FN_TYPE, get_num_elements) == LB_FATAL) then
  if (LB_Set_Num_Obj_Fn(lb_obj, get_num_elements, &
                        mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (LB_Set_Fn(lb_obj, LB_FIRST_OBJ_FN_TYPE, get_first_element, &
!               mesh_wrapper) == LB_FATAL) then
  if (LB_Set_First_Obj_Fn(lb_obj, get_first_element, &
                mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (LB_Set_Fn(lb_obj, LB_NEXT_OBJ_FN_TYPE, get_next_element, &
!               mesh_wrapper) == LB_FATAL) then
  if (LB_Set_Next_Obj_Fn(lb_obj, get_next_element, &
                mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
! if (LB_Set_Fn(lb_obj, LB_NUM_GEOM_FN_TYPE, get_num_geom) == LB_FATAL) then
  if (LB_Set_Num_Geom_Fn(lb_obj, get_num_geom, &
                         mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (LB_Set_Fn(lb_obj, LB_GEOM_FN_TYPE, get_geom, &
!                mesh_wrapper) == LB_FATAL) then
  if (LB_Set_Geom_Fn(lb_obj, get_geom, mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
! if (LB_Set_Fn(lb_obj, LB_NUM_EDGES_FN_TYPE, get_num_edges, &
!               mesh_wrapper) == LB_FATAL) then
  if (LB_Set_Num_Edges_Fn(lb_obj, get_num_edges, &
                mesh_wrapper) == LB_FATAL) then
    print *, "fatal:  error returned from LB_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (LB_Set_Fn(lb_obj, LB_EDGE_LIST_FN_TYPE, get_edge_list, &
!               mesh_wrapper) == LB_FATAL) then
  if (LB_Set_Edge_List_Fn(lb_obj, get_edge_list, mesh_wrapper) == LB_FATAL) then
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
    ierr = LB_Eval(lb_obj, .true.)
!  }

!  /*
!   * Call the load balancer
!   */
  if (LB_Balance(lb_obj, new_decomp, num_gid_entries, num_lid_entries, &
                 num_imported, import_gids,        &
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
    if (.not.migrate_elements(Proc,lb_obj, &
                          num_gid_entries, num_lid_entries, &
                          num_imported,import_gids, &
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
    ierr = LB_Eval(lb_obj, .true.)
!  }

!  /* Clean up */
  ierr = LB_Free_Data(import_gids, import_lids, import_procs, &
                      export_gids, export_lids, export_procs)

  call LB_Destroy(lb_obj)

  run_zoltan = .true.

end function run_zoltan

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_num_elements(data, ierr)
type(LB_User_Data_2), intent(in) :: data
integer(LB_INT), intent(out) :: ierr

  ierr = LB_OK !/* set error code */

  get_num_elements = data%ptr%num_elems
end function get_num_elements

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_first_element(data, &
                                          num_gid_entries, num_lid_entries, &
                                          global_id, local_id, &
                                          wdim, wgt, ierr)

  type(LB_User_Data_2), intent(in) :: data
  integer(LB_INT), intent(in) :: num_gid_entries
  integer(LB_INT), intent(in) :: num_lid_entries
  integer(LB_INT), intent(out) :: global_id(*)
  integer(LB_INT), intent(out) :: local_id(*)
  integer(LB_INT), intent(in) :: wdim
  real(LB_FLOAT), intent(out) :: wgt(*)
  integer(LB_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(MESH_INFO), pointer :: mesh_data
  integer(LB_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(LB_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_first_element = 0
    return
  endif
  
  if (num_lid_entries.gt.0) local_id(lid) = 0
  global_id(gid) = elem(0)%globalID

  if (wdim>0) then
    wgt(1) = elem(0)%cpu_wgt
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
integer(LB_INT) function get_next_element(data, &
                     num_gid_entries, num_lid_entries, global_id, local_id, &
                     next_global_id, next_local_id, wdim, next_wgt, ierr)
  type(LB_User_Data_2), intent(in) :: data
  integer(LB_INT), intent(in) :: num_gid_entries, num_lid_entries
  integer(LB_INT), intent(in) :: global_id(*), local_id(*)
  integer(LB_INT), intent(out) :: next_global_id(*), next_local_id(*)
  integer(LB_INT), intent(in) :: wdim
  real(LB_FLOAT), intent(out) :: next_wgt(*)
  integer(LB_INT), intent(out) :: ierr

  integer(LB_INT) :: found
  type(ELEM_INFO), pointer :: elem(:), current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(LB_INT) :: idx
  integer(LB_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(LB_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  found = 0
  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_next_element = 0
    return
  endif
  
  if (num_lid_entries.gt.0) then
    idx = local_id(lid);
    current_elem => elem(idx);
  else 
    !/* testing zero-length local IDs; search by global ID for current elem */
    current_elem => search_by_global_id(mesh, global_id(gid), idx)
  endif

  if (idx+1 < mesh_data%num_elems) then
    found = 1
    if (num_lid_entries.gt.0) next_local_id(lid) = idx + 1
    next_global_id(gid) = elem(idx+1)%globalID

    if (wdim>0) then
      next_wgt(1) = elem(idx+1)%cpu_wgt
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
type(LB_User_Data_2), intent(in) :: data
integer(LB_INT), intent(out) :: ierr

  ierr = LB_OK ! /* set error flag */

  get_num_geom = data%ptr%num_dims
end function get_num_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_geom(data, num_gid_entries, num_lid_entries, &
                    global_id, local_id, coor, ierr)
type (LB_User_Data_2), intent(in) :: data
integer(LB_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(LB_INT), intent(in) :: global_id(*)
integer(LB_INT), intent(in) :: local_id(*)
real(LB_DOUBLE), intent(out) :: coor(*)
integer(LB_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(ELEM_INFO), pointer :: current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(LB_INT) :: i, j
  real(LB_DOUBLE) :: tmp
  integer(LB_INT) :: idx
  integer(LB_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(LB_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => elem(local_id(lid))
  else
    current_elem => search_by_global_id(mesh_data, global_id(gid), idx)
  endif

  if (mesh_data%eb_nnodes(current_elem%elem_blk) == 0) then
    !/* No geometry info was read. */
    ierr = LB_FATAL
    return
  endif
  
!  /*
!   * calculate the geometry of the element by averaging
!   * the coordinates of the nodes in its connect table
!   */
  do i = 0, mesh_data%num_dims-1
    tmp = 0.0_LB_DOUBLE
    do j = 0, mesh_data%eb_nnodes(current_elem%elem_blk)-1
      tmp = tmp + current_elem%coord(i,j)
    end do

    coor(i+1) = tmp / mesh_data%eb_nnodes(current_elem%elem_blk)
  end do

  ierr = LB_OK
end subroutine get_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(LB_INT) function get_num_edges(data, num_gid_entries, num_lid_entries, &
                                       global_id, local_id, ierr)
type(LB_User_Data_2), intent(in) :: data
integer(LB_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(LB_INT), intent(in) :: global_id(*)
integer(LB_INT), intent(in) :: local_id(*)
integer(LB_INT), intent(out) :: ierr

type(ELEM_INFO), pointer :: elem(:), current_elem
type(MESH_INFO), pointer :: mesh_data
integer(LB_INT) :: idx
integer(LB_INT) :: gid  ! Temporary variables to change positioning of IDs.
integer(LB_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    get_num_edges = 0
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => elem(local_id(lid))
  else
    current_elem => search_by_global_id(mesh_data, global_id(gid), idx)
  endif

  ierr = LB_OK;

  get_num_edges = current_elem%nadj
end function get_num_edges

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_edge_list (data, num_gid_entries, num_lid_entries, &
                          global_id, local_id, nbor_global_id, &
                          nbor_procs, get_ewgts, nbor_ewgts, ierr)
type(LB_User_Data_2), intent(in) :: data
integer(LB_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(LB_INT), intent(in) :: global_id(*), local_id(*)
integer(LB_INT), intent(out) :: nbor_global_id(*)
integer(LB_INT), intent(out) :: nbor_procs(*)
integer(LB_INT), intent(in) :: get_ewgts
real(LB_FLOAT), intent(out) :: nbor_ewgts(*)
integer(LB_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(ELEM_INFO), pointer :: current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(LB_INT) :: i, j, proc, local_elem, mpierr
  integer(LB_INT) :: idx
  integer(LB_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(LB_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = LB_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => elem(local_id(lid))
  else
    current_elem => search_by_global_id(mesh_data, global_id(gid), idx)
  endif

!  /* get the processor number */
  call MPI_Comm_rank(MPI_COMM_WORLD, proc, mpierr)

  j = 1
  do i = 0, current_elem%adj_len-1

!    /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
    if (current_elem%adj(i) == -1) cycle

    if (current_elem%adj_proc(i) == proc) then
      local_elem = current_elem%adj(i)
      nbor_global_id(gid+(j-1)*num_gid_entries) = elem(local_elem)%globalID
    else  ! /* adjacent element on another processor */
      nbor_global_id(gid+(j-1)*num_gid_entries) = current_elem%adj(i)
    endif
    nbor_procs(j) = current_elem%adj_proc(i)

    if (get_ewgts /= 0) then
      if (.not. associated(current_elem%edge_wgt)) then
        nbor_ewgts(j) = 1 !/* uniform weights is default */
      else
        nbor_ewgts(j) = current_elem%edge_wgt(i)
      endif
    endif
    j = j+1
  end do

  ierr = LB_OK
end subroutine get_edge_list

end module dr_loadbal
