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
use zoltan_user_data
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
integer(Zoltan_INT) :: Proc
type(PROB_INFO) :: prob

!/* Local declarations. */
  type(Zoltan_Struct), pointer :: zz_obj

!  /* Variables returned by the load balancer */
  integer(Zoltan_INT),pointer :: import_gids(:)  !/* Global node nums of nodes to
                                             ! be imported
  integer(Zoltan_INT),pointer :: import_lids(:)  !/* Pointers to nodes to be
                                             ! imported
  integer(Zoltan_INT), pointer :: import_procs(:) !/* Proc IDs of procs owning
                                             ! nodes to be imported.
  integer(Zoltan_INT),pointer :: export_gids(:)  !/* Global node nums of nodes to
                                             ! be exported
  integer(Zoltan_INT),pointer :: export_lids(:)  !/* Pointers to nodes to be
                                             ! exported
  integer(Zoltan_INT), pointer :: export_procs(:) !/* Proc IDs of destination procs
                                             ! for nodes to be exported.
  integer(Zoltan_INT) :: num_imported !/* Number of nodes to be imported.
  integer(Zoltan_INT) :: num_exported !/* Number of nodes to be exported.
  logical :: new_decomp           !/* Flag indicating whether the decomposition
                                  !   has changed

  integer(Zoltan_INT) :: i            !/* Loop index
  integer(Zoltan_INT) :: ierr         !   Return code
  integer(Zoltan_INT) :: num_gid_entries  ! # of array entries in global IDs
  integer(Zoltan_INT) :: num_lid_entries  ! # of array entries in local IDs
  type(Zoltan_User_Data_2) :: mesh_wrapper ! wrapper to pass mesh to query

!/***************************** BEGIN EXECUTION ******************************/

  nullify(zz_obj, import_gids, import_lids, import_procs, &
                  export_gids, export_lids, export_procs)

! make Mesh passable to the callback functions
  mesh_wrapper%ptr => Mesh

!  /*
!   *  Create a load-balancing object.
!   */
  zz_obj => Zoltan_Create(MPI_COMM_WORLD)
  if (.not.associated(zz_obj)) then
    print *, "fatal:  NULL object returned from Zoltan_Create()"
    run_zoltan = .false.
    return
  endif

!  /* Set the user-specified parameters */
  do i = 0, prob%num_params-1
    ierr = Zoltan_Set_Param(zz_obj, trim(prob%params(i)%str(0)), &
                                trim(prob%params(i)%str(1)))
  end do


!  /* Set the method */
  if (Zoltan_Set_Param(zz_obj, "LB_METHOD", prob%method) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Param(LB_METHOD)"
    run_zoltan = .false.
    return
  endif

!  /*
!   * Set the callback functions
!   */

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE, get_num_elements) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Num_Obj_Fn(zz_obj, get_num_elements, &
                        mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_FIRST_OBJ_FN_TYPE, get_first_element, &
!               mesh_wrapper) == ZOLTAN_FATAL) then
  if (Zoltan_Set_First_Obj_Fn(zz_obj, get_first_element, &
                mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NEXT_OBJ_FN_TYPE, get_next_element, &
!               mesh_wrapper) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Next_Obj_Fn(zz_obj, get_next_element, &
                mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE, get_num_geom) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Num_Geom_Fn(zz_obj, get_num_geom, &
                         mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, get_geom, &
!                mesh_wrapper) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Geom_Fn(zz_obj, get_geom, mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Functions for geometry based algorithms */
! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_EDGES_FN_TYPE, get_num_edges, &
!               mesh_wrapper) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Num_Edges_Fn(zz_obj, get_num_edges, &
                mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_EDGE_LIST_FN_TYPE, get_edge_list, &
!               mesh_wrapper) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Edge_List_Fn(zz_obj, get_edge_list, mesh_wrapper) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    return
  endif

!  /* Evaluate the old balance */
!  if (Debug_Driver > 0) {
    if (Proc == 0) then
       print *,"BEFORE load balancing"
    endif
!    driver_eval();
    ierr = Zoltan_LB_Eval(zz_obj, .true.)
!  }

!  /*
!   * Call the load balancer
!   */
  if (Zoltan_LB_Balance(zz_obj, new_decomp, num_gid_entries, num_lid_entries, &
                 num_imported, import_gids,        &
                 import_lids, import_procs, num_exported, export_gids, &
                 export_lids, export_procs) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_LB_Balance()"
    run_zoltan = .false.
    return
  endif

!  /*
!   * Call another routine to perform the migration
!   */
  if (new_decomp) then
    if (.not.migrate_elements(Proc,zz_obj, &
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
    ierr = Zoltan_LB_Eval(zz_obj, .true.)
!  }

!  /* Clean up */
  ierr = Zoltan_LB_Free_Data(import_gids, import_lids, import_procs, &
                      export_gids, export_lids, export_procs)

  call Zoltan_Destroy(zz_obj)

  run_zoltan = .true.

end function run_zoltan

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(Zoltan_INT) function get_num_elements(data, ierr)
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), intent(out) :: ierr

  ierr = ZOLTAN_OK !/* set error code */

  get_num_elements = data%ptr%num_elems
end function get_num_elements

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(Zoltan_INT) function get_first_element(data, &
                                          num_gid_entries, num_lid_entries, &
                                          global_id, local_id, &
                                          wdim, wgt, ierr)

  type(Zoltan_User_Data_2), intent(in) :: data
  integer(Zoltan_INT), intent(in) :: num_gid_entries
  integer(Zoltan_INT), intent(in) :: num_lid_entries
  integer(Zoltan_INT), intent(out) :: global_id(*)
  integer(Zoltan_INT), intent(out) :: local_id(*)
  integer(Zoltan_INT), intent(in) :: wdim
  real(Zoltan_FLOAT), intent(out) :: wgt(*)
  integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(MESH_INFO), pointer :: mesh_data
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (mesh_data%num_elems.eq.0) then  !no elements on this processor
    ierr = ZOLTAN_OK
    get_first_element = 0
    return
  endif
    

  if (.not. associated(elem)) then
    ierr = ZOLTAN_FATAL
    get_first_element = 0
    return
  endif
  
  if (num_lid_entries.gt.0) local_id(lid) = 0
  global_id(gid) = elem(0)%globalID

  if (wdim>0) then
    wgt(1) = elem(0)%cpu_wgt
  endif

  if (wdim>1) then
    ierr = ZOLTAN_WARN ! /* we didn't expect multidimensional weights */
  else
    ierr = ZOLTAN_OK
  endif

  get_first_element = 1
end function get_first_element

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(Zoltan_INT) function get_next_element(data, &
                     num_gid_entries, num_lid_entries, global_id, local_id, &
                     next_global_id, next_local_id, wdim, next_wgt, ierr)
  type(Zoltan_User_Data_2), intent(in) :: data
  integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
  integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*)
  integer(Zoltan_INT), intent(out) :: next_global_id(*), next_local_id(*)
  integer(Zoltan_INT), intent(in) :: wdim
  real(Zoltan_FLOAT), intent(out) :: next_wgt(*)
  integer(Zoltan_INT), intent(out) :: ierr

  integer(Zoltan_INT) :: found
  type(ELEM_INFO), pointer :: elem(:), current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  found = 0
  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = ZOLTAN_FATAL
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
      ierr = ZOLTAN_WARN !/* we didn't expect multidimensional weights */
    else
      ierr = ZOLTAN_OK
    endif
  endif

  get_next_element = found
end function get_next_element

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(Zoltan_INT) function get_num_geom(data, ierr)
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), intent(out) :: ierr

  ierr = ZOLTAN_OK ! /* set error flag */

  get_num_geom = data%ptr%num_dims
end function get_num_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_geom(data, num_gid_entries, num_lid_entries, &
                    global_id, local_id, coor, ierr)
type (Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
real(Zoltan_DOUBLE), intent(out) :: coor(*)
integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(ELEM_INFO), pointer :: current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(Zoltan_INT) :: i, j
  real(Zoltan_DOUBLE) :: tmp
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => elem(local_id(lid))
  else
    current_elem => search_by_global_id(mesh_data, global_id(gid), idx)
  endif

  if (mesh_data%eb_nnodes(current_elem%elem_blk) == 0) then
    !/* No geometry info was read. */
    ierr = ZOLTAN_FATAL
    return
  endif
  
!  /*
!   * calculate the geometry of the element by averaging
!   * the coordinates of the nodes in its connect table
!   */
  do i = 0, mesh_data%num_dims-1
    tmp = 0.0_Zoltan_DOUBLE
    do j = 0, mesh_data%eb_nnodes(current_elem%elem_blk)-1
      tmp = tmp + current_elem%coord(i,j)
    end do

    coor(i+1) = tmp / mesh_data%eb_nnodes(current_elem%elem_blk)
  end do

  ierr = ZOLTAN_OK
end subroutine get_geom

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
integer(Zoltan_INT) function get_num_edges(data, num_gid_entries, num_lid_entries, &
                                       global_id, local_id, ierr)
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
integer(Zoltan_INT), intent(out) :: ierr

type(ELEM_INFO), pointer :: elem(:), current_elem
type(MESH_INFO), pointer :: mesh_data
integer(Zoltan_INT) :: idx
integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
integer(Zoltan_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = ZOLTAN_FATAL
    get_num_edges = 0
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => elem(local_id(lid))
  else
    current_elem => search_by_global_id(mesh_data, global_id(gid), idx)
  endif

  ierr = ZOLTAN_OK;

  get_num_edges = current_elem%nadj
end function get_num_edges

!/*****************************************************************************/
!/*****************************************************************************/
!/*****************************************************************************/
subroutine get_edge_list (data, num_gid_entries, num_lid_entries, &
                          global_id, local_id, nbor_global_id, &
                          nbor_procs, get_ewgts, nbor_ewgts, ierr)
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*)
integer(Zoltan_INT), intent(out) :: nbor_global_id(*)
integer(Zoltan_INT), intent(out) :: nbor_procs(*)
integer(Zoltan_INT), intent(in) :: get_ewgts
real(Zoltan_FLOAT), intent(out) :: nbor_ewgts(*)
integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: elem(:)
  type(ELEM_INFO), pointer :: current_elem
  type(MESH_INFO), pointer :: mesh_data
  integer(Zoltan_INT) :: i, j, proc, local_elem, mpierr
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries;
  lid = num_lid_entries;

  mesh_data => data%ptr
  elem => mesh_data%elements

  if (.not. associated(elem)) then
    ierr = ZOLTAN_FATAL
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

  ierr = ZOLTAN_OK
end subroutine get_edge_list

end module dr_loadbal
