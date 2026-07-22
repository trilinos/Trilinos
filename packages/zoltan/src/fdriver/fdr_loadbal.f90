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
module dr_loadbal
use zoltan
use zoltan_user_data
use mpi_h
use dr_const
use dr_input
use dr_migrate
use dr_param_file
implicit none
!private

public :: run_zoltan

!--------------------------------------------------------------------------
! Purpose: Call Zoltan to determine a new load balance.                    
!          Contains all of the callback functions that Zoltan needs        
!          for the load balancing.                                         
!--------------------------------------------------------------------------
! Author(s):  Matthew M. St.John (9226)                                    
!   Translated to Fortran by William F. Mitchell
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Revision History:                                                        
!    10 May 1999:       Date of creation.                                  
!       1 September 1999: Fortran translation
!--------------------------------------------------------------------------

contains

!***************************************************************************
!***************************************************************************
!***************************************************************************

logical function run_zoltan(Proc, prob, pio_info)
integer(Zoltan_INT) :: Proc
type(PROB_INFO) :: prob
type(PARIO_INFO) :: pio_info

! Local declarations. 
  type(Zoltan_Struct), pointer :: zz_obj
  type(Zoltan_Struct), pointer :: zz_obj_copy

!   Variables returned by the load balancer 
  integer(Zoltan_INT),pointer :: import_gids(:)  ! Global nums of elements to
                                                 ! be imported
  integer(Zoltan_INT),pointer :: import_lids(:)  ! Pointers to elements to be
                                                 ! imported
  integer(Zoltan_INT),pointer :: import_procs(:) ! Proc IDs of procs owning
                                                 ! elements to be imported.
  integer(Zoltan_INT),pointer :: import_to_part(:)! Partition to which 
                                                 ! elements are to be imported.
  integer(Zoltan_INT),pointer :: export_gids(:)  ! Global nums of elements to
                                                 ! be exported
  integer(Zoltan_INT),pointer :: export_lids(:)  ! Pointers to elements to be
                                                 ! exported
  integer(Zoltan_INT),pointer :: export_procs(:) ! Proc IDs of destination 
                                                 ! for elements to be exported.
  integer(Zoltan_INT),pointer :: export_to_part(:)! Partition to which 
                                                 ! elements are to be exported.
  integer(Zoltan_INT) :: num_imported ! Number of nodes to be imported.
  integer(Zoltan_INT) :: num_exported ! Number of nodes to be exported.
  logical :: new_decomp           ! Flag indicating whether the decomposition
                                  !   has changed

  integer(Zoltan_INT) :: i            ! Loop index
  integer(Zoltan_INT) :: ierr         !   Return code
  integer(Zoltan_INT) :: num_gid_entries  ! # of array entries in global IDs
  integer(Zoltan_INT) :: num_lid_entries  ! # of array entries in local IDs
  character(8) :: s
  real(Zoltan_FLOAT), allocatable :: psize(:)
  integer(Zoltan_INT), allocatable :: partid(:)
  integer(Zoltan_INT), allocatable :: idx(:)
  integer(Zoltan_INT), allocatable :: order(:)
  integer(Zoltan_INT), allocatable :: order_gids(:), order_lids(:)
  integer(Zoltan_INT), allocatable :: color(:)
  integer(Zoltan_INT), allocatable :: gids(:), lids(:)
  integer(Zoltan_INT) :: nprocs
  integer(Zoltan_INT) :: ndim, lid
  character(len=FILENAME_MAX+1) :: fname
  real(Zoltan_DOUBLE) :: xmin, ymin, zmin, xmax, ymax, zmax
  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: curr_idx
!**************************** BEGIN EXECUTION *****************************

  run_zoltan = .true.
  num_gid_entries = 1
  num_lid_entries = 1
  nullify(zz_obj, zz_obj_copy, import_gids, import_lids, import_procs, &
    import_to_part, export_gids, export_lids, export_procs, export_to_part)


!  Allocate space for arrays. 
  call MPI_Comm_size(zoltan_get_global_comm(), nprocs, ierr)
  allocate(psize(nprocs))
  allocate(partid(nprocs))
  allocate(idx(nprocs))

!  
!   *  Create a load-balancing object.
!   
  zz_obj => Zoltan_Create(zoltan_get_global_comm())
  if (.not.associated(zz_obj)) then
    print *, "fatal:  NULL object returned from Zoltan_Create()"
    run_zoltan = .false.
    goto 9999
  endif

!   Set the user-specified parameters 
  do i = 0, prob%num_params-1
    ierr = Zoltan_Set_Param(zz_obj, trim(prob%params(i)%str(0)), &
                                trim(prob%params(i)%str(1)))
  end do


!   Set the method 
  if (Zoltan_Set_Param(zz_obj, "LB_METHOD", prob%method) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Param(LB_METHOD)"
    run_zoltan = .false.
    goto 9999
  endif

!   if there is a paramfile specified, read it
!     note: contents of this file may override the parameters set above 
  if (prob%ztnPrm_file /= "") then  
    call ztnPrm_read_file(zz_obj, prob%ztnPrm_file, &
         zoltan_get_global_comm())
  endif

!  
!   * Set the callback functions
!   

  if (Test_Local_Partitions == 1) then
!    Compute Proc partitions for each processor 
    s(1:1) = achar(Proc/100 + iachar('0'))
    s(2:2) = achar(Proc/10 + iachar('0'))
    s(3:3) = achar(modulo(Proc,10) + iachar('0'))
    s(4:4) = '\n'
    if (Zoltan_Set_Param(zz_obj, "NUM_LOCAL_PARTS", s) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Param()"
      run_zoltan = .false.
      goto 9999
    endif

  else if (Test_Local_Partitions == 2) then
!    Compute Proc partitions for odd-ranked processors let remaining
!    * partitions be in even-ranked processors. 
   if (modulo(Proc,2) == 1) then
      s(1:1) = achar(Proc/100 + iachar('0'))
      s(2:2) = achar(Proc/10 + iachar('0'))
      s(3:3) = achar(modulo(Proc,10) + iachar('0'))
      s(4:4) = '\n'
      if (Zoltan_Set_Param(zz_obj, "NUM_LOCAL_PARTS", s) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Param()"
        run_zoltan = .false.
        goto 9999
      endif
    endif

  else if (Test_Local_Partitions == 3 .or. Test_Local_Partitions == 5) then
!    Variable partition sizes, but one partition per proc 
    partid(1) = Proc
    idx(1) = 0
    psize(1) = Proc   ! Partition size = myproc 
    if (Test_Local_Partitions == 5) psize(1) = psize(1) + 1
!    Set partition sizes using global numbers. 
    ierr = Zoltan_LB_Set_Part_Sizes(zz_obj, 1, 1, partid, idx, psize)
!    Reset partition sizes for upper half of procs. 
    if (Proc >= nprocs/2) then
      psize(1) = 0.5 + modulo(Proc,2)
      if (Test_Local_Partitions == 5) psize(1) = psize(1) + 1
      ierr = Zoltan_LB_Set_Part_Sizes(zz_obj, 1, 1, partid, idx, psize)
    endif

  else if (Test_Local_Partitions == 4) then
!    Variable number of partitions per proc and variable sizes. 
!    Request Proc partitions for each processor, of size 1/Proc.  
    s(1:1) = achar(Proc/100 + iachar('0'))
    s(2:2) = achar(Proc/10 + iachar('0'))
    s(3:3) = achar(modulo(Proc,10) + iachar('0'))
    s(4:4) = '\n'
    if (Zoltan_Set_Param(zz_obj, "NUM_LOCAL_PARTS", s) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Param()"
      run_zoltan = .false.
      goto 9999
    endif
!    Each partition size is inverse to the no. of partitions on a proc. 
    do i = 1, Proc
      partid(i) = i-1                 !   Local partition number 
      idx(i) = 0
      psize(i) = 1.0/Proc
    end do
    ierr = Zoltan_LB_Set_Part_Sizes(zz_obj, 0, Proc, partid, idx, psize)
  endif

!  Free tenmporary arrays for partition sizes. 
  deallocate(psize)
  deallocate(partid)
  deallocate(idx)

! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE, get_num_elements) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Num_Obj_Fn(zz_obj, get_num_elements) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    goto 9999
  endif

  if (Test_Multi_Callbacks .eq. 1)  then
    if (Zoltan_Set_Obj_List_Fn(zz_obj, get_elements) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  else
    if (Zoltan_Set_First_Obj_Fn(zz_obj, get_first_element) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif

    if (Zoltan_Set_Next_Obj_Fn(zz_obj, get_next_element) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  endif

!   Functions for geometry based algorithms 
! if (Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE, get_num_geom) == ZOLTAN_FATAL) then
  if (Zoltan_Set_Num_Geom_Fn(zz_obj, get_num_geom) == ZOLTAN_FATAL) then
    print *, "fatal:  error returned from Zoltan_Set_Fn()"
    run_zoltan = .false.
    goto 9999
  endif

! ZOLTAN_OBJ_SIZE_FN needed for repartitioning.
  if ((Test_Hypergraph_Callbacks .eq. 1) .or. (Test_Graph_Callbacks .eq. 1))  then
    if (Test_Multi_Callbacks.eq.1) then
      if (Zoltan_Set_Obj_Size_Multi_Fn(zz_obj, migrate_elem_size_multi) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false. 
        goto 9999
      endif
    else
      if (Zoltan_Set_Obj_Size_Fn(zz_obj, migrate_elem_size) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
    endif
  endif

  if (Test_Multi_Callbacks.eq.1) then

    if (Zoltan_Set_Geom_Multi_Fn(zz_obj, get_geom_multi) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif

  else
!   if (Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, get_geom) == ZOLTAN_FATAL) then
    if (Zoltan_Set_Geom_Fn(zz_obj, get_geom) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  endif

!   Functions for graph based algorithms 
  if (Test_Graph_Callbacks .eq. 1)  then
    if (Test_Multi_Callbacks .eq. 1)  then
      if (Zoltan_Set_Num_Edges_Multi_Fn(zz_obj, get_num_edges_multi) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
  
      if (Zoltan_Set_Edge_List_Multi_Fn(zz_obj, get_edge_list_multi) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
  
    else
      if (Zoltan_Set_Num_Edges_Fn(zz_obj, get_num_edges) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
  
      if (Zoltan_Set_Edge_List_Fn(zz_obj, get_edge_list) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
    endif
  endif

!   Functions for hypergraph based algorithms 
  if (Test_Hypergraph_Callbacks .eq. 1)  then
    if (Zoltan_Set_Hg_Size_Cs_Fn(zz_obj, get_hg_size_compressed_pins) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  
    if (Zoltan_Set_Hg_Cs_Fn(zz_obj, get_hg_compressed_pins) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  
  
    if (.true.) then 
  !   Register hypergraph edge weight query functions
  
      if (Zoltan_Set_Hg_Size_Edge_Wts_Fn(zz_obj, get_hg_size_edge_weights) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
  
      if (Zoltan_Set_Hg_Edge_Wts_Fn(zz_obj, get_hg_edge_weights) == ZOLTAN_FATAL) then
        print *, "fatal:  error returned from Zoltan_Set_Fn()"
        run_zoltan = .false.
        goto 9999
      endif
    endif
  endif


  if (Test_Multi_Callbacks .eq. 1) then
    if (Zoltan_Set_Part_Multi_Fn(zz_obj, get_part_multi) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  else
    if (Zoltan_Set_Part_Fn(zz_obj, get_part) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_Set_Fn()"
      run_zoltan = .false.
      goto 9999
    endif
  endif

  if (IAND(Driver_Action,1).gt.0) then
    if (Proc == 0) then
      print *,"BEFORE load balancing"
    endif
    ierr = Zoltan_LB_Eval(zz_obj, .true.)

    if (Test_Gen_Files .ne. 0) then
!      Write output files. 
      fname = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".before"
      ierr = Zoltan_Generate_Files(zz_obj, fname, 1, 1, 1, 0);
    endif

!  
!   * Call the load balancer
!   
    if (Zoltan_LB_Partition(zz_obj, new_decomp, &
                   num_gid_entries, num_lid_entries, &
                   num_imported, import_gids, import_lids, &
                   import_procs, import_to_part,  &
                   num_exported, export_gids, export_lids, &
                   export_procs, export_to_part) == ZOLTAN_FATAL) then
      print *, "fatal:  error returned from Zoltan_LB_Partition()"
      run_zoltan = .false.
      goto 9996
    endif

!  
!   * Call another routine to perform the migration
!   
    if (new_decomp) then
      if (.not.migrate_elements(Proc,zz_obj, &
                            num_gid_entries, num_lid_entries, &
                            num_imported,import_gids, &
                            import_lids,import_procs,import_to_part, &
                            num_exported,export_gids, &
                            export_lids,export_procs,export_to_part)) then
        print *, "fatal:  error returned from migrate_elements()"
        run_zoltan = .false.
        goto 9996
      endif
    endif

!    Test the copy function 

    zz_obj_copy => Zoltan_Copy(zz_obj)
    if (.not.associated(zz_obj_copy)) then
      print *, "fatal:  NULL object returned from Zoltan_Copy()"
      run_zoltan = .false.
      goto 9999
    endif

    ierr = Zoltan_Copy_To(zz_obj, zz_obj_copy)
    if (ierr .ne. 0) then
      print *, "fatal:  error in Zoltan_Copy_To()"
      run_zoltan = .false.
      goto 9999
    endif

    call Zoltan_Destroy(zz_obj_copy)
    

!   Evaluate the new balance 
    if (Proc == 0) then
      print *,"AFTER load balancing"
    endif
    ierr = Zoltan_LB_Eval(zz_obj, .true.)

    if (Test_Gen_Files .ne. 0) then
!      Write output files. 
      fname = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".after"
      ierr = Zoltan_Generate_Files(zz_obj, fname, 1, 1, 1, 0);
    endif

    if (Test_Drops.eq.1) then
      call test_drops_rtn(Proc, mesh, pio_info, zz_obj)
    endif 
  
    if ((prob%method .eq. "RCB").or.(prob%method .eq. "rcb")) then
      ierr = Zoltan_RCB_Box(zz_obj,Proc,ndim,xmin,ymin,zmin,xmax,ymax,zmax)
      print *, "DRIVER ", Proc, " DIM: ", ndim, " BOX: (", &
             xmin, ",", ymin, ",", zmin, ") (", xmax, ",", ymax, ",", zmax, ")"
    endif

!   Clean up 
9996 continue
    ierr = Zoltan_LB_Free_Part(import_gids,  import_lids, &
                               import_procs, import_to_part) 
    ierr = Zoltan_LB_Free_Part(export_gids,  export_lids, &
                               export_procs, export_to_part) 
  
  endif  ! End Driver_Action ==> Do balancing

  if (IAND(Driver_Action,2).gt.0) then
!    Do only ordering if this was specified in the driver input file 

    allocate(order(mesh%num_elems*num_gid_entries));
    allocate(order_gids(mesh%num_elems*num_gid_entries));
    do i = 0, mesh%num_elems-1
      order_gids(num_gid_entries*(i+1)) = mesh%elements(i)%globalID
    enddo

    if (Proc .eq. 0) print *, "BEFORE ordering"

    ierr = Zoltan_Order(zz_obj, num_gid_entries, &
        mesh%num_elems, order_gids, &
        order)
    if (ierr .ne. ZOLTAN_OK) then
      print *, "fatal:  error returned from Zoltan_Order()"
      run_zoltan = .false.
      goto 9998
    endif

!    Evaluate the new ordering 
    if (Proc == 0) print *, "AFTER ordering"

!    Copy ordering permutation into mesh structure 
    do i = 1, mesh%num_elems
      current_elem => search_by_global_id(Mesh, order_gids(num_gid_entries*i), &
                                          curr_idx)
      current_elem%perm_value = order(num_gid_entries*i);
!     TODO:  Add Zoltan fn to compute invperm; for now, do what zdrive does.
      current_elem%invperm_value = current_elem%perm_value;
    enddo

9998 continue
!    Free order data 
    deallocate(order);
    deallocate(order_gids);


  endif  ! End Driver_Action ==> Do ordering

  if (IAND(Driver_Action,4).gt.0) then
!    Do only coloring if this was specified in the driver input file 

    allocate(color(mesh%num_elems));
    allocate(gids(mesh%num_elems*num_gid_entries));
    allocate(lids(mesh%num_elems*num_lid_entries));
    do i = 0, mesh%num_elems-1
      gids(num_gid_entries*(i+1)) = mesh%elements(i)%globalID
      lids(num_lid_entries*(i+1)) = i;  ! Temp until Zoltan_Color_Test changes;
    enddo

    ierr = Zoltan_Color(zz_obj, num_gid_entries, &
        mesh%num_elems, gids, &
        color)
    if (ierr .ne. ZOLTAN_OK) then
      print *, "fatal:  error returned from Zoltan_Color()"
      run_zoltan = .false.
      goto 9994
    endif

!    Verify coloring 
    if (Proc == 0) print *, "Verifying coloring result"
    ierr = Zoltan_Color_Test(zz_obj, num_gid_entries, num_lid_entries, &
         mesh%num_elems, gids, lids, &
         color)
    if (ierr .ne. ZOLTAN_OK) then
       print *, "fatal:  error returned from Zoltan_Color_Test()"
       run_zoltan = .false.
       goto 9994
    endif

!    Copy coloring permutation into mesh structure 
    do i = 0, mesh%num_elems-1
      lid = lids(num_lid_entries * (i + 1))
      mesh%elements(lid)%perm_value = color(i+1);
    enddo

9994 continue
!    Free color data 
    deallocate(color);
    deallocate(gids);
    deallocate(lids);


  endif  ! End Driver_Action ==> Do coloring

9999 continue
  call Zoltan_Destroy(zz_obj)


end function run_zoltan

!***************************************************************************
!****** zfdrive query functions below **************************************
!***************************************************************************
integer(Zoltan_INT) function get_num_elements(data, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(out) :: ierr

  ierr = ZOLTAN_OK ! set error code 

  get_num_elements = Mesh%num_elems
end function get_num_elements

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine get_elements(data, num_gid_entries, num_lid_entries, &
                        global_id, local_id, wdim, wgt, ierr)

  INTEGER(Zoltan_INT), intent(in) :: data(*)
  integer(Zoltan_INT), intent(in) :: num_gid_entries
  integer(Zoltan_INT), intent(in) :: num_lid_entries
  integer(Zoltan_INT), intent(out) :: global_id(*)
  integer(Zoltan_INT), intent(out) :: local_id(*)
  integer(Zoltan_INT), intent(in) :: wdim
  real(Zoltan_FLOAT), intent(out) :: wgt(*)
  integer(Zoltan_INT), intent(out) :: ierr

  integer(Zoltan_INT) :: gid  ! Temp variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid
  integer(Zoltan_INT) :: i

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    return
  endif
  
  do i = 0, Mesh%num_elems-1
    if (num_lid_entries.gt.0) local_id(i*num_lid_entries + lid) = i
    global_id(i*num_gid_entries + gid) = Mesh%elements(i)%globalID

    if (wdim>0) then
      wgt(i*wdim+1) = Mesh%elements(i)%cpu_wgt
    endif

    if (wdim>1) then
      ierr = ZOLTAN_WARN !  we didn't expect multidimensional weights 
    else
      ierr = ZOLTAN_OK
    endif
  enddo

end subroutine get_elements

!***************************************************************************
!***************************************************************************
!***************************************************************************
integer(Zoltan_INT) function get_first_element(data, &
                                          num_gid_entries, num_lid_entries, &
                                          global_id, local_id, &
                                          wdim, wgt, ierr)

  INTEGER(Zoltan_INT), intent(in) :: data(*)
  integer(Zoltan_INT), intent(in) :: num_gid_entries
  integer(Zoltan_INT), intent(in) :: num_lid_entries
  integer(Zoltan_INT), intent(out) :: global_id(*)
  integer(Zoltan_INT), intent(out) :: local_id(*)
  integer(Zoltan_INT), intent(in) :: wdim
  real(Zoltan_FLOAT), intent(out) :: wgt(*)
  integer(Zoltan_INT), intent(out) :: ierr

  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (Mesh%num_elems.eq.0) then  !no elements on this processor
    ierr = ZOLTAN_OK
    get_first_element = 0
    return
  endif
    

  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    get_first_element = 0
    return
  endif
  
  if (num_lid_entries.gt.0) local_id(lid) = 0
  global_id(gid) = Mesh%elements(0)%globalID

  if (wdim>0) then
    wgt(1) = Mesh%elements(0)%cpu_wgt
  endif

  if (wdim>1) then
    ierr = ZOLTAN_WARN !  we didn't expect multidimensional weights 
  else
    ierr = ZOLTAN_OK
  endif

  get_first_element = 1
end function get_first_element

!***************************************************************************
!***************************************************************************
!***************************************************************************
integer(Zoltan_INT) function get_next_element(data, &
                     num_gid_entries, num_lid_entries, global_id, local_id, &
                     next_global_id, next_local_id, wdim, next_wgt, ierr)
  INTEGER(Zoltan_INT), intent(in) :: data(*)
  integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
  integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*)
  integer(Zoltan_INT), intent(out) :: next_global_id(*), next_local_id(*)
  integer(Zoltan_INT), intent(in) :: wdim
  real(Zoltan_FLOAT), intent(out) :: next_wgt(*)
  integer(Zoltan_INT), intent(out) :: ierr

  integer(Zoltan_INT) :: found
  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries

  found = 0

  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    get_next_element = 0
    return
  endif
  
  if (num_lid_entries.gt.0) then
    idx = local_id(lid)
    current_elem => Mesh%elements(idx)
  else 
    ! testing zero-length local IDs search by global ID for current elem 
    current_elem => search_by_global_id(mesh, global_id(gid), idx)
  endif

  if (idx+1 < Mesh%num_elems) then
    found = 1
    if (num_lid_entries.gt.0) next_local_id(lid) = idx + 1
    next_global_id(gid) = Mesh%elements(idx+1)%globalID

    if (wdim>0) then
      next_wgt(1) = Mesh%elements(idx+1)%cpu_wgt
    endif

    if (wdim>1) then
      ierr = ZOLTAN_WARN ! we didn't expect multidimensional weights 
    else
      ierr = ZOLTAN_OK
    endif
  endif

  get_next_element = found
end function get_next_element

!***************************************************************************
!***************************************************************************
!***************************************************************************
integer(Zoltan_INT) function get_num_geom(data, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(out) :: ierr

  ierr = ZOLTAN_OK !  set error flag 

  get_num_geom = Mesh%num_dims
end function get_num_geom

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine get_part_multi(data, num_gid_entries, num_lid_entries, &
                    num_obj, global_id, local_id, parts, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries, num_obj
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
integer(Zoltan_INT), intent(out) :: parts(*), ierr

  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: i
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  do i=0,num_obj-1
    if (num_lid_entries.gt.0) then
      current_elem => Mesh%elements(local_id(i*num_lid_entries+lid))
    else
      current_elem => search_by_global_id(Mesh, &
                                          global_id(i*num_gid_entries+gid), idx)
    endif

    parts(i+1) = current_elem%my_part
  end do

  ierr = ZOLTAN_OK

end subroutine get_part_multi

!***************************************************************************
!***************************************************************************
!***************************************************************************
function get_part(data, num_gid_entries, num_lid_entries, &
                    global_id, local_id, ierr)
integer(Zoltan_INT) :: get_part
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => Mesh%elements(local_id(lid))
  else
    current_elem => search_by_global_id(Mesh, global_id(gid), idx)
  endif

  get_part = current_elem%my_part

  ierr = ZOLTAN_OK

end function get_part

!***************************************************************************
!***************************************************************************
!***************************************************************************

subroutine get_geom(data, num_gid_entries, num_lid_entries, &
                    global_id, local_id, coor, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
real(Zoltan_DOUBLE), intent(out) :: coor(*)
integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: i, j
  real(Zoltan_DOUBLE) :: tmp
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => Mesh%elements(local_id(lid))
  else
    current_elem => search_by_global_id(Mesh, global_id(gid), idx)
  endif

  if (Mesh%eb_nnodes(current_elem%elem_blk) == 0) then
    ! No geometry info was read. 
    ierr = ZOLTAN_FATAL
    return
  endif
  
!  
!   * calculate the geometry of the element by averaging
!   * the coordinates of the nodes in its connect table
!   
  do i = 0, Mesh%num_dims-1
    tmp = 0.0_Zoltan_DOUBLE
    do j = 0, Mesh%eb_nnodes(current_elem%elem_blk)-1
      tmp = tmp + current_elem%coord(i,j)
    end do

    coor(i+1) = tmp / Mesh%eb_nnodes(current_elem%elem_blk)
  end do

  ierr = ZOLTAN_OK
end subroutine get_geom

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine get_geom_multi(data, num_gid_entries, num_lid_entries, &
                    num_obj, global_id, local_id, num_dim, coor, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: num_obj, num_dim
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
real(Zoltan_DOUBLE), intent(out) :: coor(*)
integer(Zoltan_INT), intent(out) :: ierr

integer(Zoltan_INT) :: i

! Not the most efficient implementation -- does not take advantage of multi.

  do i = 0, num_obj-1
    if (num_lid_entries .eq. 0) then
      call get_geom(data, num_gid_entries, num_lid_entries, &
                    global_id(i*num_gid_entries + 1),   &
                    local_id, coor(i*num_dim + 1), ierr)
    else
      call get_geom(data, num_gid_entries, num_lid_entries, &
                    global_id(i*num_gid_entries + 1),   &
                    local_id(i*num_lid_entries + 1),    &
                    coor(i*num_dim + 1), ierr)
    endif
    if (ierr.ne.ZOLTAN_OK) exit
  enddo
    
end subroutine get_geom_multi

!***************************************************************************
!***************************************************************************
!***************************************************************************
integer(Zoltan_INT) function get_num_edges(data, num_gid_entries, num_lid_entries, &
                                       global_id, local_id, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*)
integer(Zoltan_INT), intent(in) :: local_id(*)
integer(Zoltan_INT), intent(out) :: ierr

type(ELEM_INFO), pointer :: current_elem
integer(Zoltan_INT) :: idx
integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    get_num_edges = 0
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => Mesh%elements(local_id(lid))
  else
    current_elem => search_by_global_id(Mesh, global_id(gid), idx)
  endif

  ierr = ZOLTAN_OK

  get_num_edges = current_elem%nadj
end function get_num_edges

!***************************************************************************
!***************************************************************************
!***************************************************************************

subroutine get_num_edges_multi (data, num_gid_entries, num_lid_entries, &
  num_obj, global_id, local_id, num_edges, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries, num_obj
integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*)
integer(Zoltan_INT), intent(out) :: num_edges(*), ierr
integer(Zoltan_INT) :: i

! Not the most efficient implementation -- does not take advantage of multi.

  do i = 0, num_obj-1
    if (num_lid_entries .eq. 0) then
      num_edges(i+1) = get_num_edges(data, num_gid_entries, num_lid_entries, &
                                     global_id(i*num_gid_entries + 1),   &
                                     local_id,    &
                                     ierr)
    else
      num_edges(i+1) = get_num_edges(data, num_gid_entries, num_lid_entries, &
                                     global_id(i*num_gid_entries + 1),   &
                                     local_id(i*num_lid_entries + 1),    &
                                     ierr)
    endif
    if (ierr.ne.ZOLTAN_OK) exit
  enddo

end subroutine get_num_edges_multi

!***************************************************************************
!***************************************************************************
!***************************************************************************
subroutine get_edge_list (data, num_gid_entries, num_lid_entries, &
                          global_id, local_id, nbor_global_id, &
                          nbor_procs, get_ewgts, nbor_ewgts, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*)
integer(Zoltan_INT), intent(out) :: nbor_global_id(*)
integer(Zoltan_INT), intent(out) :: nbor_procs(*)
integer(Zoltan_INT), intent(in) :: get_ewgts
real(Zoltan_FLOAT), intent(out) :: nbor_ewgts(*)
integer(Zoltan_INT), intent(out) :: ierr

  type(ELEM_INFO), pointer :: current_elem
  integer(Zoltan_INT) :: i, j, proc, local_elem, mpierr
  integer(Zoltan_INT) :: idx
  integer(Zoltan_INT) :: gid  ! Temporary variables to change positioning of IDs.
  integer(Zoltan_INT) :: lid

  gid = num_gid_entries
  lid = num_lid_entries


  if (.not. associated(Mesh%elements)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  if (num_lid_entries.gt.0) then
    current_elem => Mesh%elements(local_id(lid))
  else
    current_elem => search_by_global_id(Mesh, global_id(gid), idx)
  endif

!   get the processor number 
  call MPI_Comm_rank(zoltan_get_global_comm(), proc, mpierr)

  j = 1
  do i = 0, current_elem%adj_len-1

!     Skip NULL adjacencies (sides that are not adjacent to another elem). 
    if (current_elem%adj(i) == -1) cycle

    if (current_elem%adj_proc(i) == proc) then
      local_elem = current_elem%adj(i)
      nbor_global_id(gid+(j-1)*num_gid_entries) = Mesh%elements(local_elem)%globalID
    else  !  adjacent element on another processor 
      nbor_global_id(gid+(j-1)*num_gid_entries) = current_elem%adj(i)
    endif
    nbor_procs(j) = current_elem%adj_proc(i)

    if (get_ewgts /= 0) then
      if (.not. associated(current_elem%edge_wgt)) then
        nbor_ewgts(j) = 1 ! uniform weights is default 
      else
        nbor_ewgts(j) = current_elem%edge_wgt(i)
      endif
    endif
    j = j+1
  end do

  ierr = ZOLTAN_OK
end subroutine get_edge_list

!***************************************************************************
!***************************************************************************
!***************************************************************************

subroutine get_edge_list_multi(data, num_gid_entries, num_lid_entries, &
  num_obj, global_id, local_id, num_edges, nbor_global_id, nbor_procs, &
  get_ewgts, nbor_ewgts, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, num_lid_entries, num_obj
integer(Zoltan_INT), intent(in) :: global_id(*), local_id(*), num_edges(*)
integer(Zoltan_INT), intent(out) :: nbor_global_id(*)
integer(Zoltan_INT), intent(out) :: nbor_procs(*)
integer(Zoltan_INT), intent(in) :: get_ewgts
real(Zoltan_FLOAT), intent(out) :: nbor_ewgts(*)
integer(Zoltan_INT), intent(out) :: ierr
integer(Zoltan_INT) :: sum, i

! Not the most efficient implementation -- does not take advantage of multi.

  sum = 0;
  do i = 0, num_obj-1
    if (num_lid_entries .eq. 0) then
      call get_edge_list (data, num_gid_entries, num_lid_entries, &
                          global_id(i*num_gid_entries + 1),   &
                          local_id, &
                          nbor_global_id(num_gid_entries * sum + 1), &
                          nbor_procs(sum + 1),  &
                          get_ewgts, nbor_ewgts(sum * get_ewgts + 1), ierr)
    else
      call get_edge_list (data, num_gid_entries, num_lid_entries, &
                          global_id(i*num_gid_entries + 1),   &
                          local_id(i*num_lid_entries + 1),  &
                          nbor_global_id(num_gid_entries * sum + 1), &
                          nbor_procs(sum + 1),  &
                          get_ewgts, nbor_ewgts(sum * get_ewgts + 1), ierr)
    endif
    sum = sum + num_edges(i+1)
    if (ierr.ne.ZOLTAN_OK) exit
  enddo

end subroutine get_edge_list_multi

!***************************************************************************
!******* Hypergraph query functions ****************************************
!***************************************************************************

subroutine get_hg_size_compressed_pins(data, &
  num_lists, num_pins, fmat, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(out) :: num_lists, num_pins, fmat, ierr

  num_lists = Mesh%nhedges
  fmat = 1 ! ZOLTAN_COMPRESSED_EDGE 

  if (.not. associated(Mesh%hindex)) then
    num_pins = 0
  else
    num_pins = Mesh%hindex(Mesh%nhedges)
  endif
  ierr = ZOLTAN_OK

end subroutine get_hg_size_compressed_pins

subroutine get_hg_compressed_pins(data, num_gid_entries, nedges, &
  npins, fmat, edge_GID, edge_ptr, pin_GID, ierr) 
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(in) :: num_gid_entries, nedges, npins, fmat
!integer(Zoltan_INT), intent(out), dimension(*) :: edge_GID, edge_ptr, pin_GID
integer(Zoltan_INT), intent(out) :: edge_GID(0:*), edge_ptr(0:*), pin_GID(0:*)
integer(Zoltan_INT), intent(out) :: ierr

! Local variables
  integer i, k, q


  !print *, "Hello from hypergraph query get_hg_compressed_pins"

  if (.not. associated(Mesh)) then
    ierr = ZOLTAN_FATAL
    return
  endif

  !print *, "Debug: Copying", nedges, "hyperedges"
  q = 0
  do i= 0, nedges-1
    do k= 1, num_gid_entries-1
      edge_GID(q) = 0
      q = q+1
    end do
    edge_ptr(i) = Mesh%hindex(i)
    edge_GID(q) = Mesh%hgid(i)
    !print *, "Debug: hyperedge ", i, q, " : ", Mesh%hgid(i), Mesh%hindex(i)
    q = q+1
  end do

  !print *, "Copying", npins, "pins"
  q = 0
  do i= 0, npins-1
    do k= 1, num_gid_entries-1
      pin_GID(q) = 0
      q = q+1
    end do
    pin_GID(q) = Mesh%hvertex(i)
    q = q+1
  end do

end subroutine get_hg_compressed_pins

subroutine get_hg_size_edge_weights(data, num_edge, ierr)
INTEGER(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), intent(out) :: num_edge, ierr


  num_edge = Mesh%nhedges
  ierr = ZOLTAN_OK
  
end subroutine get_hg_size_edge_weights

subroutine get_hg_edge_weights(data, num_gid_entries, num_lid_entries, &
           num_edges, edge_weight_dim, edge_GID, edge_LID, edge_weight, ierr) 
INTEGER(Zoltan_INT), INTENT(IN) :: data(*) 
INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries, num_lid_entries, num_edges, edge_weight_dim 
INTEGER(Zoltan_INT), INTENT(OUT), DIMENSION(0:*) :: edge_GID 
INTEGER(Zoltan_INT), INTENT(OUT), DIMENSION(0:*) :: edge_LID 
REAL(Zoltan_FLOAT), INTENT(OUT), DIMENSION(0:*) :: edge_weight 
INTEGER(Zoltan_INT), INTENT(OUT) :: ierr 
! Local variables
integer i, k, q

! TEST: return all unit weights just to test functionality.
!       f90 driver supports only plain MatrixMarket (no weights).


  if (edge_weight_dim > 0) then
    q = 0
    do i= 0, num_edges-1
!     First insert padded hyperedge GID
      do k= 1, num_gid_entries-1
        edge_GID(q) = 0
        q = q+1
      end do
      edge_GID(q) = Mesh%hgid(i)
      q = q+1
!     Then insert corresponding weights; always 1 for now
      do k= 0, edge_weight_dim-1
        edge_weight(i*edge_weight_dim +k) = 1.0 
      end do
    end do
  endif

  ierr = ZOLTAN_OK

end subroutine get_hg_edge_weights

!***************************************************************************
!***************************************************************************
!***************************************************************************


subroutine test_drops_rtn(Proc, mesh, pio_info, zz)
integer(Zoltan_INT) :: Proc
type(MESH_INFO), pointer :: mesh
type(PARIO_INFO) :: pio_info
type(Zoltan_Struct), pointer :: zz

integer(Zoltan_INT) :: data(0)
real(Zoltan_DOUBLE) :: xlo(3), xhi(3), x(3)
character(FILENAME_MAX+1) :: par_out_fname, ctemp
type(ELEM_INFO), pointer :: current_elem
integer :: fp
integer(Zoltan_INT) :: ierr
integer(Zoltan_INT) :: i, tmp
integer(Zoltan_INT) :: Num_Proc
integer(Zoltan_INT) :: max_part, gmax_part
integer(Zoltan_INT) :: gid(1), lid(1)
integer(Zoltan_INT) :: test_both  
              ! If true, test both Zoltan_*_Assign and Zoltan_*_PP_Assign. 
              ! If false, test only Zoltan_*_PP_Assign.                    
              ! True if # partitions == # processors.                      

  mesh => Mesh

  ! Find maximum partition number across all processors. 
  call MPI_Comm_size(zoltan_get_global_comm(), Num_Proc, ierr)
  max_part = -1
  gmax_part = -1
  do i = 0, mesh%num_elems-1
    if (mesh%elements(i)%my_part > max_part) then
      max_part = mesh%elements(i)%my_part
    endif
  end do
  call MPI_Allreduce(max_part, gmax_part, 1, MPI_INTEGER, MPI_MAX, &
                     zoltan_get_global_comm(), ierr)
  if ((gmax_part == (Num_Proc-1)) .and. (Test_Local_Partitions == 0)) then
    test_both = 1
  else
    test_both = 0
  endif

  ! generate the parallel filename for this processor 
  ctemp = pio_info%pexo_fname(1:len_trim(pio_info%pexo_fname))//".drops"
  call gen_par_filename(ctemp, par_out_fname, pio_info, Proc, Num_Proc)
  fp = 12
  open(unit=fp,file=par_out_fname,action="write")

  ! Test unit box 
  xlo(1) = 0.0
  xlo(2) = 0.0
  xlo(3) = 0.0
  xhi(1) = 1.0
  xhi(2) = 1.0
  xhi(3) = 1.0
  call test_box_drops(fp, xlo, xhi, zz, Proc, -1, -1, test_both)

  ! Test box based on this processor 
  if (mesh%num_elems > 0) then
    x(1) = 0.
    x(2) = 0.
    x(3) = 0.
    current_elem => mesh%elements(0)
    lid(1) = 0
    gid(1) = current_elem%globalID

    if (mesh%eb_nnodes(current_elem%elem_blk) == 1) then
      x(1) = current_elem%coord(0,0)
      if (mesh%num_dims > 1) x(2) = current_elem%coord(1,0)
      if (mesh%num_dims > 2) x(3) = current_elem%coord(2,0)
    else 
      call get_geom(data, 1, 1, gid, lid, x, ierr)
    endif

    xlo(1) = x(1)
    xlo(2) = x(2)
    xlo(3) = x(3)
    xhi(1) = x(1) + 1.0
    xhi(2) = x(2) + 2.0
    xhi(3) = x(3) + 3.0
    call test_box_drops(fp, xlo, xhi, zz, Proc, Proc, current_elem%my_part, &
                        test_both)
  endif

  ! Test box that (most likely) includes the entire domain. 
  ! All partitions and processors with partitions should be in the output.  
  xlo(1) = -1000000.
  xlo(2) = -1000000.
  xlo(3) = -1000000.
  xhi(1) = 1000000.
  xhi(2) = 1000000.
  xhi(3) = 1000000.
  if (max_part >= 0) then
    tmp = Proc
  else
    ! do not test for proc if proc has no partitions 
    tmp = -1
  endif
  call test_box_drops(fp, xlo, xhi, zz, Proc, tmp, -1, test_both)

  close(fp)

end subroutine test_drops_rtn

!***************************************************************************
subroutine test_point_drops(fp, x, zz, Proc, procs, proccnt, parts, partcnt, &
                            test_both)
integer :: fp 
real(Zoltan_DOUBLE) :: x(*) 
type(Zoltan_Struct), pointer :: zz
integer(Zoltan_INT) :: Proc
integer(Zoltan_INT) :: procs(*)
integer(Zoltan_INT) :: proccnt
integer(Zoltan_INT) :: parts(*)
integer(Zoltan_INT) :: partcnt
integer(Zoltan_INT) :: test_both
 
integer(Zoltan_INT) :: status
integer(Zoltan_INT) :: one_part, one_proc
integer(Zoltan_INT) :: i
integer(Zoltan_INT) :: found

  if (test_both.ne.0) then
    status = Zoltan_LB_Point_Assign(zz, x, one_proc)
    if (status .ne. ZOLTAN_OK)  then
      write(fp,*) "error returned from Zoltan_LB_Point_Assign()"
    else  
      write(fp, &
            fmt='(i2," Zoltan_LB_Point_Assign    (",es13.6,es14.6,es14.6,") on proc",i2)') &
            Proc, x(1), x(2), x(3), one_proc
      found = 0
      do i = 1,proccnt
        if (one_proc .eq. procs(i)) then
          found = 1
          exit
        endif
      end do
      if (found.eq.0) then
        write(fp, &
              fmt='(i2," Error:  processor ",i3,&
                 &" (from Zoltan_LB_Point_Assign) not in proc list from Zoltan_LB_Box_Assign")')&
              Proc, one_proc
      endif
    endif
  else 
    write(fp,*) Proc, "Zoltan_LB_Point_Assign not tested."
  endif

  status = Zoltan_LB_Point_PP_Assign(zz, x, one_proc, one_part)
  if (status .ne. ZOLTAN_OK) then
    write(fp,*) "error returned from Zoltan_LB_Point_PP_Assign()"
  else 
    write(fp, &
          fmt='(i2," Zoltan_LB_Point_PP_Assign (",es13.6,es14.6,es14.6,") on proc",&
               &i2," part ",i2)') &
          Proc, x(1), x(2), x(3), one_proc, one_part

    found = 0
    do i = 1,proccnt
      if (one_proc .eq. procs(i)) then
        found = 1
        exit
      endif
    end do
    if (found.eq.0) then
      write(fp, &
            fmt='(i2," Error:  processor ",i3,&
           &" (from Zoltan_LB_Point_PP_Assign) not in proc list from Zoltan_LB_Box_PP_Assign")')&
            Proc, one_proc
    endif

    if (partcnt .gt. 0) then
      found = 0
      do i=1,partcnt
        if (one_part .eq. parts(i)) then
          found = 1
          exit
        endif
      end do
      if (found.eq.0) then
        write(fp, &
              fmt='(i2," Error:  partition ",i3,&
           &" (from Zoltan_LB_Point_PP_Assign) not in part list from Zoltan_LB_Box_PP_Assign")')&
              Proc, one_part
      endif
    endif
  endif
end subroutine test_point_drops

!***************************************************************************

subroutine test_box_drops(fp, xlo, xhi, zz, Proc, answer_proc, answer_part, &
                          test_both)
integer :: fp 
real(Zoltan_DOUBLE) ::  xlo(*)
real(Zoltan_DOUBLE) ::  xhi(*)
type(Zoltan_Struct), pointer :: zz
integer(Zoltan_INT) :: Proc 
integer(Zoltan_INT) :: answer_proc ! If >= 0, an expected answer for proc. 
integer(Zoltan_INT) :: answer_part ! If >= 0, an expected answer for part. 
integer(Zoltan_INT) :: test_both
 
integer(Zoltan_INT) ::  status, procfound, partfound
integer(Zoltan_INT) ::  proccnt, partcnt
integer(Zoltan_INT) ::  procs(1000), parts(1000)
real(Zoltan_DOUBLE) ::  x(3)
integer(Zoltan_INT) ::  i

  write(fp,*) " "
  write(fp,*) "-------------------------------------------------------"
  if (test_both .eq. 1) then
    status = Zoltan_LB_Box_Assign(zz, xlo(1), xlo(2), xlo(3), &
                                      xhi(1), xhi(2), xhi(3), &
                                      procs, proccnt)
    if (status .ne. ZOLTAN_OK) then
      write(fp,*) "error returned from Zoltan_LB_Box_Assign()"
    else 
      write(fp, fmt='(i2," Zoltan_LB_Box_Assign    LO: (",es13.6,es14.6,es14.6,")")') &
            Proc, xlo(1), xlo(2), xlo(3)
      write(fp, fmt='(i2,"                         HI: (",es13.6,es14.6,es14.6,")")') &
            Proc, xhi(1), xhi(2), xhi(3)
  
      procfound = 0
      write(fp,fmt='("       On ",i3," Procs: ",100i3)') proccnt, (procs(i),i=1,proccnt)
      do i = 1,proccnt
        if (procs(i) .eq. answer_proc) procfound = 1
      end do
      if (answer_proc .ge. 0 .and. procfound.eq.0) then
        write(fp,*) Proc, " Zoltan_LB_Box_Assign error:  ", &
                    "expected proc ", answer_proc, " not in output proc list"
      endif
    endif
  else 
    write(fp,*) Proc, " Zoltan_LB_Box_Assign not tested."
  endif


  status = Zoltan_LB_Box_PP_Assign(zz, xlo(1), xlo(2), xlo(3),  &
                                       xhi(1), xhi(2), xhi(3),  &
                                       procs, proccnt,  &
                                       parts, partcnt)
  if (status .ne. ZOLTAN_OK) then
    write(fp,*) "error returned from Zoltan_LB_Box_PP_Assign()"
  else 
    write(fp, fmt='(i2," Zoltan_LB_Box_PP_Assign LO: (",es13.6,es14.6,es14.6,")")') &
          Proc, xlo(1), xlo(2), xlo(3)
    write(fp, fmt='(i2,"                         HI: (",es13.6,es14.6,es14.6,")")') &
          Proc, xhi(1), xhi(2), xhi(3)

    procfound = 0
    write(fp,fmt='("       On ",i3," Procs: ",100i3)') proccnt, (procs(i),i=1,proccnt)
    do i = 1,proccnt
      if (procs(i) .eq. answer_proc) procfound = 1
    end do

    partfound = 0
    write(fp,fmt='("       In ",i3," Parts: ",100i3)') partcnt, (parts(i),i=1,partcnt)
    do i=1,partcnt
      if (parts(i) .eq. answer_part) partfound = 1
    end do
    if (answer_proc .ge. 0 .and. procfound.eq.0) then
      write(fp,*) Proc, " Zoltan_LB_Box_PP_Assign error:  ", &
                   "expected proc ", answer_proc, " not in output proc list"
    endif
    if (answer_part .ge. 0 .and. partfound.eq.0) then
      write(fp,*) Proc, " Zoltan_LB_Box_PP_Assign error:  ", &
                  "expected part ", answer_part, "not in output part list"
    endif

    ! Test point assign 
    call test_point_drops(fp, xlo, zz, Proc, procs, proccnt, parts, partcnt, &
                          test_both)
    call test_point_drops(fp, xhi, zz, Proc, procs, proccnt, parts, partcnt, &
                          test_both)
    x(1) = 0.5 * (xlo(1) + xhi(1))
    x(2) = 0.5 * (xlo(2) + xhi(2))
    x(3) = 0.5 * (xlo(3) + xhi(3))
    call test_point_drops(fp, x, zz, Proc, procs, proccnt, parts, partcnt, &
                          test_both)
  endif
end subroutine test_box_drops




end module dr_loadbal
