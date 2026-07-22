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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                          //
!! File:      driver.cc                                                     //
!! Project:   Local HSFC Ordering                                           //
!! Author:    Michael Wolf                                                  //
!! Date Started:  11/02/2009                                                //
!!                                                                          //
!! Description:                                                             //
!!              File tests local HSFC ordering for simple test problem      //
!!                                                                          //
!! $Id: driver.cc 11 2009-11-10 00:15:18Z mmwolf $                          //
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module zoltanRCBex
  use mpi_h
  use zoltan

  implicit none

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Mesh data for RCB
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   integer :: numGlobObjs, numLocObjs
   integer(ZOLTAN_INT), dimension(:), allocatable :: GIDs
   real, dimension(:), allocatable :: xcoords, ycoords
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Zoltan data to store in module
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   LOGICAL :: changes 
   INTEGER(Zoltan_INT) :: numGidEntries, numLidEntries
   INTEGER(Zoltan_INT) :: numImport, numExport
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importGlobalGids, exportGlobalGids 
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importLocalGids, exportLocalGids
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importProcs, exportProcs
   INTEGER(Zoltan_INT), POINTER, DIMENSION(:) :: importToPart, exportToPart
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine partitionMeshWithRCB()
    use zoltan
    use mpi_h

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!
    !! local variables
    !!!!!!!!!!!!!!!!!!!!!!!
    type(Zoltan_Struct), pointer :: zz_obj
    integer(ZOLTAN_INT) :: ierr
    !!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Body of subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nullify(zz_obj)

    zz_obj => Zoltan_Create(MPI_COMM_WORLD)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! General Zoltan Parameters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ierr = Zoltan_Set_Param(zz_obj, "LB_METHOD", "RCB")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! register query functions
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_OBJ_FN_TYPE,zoltNumObjs)
    ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_OBJ_LIST_FN_TYPE,zoltGetObjs)
    ierr = Zoltan_Set_Fn(zz_obj, ZOLTAN_NUM_GEOM_FN_TYPE,zoltNumGeom)
    ierr =  Zoltan_Set_Fn(zz_obj, ZOLTAN_GEOM_FN_TYPE, zoltGeom)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Use Zoltan to partition the vertices in the simple mesh.
    !!
    !! Params:
    !!     zz_obj           -- input (all remaining fields are output)
    !!     changes          -- 1 if partition was changed, 0 otherwise 
    !!     numGidEntries    -- Number of integers used for a global ID 
    !!     numLidEntries    -- Number of integers used for a local ID 
    !!     numImport        -- Number of vertices to be sent to me 
    !!     importGlobalGids -- Global IDs of vertices to be sent to me 
    !!     importLocalGids  -- Local IDs of vertices to be sent to me 
    !!     importProcs      -- Process rank for source of each incoming vertex 
    !!     importToPart     -- New part for each incoming vertex 
    !!     numExport        -- Number of vertices I must send to other processes
    !!     exportGlobalGids -- Global IDs of the vertices I must send 
    !!     exportLocalGids  -- Local IDs of the vertices I must send 
    !!     exportProcs      -- Process to which I send each of the vertices 
    !!     exportToPart     -- Part to which each vertex will belong 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ierr = Zoltan_LB_Partition(zz_obj, changes, numGidEntries, numLidEntries, &
                               numImport, importGlobalGids, importLocalGids, importProcs, importToPart, &
                               numExport, exportGlobalGids, exportLocalGids, exportProcs, exportToPart)


    call Zoltan_Destroy(zz_obj)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine partitionMeshWithRCB
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Frees arrays allocated by Zoltan_LB_Partition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltanCleanup()
    use zoltan
    implicit none

    integer :: error

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    error = Zoltan_LB_Free_Part(importGlobalGids, importLocalGids, importProcs, importToPart)
    error = Zoltan_LB_Free_Part(exportGlobalGids, exportLocalGids, exportProcs, exportToPart)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine zoltanCleanup
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumObjs(data, ierr)
    use zoltan 
    implicit none

    ! Local declarations
    INTEGER(Zoltan_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zoltNumObjs = numLocObjs
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end function zoltNumObjs
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGetObjs (data, num_gid_entries, num_lid_entries, global_ids, & 
                        local_ids, wgt_dim, obj_wgts, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
    integer(ZOLTAN_INT), intent(out) :: global_ids(*)
    integer(ZOLTAN_INT), intent(out) :: local_ids(*)
    integer(ZOLTAN_INT), intent(in) :: wgt_dim 
    real(ZOLTAN_FLOAT), intent(out) :: obj_wgts(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    ! local declarations
    integer :: i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i= 1, numLocObjs
      global_ids(i) = GIDs(i)
      local_ids(i) = i
    end do
    
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer function zoltNumGeom(data, ierr)
    use zoltan 
    implicit none
    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    zoltNumGeom = 2
    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end function zoltNumGeom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! User defined query function to register with Zoltan
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine zoltGeom(data, num_gid_entries, num_lid_entries, global_id, &
                    local_id, geom_vec, ierr)
    use zoltan
    implicit none

    integer(ZOLTAN_INT), intent(in) :: data(*)
    integer(ZOLTAN_INT), intent(in) :: num_gid_entries 
    integer(ZOLTAN_INT), intent(in) ::  num_lid_entries
    integer(ZOLTAN_INT), intent(in) :: global_id
    integer(ZOLTAN_INT), intent(in) :: local_id
    real(ZOLTAN_DOUBLE), intent(out) :: geom_vec(*)
    integer(ZOLTAN_INT), intent(out) :: ierr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    geom_vec(1) =  xcoords(local_id)
    geom_vec(2) =  ycoords(local_id)

    ierr = ZOLTAN_OK
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  end subroutine zoltGeom
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zoltanRCBex
