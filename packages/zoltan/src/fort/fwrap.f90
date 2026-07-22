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

!--------------------------------------------------------------------------
! preprocessor directives to handle special case compilers

module zoltan
use zoltan_types
use zoltan_user_data
implicit none
private

!--------------------------------------------------------------------------
! public entities

public :: &
   Zoltan_INT, &
   Zoltan_FLOAT, &
   Zoltan_DOUBLE, &
   Zoltan_User_Data_1, &
   Zoltan_User_Data_2, &
   Zoltan_User_Data_3, &
   Zoltan_User_Data_4

public :: &
   Zoltan_Struct, &
   ZOLTAN_FN_TYPEF, &
   ZOLTAN_FN_TYPES

public :: &
   ZOLTAN_PART_FN_TYPE, &
   ZOLTAN_PART_MULTI_FN_TYPE, &
   ZOLTAN_NUM_EDGES_FN_TYPE, &
   ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, &
   ZOLTAN_EDGE_LIST_FN_TYPE, &
   ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, &
   ZOLTAN_NUM_GEOM_FN_TYPE, &
   ZOLTAN_GEOM_MULTI_FN_TYPE, &
   ZOLTAN_GEOM_FN_TYPE, &
   ZOLTAN_NUM_OBJ_FN_TYPE, &
   ZOLTAN_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_OBJ_FN_TYPE, &
   ZOLTAN_NUM_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_BORDER_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_PRE_MIGRATE_PP_FN_TYPE, &
   ZOLTAN_MID_MIGRATE_PP_FN_TYPE, &
   ZOLTAN_POST_MIGRATE_PP_FN_TYPE, &
   ZOLTAN_PRE_MIGRATE_FN_TYPE, &
   ZOLTAN_MID_MIGRATE_FN_TYPE, &
   ZOLTAN_POST_MIGRATE_FN_TYPE, &
   ZOLTAN_OBJ_SIZE_FN_TYPE, &
   ZOLTAN_PACK_OBJ_FN_TYPE, &
   ZOLTAN_UNPACK_OBJ_FN_TYPE, &
   ZOLTAN_HIER_NUM_LEVELS_FN_TYPE, &
   ZOLTAN_HIER_PART_FN_TYPE

! Backward compatibility with v3.0
public:: &
   ZOLTAN_PARTITION_FN_TYPE, &    
   ZOLTAN_PARTITION_MULTI_FN_TYPE, &
   ZOLTAN_HIER_PARTITION_FN_TYPE

public:: &
   ZOLTAN_NUM_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_COARSE_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_NUM_CHILD_FN_TYPE, &
   ZOLTAN_CHILD_LIST_FN_TYPE, &
   ZOLTAN_CHILD_WEIGHT_FN_TYPE, &
   ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, &
   ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, &
   ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE, &
   ZOLTAN_HG_SIZE_CS_FN_TYPE, &
   ZOLTAN_HG_CS_FN_TYPE, &
   ZOLTAN_HG_SIZE_EDGE_WTS_FN_TYPE, &
   ZOLTAN_HG_EDGE_WTS_FN_TYPE, &
   ZOLTAN_NUM_FIXED_OBJ_FN_TYPE, &
   ZOLTAN_FIXED_OBJ_LIST_FN_TYPE, &
   ZOLTAN_HIER_METHOD_FN_TYPE

public :: &
   ZOLTAN_OTHER_REF, &
   ZOLTAN_IN_ORDER, &
   ZOLTAN_TRI_BISECT, &
   ZOLTAN_QUAD_QUAD, &
   ZOLTAN_HEX3D_OCT

public :: &
   ZOLTAN_OK, &
   ZOLTAN_WARN, &
   ZOLTAN_FATAL, &
   ZOLTAN_MEMERR

public :: &
   ZOLTAN_COMPRESSED_EDGE, &
   ZOLTAN_COMPRESSED_VERTEX

public :: &
   Zoltan_Initialize, &
   Zoltan_Create, &
   Zoltan_Copy, &
   Zoltan_Copy_To, &
   Zoltan_Destroy, &
   Zoltan_Get_Struct_Addr, &
   Zoltan_Align, &
   Zoltan_Memory_Stats, &
   Zoltan_Set_Fn, &
   Zoltan_Set_Param, &
   Zoltan_Set_Param_Vec, &
   Zoltan_LB_Partition, &
   Zoltan_LB_Eval, &
   Zoltan_LB_Free_Part, &
   Zoltan_LB_Free_Data, &
   Zoltan_LB_Set_Part_Sizes, &
   Zoltan_LB_Point_Assign, &
   Zoltan_LB_Point_PP_Assign, &
   Zoltan_LB_Box_Assign, &
   Zoltan_LB_Box_PP_Assign, &
   Zoltan_LB_Balance, &
   Zoltan_Invert_Lists, &
   Zoltan_Compute_Destinations, &
   Zoltan_Migrate, &
   Zoltan_Help_Migrate, &
   Zoltan_Order, &
   Zoltan_Color, &
   Zoltan_Color_Test, &
   Zoltan_Generate_Files, &
   Zoltan_RCB_Box

! Registration functions with strict type checking.
public :: &
   Zoltan_Set_Num_Obj_Fn, Zoltan_Set_Obj_List_Fn, &
   Zoltan_Set_First_Obj_Fn, Zoltan_Set_Next_Obj_Fn, &
   Zoltan_Set_Num_Border_Obj_Fn, Zoltan_Set_Border_Obj_List_Fn, &
   Zoltan_Set_First_Border_Obj_Fn, Zoltan_Set_Next_Border_Obj_Fn, &
   Zoltan_Set_Num_Geom_Fn, Zoltan_Set_Geom_Multi_Fn, Zoltan_Set_Geom_Fn, &
   Zoltan_Set_Part_Fn, Zoltan_Set_Part_Multi_Fn, &
   Zoltan_Set_Num_Edges_Fn, Zoltan_Set_Num_Edges_Multi_Fn, &
   Zoltan_Set_Edge_List_Fn, Zoltan_Set_Edge_List_Multi_Fn, &
   Zoltan_Set_Num_Coarse_Obj_Fn, Zoltan_Set_Coarse_Obj_List_Fn, &
   Zoltan_Set_First_Coarse_Obj_Fn, Zoltan_Set_Next_Coarse_Obj_Fn, &
   Zoltan_Set_Num_Child_Fn, Zoltan_Set_Child_List_Fn, &
   Zoltan_Set_Child_Weight_Fn, &
   Zoltan_Set_Obj_Size_Fn, Zoltan_Set_Pack_Obj_Fn, Zoltan_Set_Unpack_Obj_Fn, &
   Zoltan_Set_Pre_Migrate_PP_Fn, Zoltan_Set_Mid_Migrate_PP_Fn, &
   Zoltan_Set_Post_Migrate_PP_Fn, &
   Zoltan_Set_Pre_Migrate_Fn, Zoltan_Set_Mid_Migrate_Fn, &
   Zoltan_Set_Post_Migrate_Fn, &
   Zoltan_Set_Obj_Size_Multi_Fn, &
   Zoltan_Set_Pack_Obj_Multi_Fn, Zoltan_Set_Unpack_Obj_Multi_Fn, &
   Zoltan_Set_HG_Size_CS_Fn, Zoltan_Set_HG_CS_Fn, &
   Zoltan_Set_HG_Size_Edge_Wts_Fn, Zoltan_Set_HG_Edge_Wts_Fn,  &
   Zoltan_Set_Num_Fixed_Obj_Fn, Zoltan_Set_Fixed_Obj_List_Fn, &
   Zoltan_Set_Hier_Num_Levels_Fn, Zoltan_Set_Hier_Part_Fn, &
   Zoltan_Set_Hier_Method_Fn

! Backward compatibility with v3.0
public :: &
   Zoltan_Set_Partition_Fn, Zoltan_Set_Partition_Multi_Fn, &
   Zoltan_Set_Hier_Partition_Fn

public :: &
   Zoltan_Get_Child_Order

!--------------------------------------------------------------------------
! defined constants corresponding to Zoltan enumerated types

! Enumerated type used to indicate which function is to be set by Zoltan_Set_Fn.
! These values must agree with those in the Zoltan_Set_Fn wrapper in cwrap.c

type ZOLTAN_FN_TYPEF
   private
   integer(Zoltan_INT) :: choice
end type ZOLTAN_FN_TYPEF

type ZOLTAN_FN_TYPES
   private
   integer(Zoltan_INT) :: choice
end type ZOLTAN_FN_TYPES

type(ZOLTAN_FN_TYPEF), parameter :: &
   ZOLTAN_NUM_EDGES_FN_TYPE        = ZOLTAN_FN_TYPEF(0_Zoltan_INT), &
   ZOLTAN_NUM_GEOM_FN_TYPE         = ZOLTAN_FN_TYPEF(4_Zoltan_INT), &
   ZOLTAN_NUM_OBJ_FN_TYPE          = ZOLTAN_FN_TYPEF(7_Zoltan_INT), &
   ZOLTAN_FIRST_OBJ_FN_TYPE        = ZOLTAN_FN_TYPEF(9_Zoltan_INT), &
   ZOLTAN_NEXT_OBJ_FN_TYPE         = ZOLTAN_FN_TYPEF(10_Zoltan_INT), &
   ZOLTAN_NUM_BORDER_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(11_Zoltan_INT), &
   ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(13_Zoltan_INT), &
   ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(14_Zoltan_INT), &
   ZOLTAN_OBJ_SIZE_FN_TYPE         = ZOLTAN_FN_TYPEF(21_Zoltan_INT), &
   ZOLTAN_NUM_COARSE_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(24_Zoltan_INT), &
   ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(26_Zoltan_INT), &
   ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(27_Zoltan_INT), &
   ZOLTAN_NUM_CHILD_FN_TYPE        = ZOLTAN_FN_TYPEF(28_Zoltan_INT), &
   ZOLTAN_PART_FN_TYPE             = ZOLTAN_FN_TYPEF(34_Zoltan_INT), &
   ZOLTAN_HIER_NUM_LEVELS_FN_TYPE  = ZOLTAN_FN_TYPEF(43_Zoltan_INT), &
   ZOLTAN_HIER_PART_FN_TYPE        = ZOLTAN_FN_TYPEF(44_Zoltan_INT)

type(ZOLTAN_FN_TYPES), parameter :: &
   ZOLTAN_NUM_EDGES_MULTI_FN_TYPE  = ZOLTAN_FN_TYPES(1_Zoltan_INT), &
   ZOLTAN_EDGE_LIST_FN_TYPE        = ZOLTAN_FN_TYPES(2_Zoltan_INT), &
   ZOLTAN_EDGE_LIST_MULTI_FN_TYPE  = ZOLTAN_FN_TYPES(3_Zoltan_INT), &
   ZOLTAN_GEOM_MULTI_FN_TYPE       = ZOLTAN_FN_TYPES(5_Zoltan_INT), &
   ZOLTAN_GEOM_FN_TYPE             = ZOLTAN_FN_TYPES(6_Zoltan_INT), &
   ZOLTAN_OBJ_LIST_FN_TYPE         = ZOLTAN_FN_TYPES(8_Zoltan_INT), &
   ZOLTAN_BORDER_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(12_Zoltan_INT), &
   ZOLTAN_PRE_MIGRATE_PP_FN_TYPE   = ZOLTAN_FN_TYPES(15_Zoltan_INT), &
   ZOLTAN_MID_MIGRATE_PP_FN_TYPE   = ZOLTAN_FN_TYPES(16_Zoltan_INT), &
   ZOLTAN_POST_MIGRATE_PP_FN_TYPE  = ZOLTAN_FN_TYPES(17_Zoltan_INT), &
   ZOLTAN_PRE_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(18_Zoltan_INT), &
   ZOLTAN_MID_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(19_Zoltan_INT), &
   ZOLTAN_POST_MIGRATE_FN_TYPE     = ZOLTAN_FN_TYPES(20_Zoltan_INT), &
   ZOLTAN_PACK_OBJ_FN_TYPE         = ZOLTAN_FN_TYPES(22_Zoltan_INT), &
   ZOLTAN_UNPACK_OBJ_FN_TYPE       = ZOLTAN_FN_TYPES(23_Zoltan_INT), &
   ZOLTAN_COARSE_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(25_Zoltan_INT), &
   ZOLTAN_CHILD_LIST_FN_TYPE       = ZOLTAN_FN_TYPES(29_Zoltan_INT), &
   ZOLTAN_CHILD_WEIGHT_FN_TYPE     = ZOLTAN_FN_TYPES(30_Zoltan_INT), &
   ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE   = ZOLTAN_FN_TYPES(31_Zoltan_INT), &
   ZOLTAN_PACK_OBJ_MULTI_FN_TYPE   = ZOLTAN_FN_TYPES(32_Zoltan_INT), &
   ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE = ZOLTAN_FN_TYPES(33_Zoltan_INT), &
   ZOLTAN_PART_MULTI_FN_TYPE       = ZOLTAN_FN_TYPES(35_Zoltan_INT), &
   ZOLTAN_PROC_NAME_FN_TYPE        = ZOLTAN_FN_TYPES(36_Zoltan_INT), &
   ZOLTAN_HG_SIZE_CS_FN_TYPE       = ZOLTAN_FN_TYPES(37_Zoltan_INT), &
   ZOLTAN_HG_CS_FN_TYPE            = ZOLTAN_FN_TYPES(38_Zoltan_INT), &
   ZOLTAN_HG_SIZE_EDGE_WTS_FN_TYPE = ZOLTAN_FN_TYPES(39_Zoltan_INT), &
   ZOLTAN_HG_EDGE_WTS_FN_TYPE      = ZOLTAN_FN_TYPES(40_Zoltan_INT), &
   ZOLTAN_NUM_FIXED_OBJ_FN_TYPE    = ZOLTAN_FN_TYPES(41_Zoltan_INT), &
   ZOLTAN_FIXED_OBJ_LIST_FN_TYPE   = ZOLTAN_FN_TYPES(42_Zoltan_INT), &
   ZOLTAN_HIER_METHOD_FN_TYPE      = ZOLTAN_FN_TYPES(45_Zoltan_INT)

! Backward compatibility with v3.0
type(ZOLTAN_FN_TYPEF), parameter :: &
   ZOLTAN_PARTITION_FN_TYPE        = ZOLTAN_FN_TYPEF(34_Zoltan_INT), &
   ZOLTAN_HIER_PARTITION_FN_TYPE   = ZOLTAN_FN_TYPEF(44_Zoltan_INT)
type(ZOLTAN_FN_TYPES), parameter :: &
   ZOLTAN_PARTITION_MULTI_FN_TYPE  = ZOLTAN_FN_TYPES(35_Zoltan_INT)

! Type of refinement used when building a refinement tree
! These values must agree with the values in zoltan.h

integer(Zoltan_INT), parameter :: &
  ZOLTAN_OTHER_REF     = 0_Zoltan_INT, &
  ZOLTAN_IN_ORDER      = 1_Zoltan_INT, &
  ZOLTAN_TRI_BISECT    = 2_Zoltan_INT, &
  ZOLTAN_QUAD_QUAD     = 3_Zoltan_INT, &
  ZOLTAN_HEX3D_OCT     = 4_Zoltan_INT

! Error codes for LB library
! These values must agree with the values in zoltan.h

integer(Zoltan_INT), parameter :: &
   ZOLTAN_OK     =  0_Zoltan_INT, &
   ZOLTAN_WARN   =  1_Zoltan_INT, &
   ZOLTAN_FATAL  = -1_Zoltan_INT, &
   ZOLTAN_MEMERR = -2_Zoltan_INT

integer(Zoltan_INT), parameter :: &
   ZOLTAN_COMPRESSED_EDGE   = 1_Zoltan_INT, &
   ZOLTAN_COMPRESSED_VERTEX = 2_Zoltan_INT

!--------------------------------------------------------------------------
! defined constants for internal use

integer, parameter :: stderr = 6

!--------------------------------------------------------------------------
! interface blocks for the C wrapper functions

interface
subroutine Zfw_Get_Address_int(arg,ret_addr)
use zoltan_types
use zoltan_user_data
integer(Zoltan_INT) :: arg
integer(Zoltan_INT_PTR), intent(out) :: ret_addr
end subroutine Zfw_Get_Address_int
end interface

interface
subroutine Zfw_Get_Address_struct(arg,ret_addr)
use zoltan_types
use zoltan_user_data
type(Zoltan_Struct) :: arg
integer(Zoltan_INT_PTR), intent(out) :: ret_addr
end subroutine Zfw_Get_Address_struct
end interface

interface
function Zfw_Initialize(ver)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Initialize
real(Zoltan_FLOAT), intent(out) :: ver
end function Zfw_Initialize
end interface

interface
function Zfw_Initialize1(argc,argv,starts,ver)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Initialize1
integer(Zoltan_INT), intent(in) :: argc
integer(Zoltan_INT), dimension(*), intent(in) :: argv, starts
real(Zoltan_FLOAT), intent(out) :: ver
end function Zfw_Initialize1
end interface

interface
subroutine Zfw_Create(communicator,zz,nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer, intent(in) :: communicator
integer(Zoltan_INT), dimension(*), intent(out) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
end subroutine Zfw_Create
end interface

interface
subroutine Zfw_Copy(zzIn, zzOut, nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT), dimension(*), intent(in) :: zzIn
integer(Zoltan_INT), dimension(*), intent(out) :: zzOut
integer(Zoltan_INT), intent(in) :: nbytes
end subroutine Zfw_Copy
end interface

interface
function Zfw_Copy_To(zz1, zz2, nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Copy_To
integer(Zoltan_INT), dimension(*), intent(in) :: zz1, zz2
integer(Zoltan_INT), intent(in) :: nbytes
end function Zfw_Copy_To
end interface

interface
subroutine Zfw_Destroy(zz,nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
end subroutine Zfw_Destroy
end interface

interface
function Zfw_Align(size)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Align
integer(Zoltan_INT) :: size
end function Zfw_Align
end interface

interface
subroutine Zfw_Memory_Stats()
use zoltan_types
use zoltan_user_data
implicit none
end subroutine Zfw_Memory_Stats
end interface

interface
function Zfw_Set_Fn0f(zz,nbytes,fn_type,fn_ptr)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn0f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
end function Zfw_Set_Fn0f
end interface

interface
function Zfw_Set_Fn0s(zz,nbytes,fn_type,fn_ptr)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn0s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
end function Zfw_Set_Fn0s
end interface

interface
function Zfw_Set_Fn1f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn1f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT), dimension(*), intent(in) :: data
end function Zfw_Set_Fn1f
end interface

interface
function Zfw_Set_Fn1s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn1s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
integer(Zoltan_INT), dimension(*), intent(in) :: data
end function Zfw_Set_Fn1s
end interface

interface
function Zfw_Set_Fn2f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn2f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_FLOAT), dimension(*), intent(in) :: data
end function Zfw_Set_Fn2f
end interface

interface
function Zfw_Set_Fn2s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn2s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
real(Zoltan_FLOAT), dimension(*), intent(in) :: data
end function Zfw_Set_Fn2s
end interface

interface
function Zfw_Set_Fn3f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn3f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_DOUBLE), dimension(*), intent(in) :: data
end function Zfw_Set_Fn3f
end interface

interface
function Zfw_Set_Fn3s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn3s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
real(Zoltan_DOUBLE), dimension(*), intent(in) :: data
end function Zfw_Set_Fn3s
end interface

interface
function Zfw_Set_Fn4f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn4f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_1), intent(in) :: data
end function Zfw_Set_Fn4f
end interface

interface
function Zfw_Set_Fn4s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn4s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
type(Zoltan_User_Data_1), intent(in) :: data
end function Zfw_Set_Fn4s
end interface

interface
function Zfw_Set_Fn5f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn5f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_2), intent(in) :: data
end function Zfw_Set_Fn5f
end interface

interface
function Zfw_Set_Fn5s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn5s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
type(Zoltan_User_Data_2), intent(in) :: data
end function Zfw_Set_Fn5s
end interface

interface
function Zfw_Set_Fn6f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn6f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_3), intent(in) :: data
end function Zfw_Set_Fn6f
end interface

interface
function Zfw_Set_Fn6s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn6s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
type(Zoltan_User_Data_3), intent(in) :: data
end function Zfw_Set_Fn6s
end interface

interface
function Zfw_Set_Fn7f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn7f
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_4), intent(in) :: data
end function Zfw_Set_Fn7f
end interface

interface
function Zfw_Set_Fn7s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn7s
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, fn_type
external fn_ptr
type(Zoltan_User_Data_4), intent(in) :: data
end function Zfw_Set_Fn7s
end interface

interface
function Zfw_Set_Param(zz,nbytes,param_name,param_name_len, &
                          new_value,new_value_len)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Param
integer(Zoltan_INT), dimension(*), intent(in) :: zz, param_name, new_value
integer(Zoltan_INT), intent(in) :: nbytes, param_name_len, new_value_len
end function Zfw_Set_Param
end interface

interface
function Zfw_Set_Param_Vec(zz,nbytes,param_name,param_name_len, &
                          new_value,new_value_len,index)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Param_Vec
integer(Zoltan_INT), dimension(*), intent(in) :: zz, param_name, new_value
integer(Zoltan_INT), intent(in) :: nbytes, param_name_len, new_value_len
integer(Zoltan_INT), intent(in) :: index
end function Zfw_Set_Param_Vec
end interface

interface
function Zfw_LB_Partition(zz,nbytes,changes,num_gid_entries,num_lid_entries, &
                num_import,import_global_ids, &
                import_local_ids,import_procs,import_to_part,num_export, &
                export_global_ids,export_local_ids,export_procs,export_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Partition
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), pointer, dimension(:) :: import_to_part, export_to_part
end function Zfw_LB_Partition
end interface

interface
function Zfw_LB_Eval(zz,nbytes,print_stats)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Eval
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes, print_stats
end function Zfw_LB_Eval
end interface

interface
function Zfw_LB_Set_Part_Sizes(zz,nbytes,global_part,len,partids,&
                               wgtidx,partsizes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Set_Part_Sizes
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes,global_part,len,partids(*),wgtidx(*)
real(Zoltan_FLOAT), intent(in) :: partsizes(*)
end function Zfw_LB_Set_Part_Sizes
end interface

interface
function Zfw_LB_Point_Assign(zz,nbytes,coords,proc)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Point_Assign
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
real(Zoltan_DOUBLE), dimension(*), intent(in) :: coords
integer(Zoltan_INT), intent(out) :: proc
end function Zfw_LB_Point_Assign
end interface

interface
function Zfw_LB_Point_PP_Assign(zz,nbytes,coords,proc,part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Point_PP_Assign
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
real(Zoltan_DOUBLE), dimension(*), intent(in) :: coords
integer(Zoltan_INT), intent(out) :: proc
integer(Zoltan_INT), intent(out) :: part
end function Zfw_LB_Point_PP_Assign
end interface

interface
function Zfw_LB_Box_Assign(zz,nbytes,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Box_Assign
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
real(Zoltan_DOUBLE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), dimension(*), intent(out) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
end function Zfw_LB_Box_Assign
end interface

interface
function Zfw_LB_Box_PP_Assign(zz,nbytes,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs,parts,numparts)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Box_PP_Assign
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
real(Zoltan_DOUBLE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), dimension(*), intent(out) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
integer(Zoltan_INT), dimension(*), intent(out) :: parts
integer(Zoltan_INT), intent(out) :: numparts
end function Zfw_LB_Box_PP_Assign
end interface

interface
function Zfw_Invert_Lists(zz,nbytes, &
                       num_input,input_global_ids,input_local_ids, &
                       input_procs,input_to_part, &
                       num_output,output_global_ids,output_local_ids, &
                       output_procs,output_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Invert_Lists
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(in) :: num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), dimension(*), intent(in) :: input_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_global_ids
integer(Zoltan_INT), dimension(*), intent(in) :: input_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_local_ids
integer(Zoltan_INT), dimension(*), intent(in) :: input_procs
integer(Zoltan_INT), pointer, dimension(:) :: output_procs
integer(Zoltan_INT), dimension(*), intent(in) :: input_to_part
integer(Zoltan_INT), pointer, dimension(:) :: output_to_part
end function Zfw_Invert_Lists
end interface

interface
function Zfw_Compute_Destinations(zz,nbytes, &
                       num_input,input_global_ids, &
                       input_local_ids,input_procs,num_output, &
                       output_global_ids,output_local_ids,output_procs)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Compute_Destinations
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(in) :: num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), dimension(*), intent(in) :: input_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_global_ids
integer(Zoltan_INT), dimension(*), intent(in) :: input_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_local_ids
integer(Zoltan_INT), dimension(*), intent(in) :: input_procs
integer(Zoltan_INT), pointer, dimension(:) :: output_procs
end function Zfw_Compute_Destinations
end interface

interface
function Zfw_Migrate(zz,nbytes, &
                     num_import,import_global_ids,import_local_ids, &
                     import_procs,import_to_part, &
                     num_export,export_global_ids,export_local_ids, &
                     export_procs,export_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Migrate
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(in) :: num_import, num_export
integer(Zoltan_INT), dimension(*), intent(in) :: import_global_ids, export_global_ids
integer(Zoltan_INT), dimension(*), intent(in) :: import_local_ids, export_local_ids
integer(Zoltan_INT), dimension(*), intent(in) :: import_procs, export_procs
integer(Zoltan_INT), dimension(*), intent(in) :: import_to_part, export_to_part
end function Zfw_Migrate
end interface

interface
function Zfw_Help_Migrate(zz,nbytes, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Help_Migrate
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(in) :: num_import, num_export
integer(Zoltan_INT), dimension(*), intent(in) :: import_global_ids, export_global_ids
integer(Zoltan_INT), dimension(*), intent(in) :: import_local_ids, export_local_ids
integer(Zoltan_INT), dimension(*), intent(in) :: import_procs, export_procs
end function Zfw_Help_Migrate
end interface

interface
function Zfw_Order(zz,nbytes,num_gid_entries,num_obj, &
                   gids,perm)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Order
INTEGER(Zoltan_INT), dimension(*), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: nbytes
INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*)
INTEGER(Zoltan_INT) :: perm(*)
end function Zfw_Order
end interface

interface
function Zfw_Color(zz,nbytes,num_gid_entries,num_obj, &
                   gids,color_exp)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Color
INTEGER(Zoltan_INT), dimension(*), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: nbytes
INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*)
INTEGER(Zoltan_INT) :: color_exp(*)
end function Zfw_Color
end interface

interface
function Zfw_Color_Test(zz,nbytes,num_gid_entries,num_lid_entries,num_obj, &
                   gids,lids,color_exp)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Color_Test
INTEGER(Zoltan_INT), dimension(*), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: nbytes
INTEGER(Zoltan_INT), INTENT(OUT) :: num_gid_entries, num_lid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*), lids(*)
INTEGER(Zoltan_INT) :: color_exp(*)
end function Zfw_Color_Test
end interface

interface
function Zfw_Generate_Files(zz,nbytes,filename,filename_len, &
                          base_index, gen_geom, gen_graph, gen_hg)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Generate_Files
integer(Zoltan_INT), dimension(*), intent(in) :: zz, filename
integer(Zoltan_INT), intent(in) :: nbytes, filename_len, base_index
integer(Zoltan_INT), intent(in) :: gen_geom, gen_graph, gen_hg
end function Zfw_Generate_Files
end interface

interface
function Zfw_RCB_Box(zz,nbytes,part,ndim,xmin,ymin,zmin,xmax,ymax,zmax)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_RCB_Box
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(in) :: part
integer(Zoltan_INT), intent(out) :: ndim
real(Zoltan_DOUBLE), intent(out) :: xmin,ymin,zmin,xmax,ymax,zmax
end function Zfw_RCB_Box
end interface

interface
subroutine Zfw_Register_Fort_Malloc(malloc_int,free_int,&
      fort_malloc_set_struct)
use zoltan_types
use zoltan_user_data
implicit none
external malloc_int,free_int,fort_malloc_set_struct
end subroutine Zfw_Register_Fort_Malloc
end interface

interface
function Zfw_Get_Wgt_Dim(zz,nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Get_Wgt_Dim
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
end function Zfw_Get_Wgt_Dim
end interface

interface
function Zfw_Get_Comm_Dim(zz,nbytes)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Get_Comm_Dim
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
end function Zfw_Get_Comm_Dim
end interface

interface
subroutine Zfw_Reftree_Get_Child_Order(zz,nbytes,order,ierr)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT), dimension(*), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: nbytes
integer(Zoltan_INT), intent(inout), dimension(*) :: order
integer(Zoltan_INT), intent(out) :: ierr
end subroutine Zfw_Reftree_Get_Child_Order
end interface

!--------------------------------------------------------------------------
! generic names for the Fortran wrapper procedures

interface Zoltan_Initialize
   module procedure Zf90_Initialize
   module procedure Zf90_Initialize1
end interface

interface Zoltan_Create
   module procedure Zf90_Create
end interface

interface Zoltan_Copy
   module procedure Zf90_Copy
end interface

interface Zoltan_Copy_To
   module procedure Zf90_Copy_To
end interface

interface Zoltan_Destroy
   module procedure Zf90_Destroy
end interface

interface Zoltan_Get_Struct_Addr
   module procedure Zf90_Get_Struct_Addr
end interface

interface Zoltan_Align
   module procedure Zf90_Align
end interface

interface Zoltan_Memory_Stats
   module procedure Zf90_Memory_Stats
end interface

interface Zoltan_Set_Fn
   module procedure Zf90_Set_Fn0f
   module procedure Zf90_Set_Fn1f
   module procedure Zf90_Set_Fn2f
   module procedure Zf90_Set_Fn3f
   module procedure Zf90_Set_Fn4f
   module procedure Zf90_Set_Fn5f
   module procedure Zf90_Set_Fn6f
   module procedure Zf90_Set_Fn7f
   module procedure Zf90_Set_Fn0s
   module procedure Zf90_Set_Fn1s
   module procedure Zf90_Set_Fn2s
   module procedure Zf90_Set_Fn3s
   module procedure Zf90_Set_Fn4s
   module procedure Zf90_Set_Fn5s
   module procedure Zf90_Set_Fn6s
   module procedure Zf90_Set_Fn7s
end interface

interface Zoltan_Set_Param
   module procedure Zf90_Set_Param
end interface

interface Zoltan_Set_Param_Vec
   module procedure Zf90_Set_Param_Vec
end interface

interface Zoltan_LB_Partition
   module procedure Zf90_LB_Partition
end interface

interface Zoltan_LB_Balance
   module procedure Zf90_LB_Balance
end interface

interface Zoltan_LB_Eval
   module procedure Zf90_LB_Eval
end interface

interface Zoltan_LB_Free_Part
   module procedure Zf90_LB_Free_Part
end interface

interface Zoltan_LB_Free_Data
   module procedure Zf90_LB_Free_Data
end interface

interface Zoltan_LB_Set_Part_Sizes
   module procedure Zf90_LB_Set_Part_Sizes
end interface

interface Zoltan_LB_Point_Assign
   module procedure Zf90_LB_Point_Assign
end interface

interface Zoltan_LB_Point_PP_Assign
   module procedure Zf90_LB_Point_PP_Assign
end interface

interface Zoltan_LB_Box_Assign
   module procedure Zf90_LB_Box_Assign
end interface

interface Zoltan_LB_Box_PP_Assign
   module procedure Zf90_LB_Box_PP_Assign
end interface

interface Zoltan_Invert_Lists
   module procedure Zf90_Invert_Lists
end interface

interface Zoltan_Compute_Destinations
   module procedure Zf90_Compute_Destinations
end interface

interface Zoltan_Migrate
   module procedure Zf90_Migrate
end interface

interface Zoltan_Help_Migrate
   module procedure Zf90_Help_Migrate
end interface

interface Zoltan_RCB_Box
   module procedure Zf90_RCB_Box
end interface

interface Zoltan_Order
   module procedure Zf90_Order
end interface

interface Zoltan_Color
   module procedure Zf90_Color
end interface

interface Zoltan_Color_Test
   module procedure Zf90_Color_Test
end interface

interface Zoltan_Generate_Files
   module procedure Zf90_Generate_Files
end interface

interface Zoltan_Get_Child_Order
   module procedure Zf90_Reftree_Get_Child_Order
end interface

INCLUDE "set_numgeom.if"
INCLUDE "set_geommulti.if"
INCLUDE "set_geom.if"
INCLUDE "set_partition.if"
INCLUDE "set_partitionmulti.if"
INCLUDE "set_numedges.if"
INCLUDE "set_numedgesmulti.if"
INCLUDE "set_edgelist.if"
INCLUDE "set_edgelistmulti.if"
INCLUDE "set_numobj.if"
INCLUDE "set_objlist.if"
INCLUDE "set_firstobj.if"
INCLUDE "set_nextobj.if"
INCLUDE "set_numborderobj.if"
INCLUDE "set_borderobjlist.if"
INCLUDE "set_firstborderobj.if"
INCLUDE "set_nextborderobj.if"
INCLUDE "set_premigratepp.if"
INCLUDE "set_midmigratepp.if"
INCLUDE "set_postmigratepp.if"
INCLUDE "set_premigrate.if"
INCLUDE "set_midmigrate.if"
INCLUDE "set_postmigrate.if"
INCLUDE "set_objsize.if"
INCLUDE "set_packobj.if"
INCLUDE "set_unpackobj.if"
INCLUDE "set_objsizemulti.if"
INCLUDE "set_packobjmulti.if"
INCLUDE "set_unpackobjmulti.if"
INCLUDE "set_numcoarseobj.if"
INCLUDE "set_coarseobjlist.if"
INCLUDE "set_firstcoarseobj.if"
INCLUDE "set_nextcoarseobj.if"
INCLUDE "set_numchild.if"
INCLUDE "set_childlist.if"
INCLUDE "set_childweight.if"
INCLUDE "set_hgsizecs.if"
INCLUDE "set_hgsizeedgeweights.if"
INCLUDE "set_hgcs.if"
INCLUDE "set_hgedgeweights.if"
INCLUDE "set_numfixedobj.if"
INCLUDE "set_fixedobjlist.if"
INCLUDE "set_hiernumlevels.if"
INCLUDE "set_hierpartition.if"
INCLUDE "set_hiermethod.if"

contains

!--------------------------------------------------------------------------
! Utilities
!--------------------------------------------------------------------------

subroutine fort_malloc_int(array,n,ret_addr)
! This gets called by the C special_malloc to do the allocation
integer(Zoltan_INT), pointer :: array(:)
integer(Zoltan_INT), intent(in) :: n
integer(Zoltan_INT_PTR), intent(out) :: ret_addr
integer :: stat
! Allocate the space
allocate(array(n),stat=stat)
if (stat==0) then
! Send the address of the allocated space to C
   call Zfw_Get_Address_int(array(1),ret_addr)
else
   write(stderr,*) "Error: out of memory during allocation from Fortran"
   ret_addr = 0
endif
end subroutine fort_malloc_int

subroutine fort_free_int(array)
! This gets called by the C special_free to do the deallocation
integer(Zoltan_INT), pointer :: array(:)
integer :: stat
deallocate(array,stat=stat)
if (stat /= 0) then
   write(stderr,*) "Warning: failed to deallocate memory from Fortran"
endif
end subroutine fort_free_int

subroutine fort_malloc_set_struct(struct_addr,ret_addr)
! This routine is called from C to allocate a type(Zoltan_Struct) variable
! and set it to correspond to a C Zoltan_Struct.  The address of the C
! Zoltan_Struct is passed in through struct_addr as an array of integers,
! each containing one byte of the address.  The address of the Fortran
! type(Zoltan_Struct) is returned in ret_addr.
integer(Zoltan_INT), intent(in) :: struct_addr(*)
integer(Zoltan_INT_PTR), intent(out) :: ret_addr
type(Zoltan_Struct), save :: new_struct
integer :: i

! copy the address of the C structure into the Fortran structure
   do i=1,Zoltan_PTR_LENGTH
      new_struct%addr%addr(i:i) = char(struct_addr(i))
   end do
! send the address of the allocated space to C.  I don't think we need a
! different routine than the one used for integers, because it is only
! using the address of the argument
   call Zfw_Get_Address_struct(new_struct,ret_addr)

end subroutine fort_malloc_set_struct

!--------------------------------------------------------------------------
! Fortran wrapper procedures
!--------------------------------------------------------------------------

function Zf90_Initialize(ver)
integer(Zoltan_INT) :: Zf90_Initialize
real(Zoltan_FLOAT), intent(out) :: ver
call Zfw_Register_Fort_Malloc(fort_malloc_int,fort_free_int, &
fort_malloc_set_struct)
Zf90_Initialize = Zfw_Initialize(ver)
end function Zf90_Initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Initialize1(argc,argv,ver)
integer(Zoltan_INT) :: Zf90_Initialize1
integer(Zoltan_INT), intent(in) :: argc
character(len=*), dimension(*), intent(in) :: argv
real(Zoltan_FLOAT), intent(out) :: ver
integer(Zoltan_INT), allocatable, dimension(:) :: int_argv,starts
integer(Zoltan_INT) :: i, j, leng
call Zfw_Register_Fort_Malloc(fort_malloc_int,fort_free_int, &
fort_malloc_set_struct)
allocate(starts(argc+1), int_argv(len(argv(1))*argc))
starts(1) = 1
do i=1,argc
   leng = len_trim(argv(i))
   do j=1,leng
      int_argv(j+starts(i)-1) = ichar(argv(i)(j:j))
   end do
   starts(i+1) = starts(i) + leng
end do
Zf90_Initialize1 = Zfw_Initialize1(argc,int_argv,starts,ver)
deallocate(starts,int_argv)
end function Zf90_Initialize1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Create(communicator)
type(Zoltan_Struct), pointer :: Zf90_Create
integer, intent(in) :: communicator
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz
integer(Zoltan_INT) :: nbytes
integer :: i
logical :: isnull
allocate(Zf90_Create)
nbytes = Zoltan_PTR_LENGTH
call Zfw_Create(communicator,zz,nbytes)
do i=1,Zoltan_PTR_LENGTH
   Zf90_Create%addr%addr(i:i) = char(zz(i))
end do
isnull = (Zf90_Create%addr == Zoltan_NULL_PTR)
if (isnull) then
   deallocate(Zf90_Create)
   nullify(Zf90_Create)
endif
end function Zf90_Create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Copy(zz_from)
type(Zoltan_Struct), pointer :: Zf90_Copy
type(Zoltan_Struct), intent(in) ::  zz_from
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_to, zz_addr
integer(Zoltan_INT) :: nbytes, i
logical :: isnull
allocate(Zf90_Copy)
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz_from%addr%addr(i:i))
end do
call Zfw_Copy(zz_addr, zz_to, nbytes)
do i=1,Zoltan_PTR_LENGTH
   Zf90_Copy%addr%addr(i:i) = char(zz_to(i))
end do
isnull = (Zf90_Copy%addr == Zoltan_NULL_PTR)
if (isnull) then
   deallocate(Zf90_Copy)
   nullify(Zf90_Copy)
endif
end function Zf90_Copy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Copy_To(zz_to, zz_from)
integer(Zoltan_INT) :: Zf90_Copy_To
type(Zoltan_Struct), pointer :: zz_to, zz_from
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr_to, zz_addr_from
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr_to(i)   = ichar(zz_to%addr%addr(i:i))
   zz_addr_from(i) = ichar(zz_from%addr%addr(i:i))
end do
Zf90_Copy_To = Zfw_Copy_To(zz_addr_to,zz_addr_from,nbytes)
end function Zf90_Copy_To

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Zf90_Destroy(zz)
type(Zoltan_Struct), pointer :: zz
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
call Zfw_Destroy(zz_addr,nbytes)
deallocate(zz)
nullify(zz)
end subroutine Zf90_Destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Zf90_Get_Struct_Addr(zz,zz_addr)
type(Zoltan_Struct), pointer :: zz
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
end subroutine Zf90_Get_Struct_Addr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Align(size)
integer(Zoltan_INT) :: Zf90_Align
integer(Zoltan_INT) :: size
Zf90_Align = Zfw_Align(size)
end function Zf90_Align

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Zf90_Memory_Stats()
call Zfw_Memory_Stats()
end subroutine Zf90_Memory_Stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn0f(zz,fn_type,fn_ptr)
integer(Zoltan_INT) :: Zf90_Set_Fn0f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn0f = Zfw_Set_Fn0f(zz_addr,nbytes,fn_type%choice,fn_ptr)
end function Zf90_Set_Fn0f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn0s(zz,fn_type,fn_ptr)
integer(Zoltan_INT) :: Zf90_Set_Fn0s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn0s = Zfw_Set_Fn0s(zz_addr,nbytes,fn_type%choice,fn_ptr)
end function Zf90_Set_Fn0s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn1f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn1f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn1f = Zfw_Set_Fn1f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn1f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn1s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn1s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
integer(Zoltan_INT), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn1s = Zfw_Set_Fn1s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn1s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn2f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn2f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_FLOAT), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn2f = Zfw_Set_Fn2f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn2f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn2s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn2s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
real(Zoltan_FLOAT), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn2s = Zfw_Set_Fn2s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn2s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn3f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn3f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_DOUBLE), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn3f = Zfw_Set_Fn3f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn3f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn3s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn3s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
real(Zoltan_DOUBLE), intent(in) :: data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn3s = Zfw_Set_Fn3s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn3s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn4f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn4f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_1), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn4f = Zfw_Set_Fn4f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn4f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn4s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn4s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
type(Zoltan_User_Data_1), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn4s = Zfw_Set_Fn4s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn4s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn5f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn5f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn5f = Zfw_Set_Fn5f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn5f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn5s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn5s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
type(Zoltan_User_Data_2), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn5s = Zfw_Set_Fn5s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn5s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn6f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn6f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_3), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn6f = Zfw_Set_Fn6f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn6f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn6s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn6s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
type(Zoltan_User_Data_3), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn6s = Zfw_Set_Fn6s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn6s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn7f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn7f
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPEF), intent(in) :: fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_4), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn7f = Zfw_Set_Fn7f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn7f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Fn7s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn7s
type(Zoltan_Struct), intent(in) :: zz
type(ZOLTAN_FN_TYPES), intent(in) :: fn_type
external fn_ptr
type(Zoltan_User_Data_4), intent(in) :: data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Set_Fn7s = Zfw_Set_Fn7s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
end function Zf90_Set_Fn7s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Param(zz,param_name,new_value)
integer(Zoltan_INT) :: Zf90_Set_Param
type(Zoltan_Struct), intent(in) :: zz
character(len=*), intent(in) :: param_name, new_value
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT), dimension(len_trim(param_name)) :: int_param_name
integer(Zoltan_INT), dimension(len_trim(new_value)) :: int_new_value
integer(Zoltan_INT) :: nbytes, param_name_len, new_value_len, i
nbytes = Zoltan_PTR_LENGTH
param_name_len = len_trim(param_name)
new_value_len = len_trim(new_value)
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
do i=1,param_name_len
   int_param_name(i) = ichar(param_name(i:i))
end do
do i=1,new_value_len
   int_new_value(i) = ichar(new_value(i:i))
end do
Zf90_Set_Param = Zfw_Set_Param(zz_addr,nbytes,int_param_name, &
                                    param_name_len,int_new_value,new_value_len)
end function Zf90_Set_Param

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Param_Vec(zz,param_name,new_value,index)
integer(Zoltan_INT) :: Zf90_Set_Param_Vec
type(Zoltan_Struct), intent(in) :: zz
character(len=*), intent(in) :: param_name, new_value
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT), dimension(len_trim(param_name)) :: int_param_name
integer(Zoltan_INT), dimension(len_trim(new_value)) :: int_new_value
integer(Zoltan_INT) :: index
integer(Zoltan_INT) :: nbytes, param_name_len, new_value_len, i
nbytes = Zoltan_PTR_LENGTH
param_name_len = len_trim(param_name)
new_value_len = len_trim(new_value)
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
do i=1,param_name_len
   int_param_name(i) = ichar(param_name(i:i))
end do
do i=1,new_value_len
   int_new_value(i) = ichar(new_value(i:i))
end do
Zf90_Set_Param_Vec = Zfw_Set_Param_Vec(zz_addr,nbytes,int_param_name, &
                     param_name_len,int_new_value,new_value_len,index)
end function Zf90_Set_Param_Vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Partition(zz,changes,num_gid_entries,num_lid_entries, &
                 num_import,import_global_ids, &
                 import_local_ids,import_procs,import_to_part,num_export, &
                 export_global_ids,export_local_ids,export_procs,export_to_part)
integer(Zoltan_INT) :: Zf90_LB_Partition
type(Zoltan_Struct), intent(in) :: zz
logical, intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), pointer, dimension(:) :: import_to_part, export_to_part
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i, int_changes
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Partition = Zfw_LB_Partition(zz_addr,nbytes,int_changes, &
                             num_gid_entries, num_lid_entries, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,import_to_part, &
                             num_export,export_global_ids, &
                             export_local_ids,export_procs,export_to_part)
changes = .not.(int_changes==0)
end function Zf90_LB_Partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Balance(zz,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_LB_Balance
type(Zoltan_Struct), intent(in) :: zz
logical, intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), pointer, dimension(:) :: import_to_part, export_to_part
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i, int_changes
integer :: stat
stat = 0
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
nullify(import_to_part, export_to_part)
Zf90_LB_Balance = Zfw_LB_Partition(zz_addr,nbytes,int_changes, &
                             num_gid_entries, num_lid_entries, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,import_to_part, &
                             num_export,export_global_ids, &
                             export_local_ids,export_procs,export_to_part)

! Do not return import_to_part, export_to_part.
! Deallocate them if they were allocated.
if (associated(import_to_part)) deallocate(import_to_part,stat=stat)
if (associated(export_to_part)) deallocate(export_to_part,stat=stat)

changes = .not.(int_changes==0)
end function Zf90_LB_Balance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Eval(zz,print_stats)
integer(Zoltan_INT) :: Zf90_LB_Eval
type(Zoltan_Struct), intent(in) :: zz
logical, intent(in) :: print_stats
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i, int_print_stats, dim, edim
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
if (print_stats) then
   int_print_stats = 1
else
   int_print_stats = 0
endif
Zf90_LB_Eval = Zfw_LB_Eval(zz_addr,nbytes,int_print_stats)
end function Zf90_LB_Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Free_Part(global_ids, local_ids, procs, part)
integer(Zoltan_INT) :: Zf90_LB_Free_Part
integer(Zoltan_INT), pointer, dimension(:) :: global_ids
integer(Zoltan_INT), pointer, dimension(:) :: local_ids
integer(Zoltan_INT), pointer, dimension(:) :: procs, part
integer :: stat
stat = 0
Zf90_LB_Free_Part = ZOLTAN_OK
if (associated(global_ids)) deallocate(global_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Part = ZOLTAN_WARN
nullify(global_ids)
if (associated(local_ids)) deallocate(local_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Part = ZOLTAN_WARN
nullify(local_ids)
if (associated(procs)) deallocate(procs,stat=stat)
if (stat /= 0) Zf90_LB_Free_Part = ZOLTAN_WARN
nullify(procs)
if (associated(part)) deallocate(part,stat=stat)
if (stat /= 0) Zf90_LB_Free_Part = ZOLTAN_WARN
nullify(part)
end function Zf90_LB_Free_Part

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Free_Data(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_LB_Free_Data
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
Zf90_LB_Free_Data = ZOLTAN_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) Zf90_LB_Free_Data = ZOLTAN_WARN
nullify(export_procs)
end function Zf90_LB_Free_Data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Set_Part_Sizes(zz,global_part,len,partids,wgtidx,partsizes)
integer(Zoltan_INT) :: Zf90_LB_Set_Part_Sizes
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: global_part,len,partids(*),wgtidx(*)
real(Zoltan_FLOAT), intent(in) :: partsizes(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Set_Part_Sizes = Zfw_LB_Set_Part_Sizes(zz_addr,nbytes,global_part,len,&
                                               partids,wgtidx,partsizes)
end function Zf90_LB_Set_Part_Sizes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Point_Assign(zz,coords,proc)
integer(Zoltan_INT) :: Zf90_LB_Point_Assign
type(Zoltan_Struct), intent(in) :: zz
real(Zoltan_DOUBLE), dimension(*), intent(in) :: coords
integer(Zoltan_INT), intent(out) :: proc
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Point_Assign = Zfw_LB_Point_Assign(zz_addr,nbytes,coords,proc)
end function Zf90_LB_Point_Assign

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Point_PP_Assign(zz,coords,proc,part)
integer(Zoltan_INT) :: Zf90_LB_Point_PP_Assign
type(Zoltan_Struct), intent(in) :: zz
real(Zoltan_DOUBLE), dimension(*), intent(in) :: coords
integer(Zoltan_INT), intent(out) :: proc
integer(Zoltan_INT), intent(out) :: part
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Point_PP_Assign = Zfw_LB_Point_PP_Assign(zz_addr,nbytes,coords,proc,part)
end function Zf90_LB_Point_PP_Assign

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Box_Assign(zz,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
integer(Zoltan_INT) :: Zf90_LB_Box_Assign
type(Zoltan_Struct), intent(in) :: zz
real(Zoltan_DOUBLE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), intent(out), dimension(*) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Box_Assign = Zfw_LB_Box_Assign(zz_addr,nbytes,xmin,ymin,zmin,xmax,ymax, &
                                    zmax,procs,numprocs)
end function Zf90_LB_Box_Assign

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Box_PP_Assign(zz,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs,parts,numparts)
integer(Zoltan_INT) :: Zf90_LB_Box_PP_Assign
type(Zoltan_Struct), intent(in) :: zz
real(Zoltan_DOUBLE), intent(in) :: xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), intent(out), dimension(*) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
integer(Zoltan_INT), intent(out), dimension(*) :: parts
integer(Zoltan_INT), intent(out) :: numparts
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Box_PP_Assign = Zfw_LB_Box_PP_Assign(zz_addr,nbytes,xmin,ymin,zmin,xmax,ymax, &
                                    zmax,procs,numprocs,parts,numparts)
end function Zf90_LB_Box_PP_Assign

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Invert_Lists(zz, &
                       num_input,input_global_ids,input_local_ids, &
                       input_procs,input_to_part, &
                       num_output,output_global_ids,output_local_ids, &
                       output_procs,output_to_part)
integer(Zoltan_INT) :: Zf90_Invert_Lists
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), pointer, dimension(:) :: input_global_ids,output_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: input_local_ids, output_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: input_procs, output_procs
integer(Zoltan_INT), pointer, dimension(:) :: input_to_part, output_to_part
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
if (.not.associated(input_global_ids) .or. .not.associated(input_local_ids) &
    .or. .not.associated(input_procs) .or. .not.associated(input_to_part)) then
   write(stderr,*) "Error from Zoltan_Invert_Lists: input pointers are not associated"
   Zf90_Invert_Lists = ZOLTAN_WARN
   return
endif
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Invert_Lists = Zfw_Invert_Lists(zz_addr,nbytes, &
                             num_input,input_global_ids,input_local_ids, &
                             input_procs,input_to_part, &
                             num_output,output_global_ids, &
                             output_local_ids,output_procs,output_to_part)
end function Zf90_Invert_Lists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Compute_Destinations(zz, &
                       num_input,input_global_ids, &
                       input_local_ids,input_procs,num_output, &
                       output_global_ids,output_local_ids,output_procs)
integer(Zoltan_INT) :: Zf90_Compute_Destinations
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), pointer, dimension(:) :: input_global_ids, output_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: input_local_ids, output_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: input_procs, output_procs
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
if (.not.associated(input_global_ids) .or. .not.associated(input_local_ids) &
    .or. .not.associated(input_procs)) then
   write(stderr,*) "Error from Zoltan_Compute_Destinations: input pointers are not associated"
   Zf90_Compute_Destinations = ZOLTAN_WARN
   return
endif
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Compute_Destinations = Zfw_Compute_Destinations(zz_addr,nbytes, &
                             num_input,input_global_ids,input_local_ids, &
                             input_procs,num_output,output_global_ids, &
                             output_local_ids,output_procs)
end function Zf90_Compute_Destinations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Migrate(zz, &
                       num_import,import_global_ids,import_local_ids, &
                       import_procs,import_to_part, &
                       num_export,export_global_ids,export_local_ids, &
                       export_procs,export_to_part)
integer(Zoltan_INT) :: Zf90_Migrate
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), pointer, dimension(:) :: import_to_part, export_to_part
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
logical :: free_import_global_ids, free_import_local_ids, free_import_procs
logical :: free_export_global_ids, free_export_local_ids, free_export_procs
logical :: free_import_to_part, free_export_to_part

if ((num_import.gt.0).and.(.not.associated(import_global_ids) .or. &
                           .not.associated(import_local_ids)  .or. &
                           .not.associated(import_procs))) then
   ! OK if import_to_part is not associated; some methods don't return parts.
   write(stderr,*) "Error from Zoltan_Migrate: import pointers are not associated"
   Zf90_Migrate = ZOLTAN_WARN
   return
endif
if ((num_export.gt.0).and.(.not.associated(export_procs) .or. &
                           .not.associated(export_global_ids) .or. &
                           .not.associated(export_local_ids))) then
   ! OK if export_to_part is not associated; some methods don't return parts.
   write(stderr,*) "Error from Zoltan_Migrate: export pointers are not associated"

   Zf90_Migrate = ZOLTAN_WARN
   return
endif

! generate place-holders to make call to Zfw_Migrate valid;
! can't call it with non-associated arrays, even if we aren't importing
! or exporting items.
free_import_global_ids = .false.
free_import_local_ids  = .false.
free_import_procs      = .false.
free_import_to_part    = .false.
free_export_global_ids = .false.
free_export_local_ids  = .false.
free_export_procs      = .false.
free_export_to_part    = .false.

if (.not.associated(import_global_ids)) then
   free_import_global_ids = .true.
   allocate(import_global_ids(0)) 
endif
if (.not.associated(import_local_ids)) then
   free_import_local_ids = .true.
   allocate(import_local_ids(0)) 
endif
if (.not.associated(import_procs)) then
   free_import_procs = .true.
   allocate(import_procs(0)) 
endif
if (.not.associated(import_to_part)) then
   free_import_to_part = .true.
   allocate(import_to_part(0)) 
endif
if (.not.associated(export_global_ids)) then
   free_export_global_ids = .true.
   allocate(export_global_ids(0)) 
endif
if (.not.associated(export_local_ids)) then
   free_export_local_ids = .true.
   allocate(export_local_ids(0)) 
endif
if (.not.associated(export_procs)) then
   free_export_procs = .true.
   allocate(export_procs(0)) 
endif
if (.not.associated(export_to_part)) then
   free_export_to_part = .true.
   allocate(export_to_part(0)) 
endif

nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Migrate = Zfw_Migrate(zz_addr,nbytes, &
                           num_import,import_global_ids,import_local_ids, &
                           import_procs,import_to_part, &
                           num_export,export_global_ids, &
                           export_local_ids,export_procs,export_to_part)

! clean up the place holders
if (free_import_global_ids) deallocate(import_global_ids)
if (free_import_local_ids) deallocate(import_local_ids)
if (free_import_procs) deallocate(import_procs)
if (free_import_to_part) deallocate(import_to_part)
if (free_export_global_ids) deallocate(export_global_ids)
if (free_export_local_ids) deallocate(export_local_ids)
if (free_export_procs) deallocate(export_procs)
if (free_export_to_part) deallocate(export_to_part)

end function Zf90_Migrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Help_Migrate(zz, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_Help_Migrate
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
logical :: free_import_global_ids, free_import_local_ids, free_import_procs
logical :: free_export_global_ids, free_export_local_ids, free_export_procs

if ((num_import.gt.0).and.(.not.associated(import_global_ids) .or. &
                           .not.associated(import_local_ids)  .or. &
                           .not.associated(import_procs))) then
   write(stderr,*) "Error from Zoltan_Help_Migrate: import pointers are not associated"
   Zf90_Help_Migrate = ZOLTAN_WARN
   return
endif
if ((num_export.gt.0).and.(.not.associated(export_procs) .or. &
                           .not.associated(export_global_ids) .or. &
                           .not.associated(export_local_ids))) then
   write(stderr,*) "Error from Zoltan_Help_Migrate: export pointers are not associated"

   Zf90_Help_Migrate = ZOLTAN_WARN
   return
endif

! generate place-holders to make call to Zfw_Help_Migrate valid;
! can't call it with non-associated arrays, even if we aren't importing
! or exporting items.
free_import_global_ids = .false.
free_import_local_ids  = .false.
free_import_procs      = .false.
free_export_global_ids = .false.
free_export_local_ids  = .false.
free_export_procs      = .false.

if (.not.associated(import_global_ids)) then
   free_import_global_ids = .true.
   allocate(import_global_ids(0)) 
endif
if (.not.associated(import_local_ids)) then
   free_import_local_ids = .true.
   allocate(import_local_ids(0)) 
endif
if (.not.associated(import_procs)) then
   free_import_procs = .true.
   allocate(import_procs(0)) 
endif
if (.not.associated(export_global_ids)) then
   free_export_global_ids = .true.
   allocate(export_global_ids(0)) 
endif
if (.not.associated(export_local_ids)) then
   free_export_local_ids = .true.
   allocate(export_local_ids(0)) 
endif
if (.not.associated(export_procs)) then
   free_export_procs = .true.
   allocate(export_procs(0)) 
endif

nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Help_Migrate = Zfw_Help_Migrate(zz_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)

! clean up the place holders
if (free_import_global_ids) deallocate(import_global_ids)
if (free_import_local_ids) deallocate(import_local_ids)
if (free_import_procs) deallocate(import_procs)
if (free_export_global_ids) deallocate(export_global_ids)
if (free_export_local_ids) deallocate(export_local_ids)
if (free_export_procs) deallocate(export_procs)

end function Zf90_Help_Migrate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Order(zz,num_gid_entries,num_obj,gids,perm)
integer(Zoltan_INT) :: Zf90_Order
TYPE(Zoltan_Struct), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*)
INTEGER(Zoltan_INT) :: perm(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Order = Zfw_Order(zz_addr,nbytes,num_gid_entries,num_obj,&
                       gids,perm)
end function Zf90_Order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Color(zz,num_gid_entries,num_obj,gids,color_exp)
integer(Zoltan_INT) :: Zf90_Color
TYPE(Zoltan_Struct), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: num_gid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*)
INTEGER(Zoltan_INT) :: color_exp(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Color = Zfw_Color(zz_addr,nbytes,num_gid_entries,num_obj,&
                       gids,color_exp)
end function Zf90_Color

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Color_Test(zz,num_gid_entries,num_lid_entries,num_obj,gids,lids,color_exp)
integer(Zoltan_INT) :: Zf90_Color_Test
TYPE(Zoltan_Struct), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(OUT) :: num_gid_entries, num_lid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*), lids(*)
INTEGER(Zoltan_INT) :: color_exp(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Color_Test = Zfw_Color_Test(zz_addr,nbytes,num_gid_entries,num_lid_entries,num_obj,&
                       gids,lids,color_exp)
end function Zf90_Color_Test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Generate_Files(zz,filename,base_index,gen_geom,gen_graph,gen_hg)
integer(Zoltan_INT) :: Zf90_Generate_Files
type(Zoltan_Struct), intent(in) :: zz
character(len=*), intent(in) :: filename
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT), dimension(len_trim(filename)) :: int_filename
integer(Zoltan_INT) :: nbytes, filename_len,base_index,gen_geom,gen_graph,gen_hg
integer(Zoltan_INT) :: i
nbytes = Zoltan_PTR_LENGTH
filename_len = len_trim(filename)
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
do i=1,filename_len
   int_filename(i) = ichar(filename(i:i))
end do
Zf90_Generate_Files = Zfw_Generate_Files(zz_addr,nbytes,int_filename, &
                          filename_len,base_index,gen_geom,gen_graph,gen_hg)
end function Zf90_Generate_Files

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_RCB_Box(zz,part,ndim,xmin,ymin,zmin,xmax,ymax,zmax)
integer(Zoltan_INT) :: Zf90_RCB_Box
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(in) :: part
integer(Zoltan_INT), intent(out) :: ndim
real(Zoltan_DOUBLE), intent(out) :: xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_RCB_Box = Zfw_RCB_Box(zz_addr,nbytes,part,ndim,xmin,ymin,zmin,xmax,ymax, &
                           zmax)
end function Zf90_RCB_Box

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Zf90_Reftree_Get_Child_Order(zz,order,ierr)
type(Zoltan_Struct), intent(in) :: zz
integer(Zoltan_INT), intent(inout), dimension(*) :: order
integer(Zoltan_INT), intent(out) :: ierr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
call Zfw_Reftree_Get_Child_Order(zz_addr,nbytes,order,ierr)
end subroutine Zf90_Reftree_Get_Child_Order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE "set_numgeom.fn"
INCLUDE "set_geommulti.fn"
INCLUDE "set_geom.fn"
INCLUDE "set_partition.fn"
INCLUDE "set_partitionmulti.fn"
INCLUDE "set_numedges.fn"
INCLUDE "set_numedgesmulti.fn"
INCLUDE "set_edgelist.fn"
INCLUDE "set_edgelistmulti.fn"
INCLUDE "set_numobj.fn"
INCLUDE "set_objlist.fn"
INCLUDE "set_firstobj.fn"
INCLUDE "set_nextobj.fn"
INCLUDE "set_numborderobj.fn"
INCLUDE "set_borderobjlist.fn"
INCLUDE "set_firstborderobj.fn"
INCLUDE "set_nextborderobj.fn"
INCLUDE "set_premigratepp.fn"
INCLUDE "set_midmigratepp.fn"
INCLUDE "set_postmigratepp.fn"
INCLUDE "set_premigrate.fn"
INCLUDE "set_midmigrate.fn"
INCLUDE "set_postmigrate.fn"
INCLUDE "set_objsize.fn"
INCLUDE "set_packobj.fn"
INCLUDE "set_unpackobj.fn"
INCLUDE "set_objsizemulti.fn"
INCLUDE "set_packobjmulti.fn"
INCLUDE "set_unpackobjmulti.fn"
INCLUDE "set_numcoarseobj.fn"
INCLUDE "set_coarseobjlist.fn"
INCLUDE "set_firstcoarseobj.fn"
INCLUDE "set_nextcoarseobj.fn"
INCLUDE "set_numchild.fn"
INCLUDE "set_childlist.fn"
INCLUDE "set_childweight.fn"
INCLUDE "set_hgsizecs.fn"
INCLUDE "set_hgsizeedgeweights.fn"
INCLUDE "set_hgcs.fn"
INCLUDE "set_hgedgeweights.fn"
INCLUDE "set_numfixedobj.fn"
INCLUDE "set_fixedobjlist.fn"
INCLUDE "set_hiernumlevels.fn"
INCLUDE "set_hierpartition.fn"
INCLUDE "set_hiermethod.fn"

end module zoltan
