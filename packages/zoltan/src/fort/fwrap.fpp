!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zoltan Dynamic Load-Balancing Library for Parallel Applications            !
! For more info, see the README file in the top-level Zoltan directory.      ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CVS File Information :
!     $RCSfile$
!     $Author$
!     $Date$
!     $Revision$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "sppr_header"
!--------------------------------------------------------------------------
! preprocessor directives to handle special case compilers

#ifdef NASOFTWARE
#define INTENT_IN ::
#define C_BINDING bind(c)
#define PASS_BY_REF ,pass_by("*")
#else
#define INTENT_IN ,intent(in) ::
#define C_BINDING
#define PASS_BY_REF
#endif

module zoltan
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
private

!--------------------------------------------------------------------------
! public entities

public :: &
   LB_INT, &
   LB_FLOAT, &
   LB_DOUBLE, &
   LB_User_Data_1, &
   LB_User_Data_2, &
   LB_User_Data_3, &
   LB_User_Data_4

public :: &
   LB_Struct, &
   ZOLTAN_FN_TYPEF, &
   ZOLTAN_FN_TYPES

public :: &
   ZOLTAN_NUM_EDGES_FN_TYPE, &
   ZOLTAN_EDGE_LIST_FN_TYPE, &
   ZOLTAN_NUM_GEOM_FN_TYPE, &
   ZOLTAN_GEOM_FN_TYPE, &
   ZOLTAN_NUM_OBJ_FN_TYPE, &
   ZOLTAN_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_OBJ_FN_TYPE, &
   ZOLTAN_NUM_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_BORDER_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE, &
   ZOLTAN_PRE_MIGRATE_FN_TYPE, &
   ZOLTAN_MID_MIGRATE_FN_TYPE, &
   ZOLTAN_POST_MIGRATE_FN_TYPE, &
   ZOLTAN_OBJ_SIZE_FN_TYPE, &
   ZOLTAN_PACK_OBJ_FN_TYPE, &
   ZOLTAN_UNPACK_OBJ_FN_TYPE, &
   ZOLTAN_NUM_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_COARSE_OBJ_LIST_FN_TYPE, &
   ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE, &
   ZOLTAN_NUM_CHILD_FN_TYPE, &
   ZOLTAN_CHILD_LIST_FN_TYPE, &
   ZOLTAN_CHILD_WEIGHT_FN_TYPE

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
   Zoltan_Initialize, &
   Zoltan_Create, &
   Zoltan_Destroy, &
   Zoltan_Memory_Stats, &
   Zoltan_Set_Fn, &
   Zoltan_LB_Set_Method, &
   Zoltan_Set_Param, &
   Zoltan_LB_Balance, &
   Zoltan_LB_Eval, &
   Zoltan_LB_Free_Data, &
   Zoltan_LB_Point_Assign, &
   Zoltan_LB_Box_Assign, &
   Zoltan_Compute_Destinations, &
   Zoltan_Help_Migrate

! Registration functions with strict type checking.
public :: &
   Zoltan_Set_Num_Obj_Fn, Zoltan_Set_Obj_List_Fn, &
   Zoltan_Set_First_Obj_Fn, Zoltan_Set_Next_Obj_Fn, &
   Zoltan_Set_Num_Border_Obj_Fn, Zoltan_Set_Border_Obj_List_Fn, &
   Zoltan_Set_First_Border_Obj_Fn, Zoltan_Set_Next_Border_Obj_Fn, &
   Zoltan_Set_Num_Geom_Fn, Zoltan_Set_Geom_Fn, &
   Zoltan_Set_Num_Edges_Fn, Zoltan_Set_Edge_List_Fn, &
   Zoltan_Set_Num_Coarse_Obj_Fn, Zoltan_Set_Coarse_Obj_List_Fn, &
   Zoltan_Set_First_Coarse_Obj_Fn, Zoltan_Set_Next_Coarse_Obj_Fn, &
   Zoltan_Set_Num_Child_Fn, Zoltan_Set_Child_List_Fn, &
   Zoltan_Set_Child_Weight_Fn, &
   Zoltan_Set_Pre_Migrate_Fn, Zoltan_Set_Mid_Migrate_Fn, &
   Zoltan_Set_Post_Migrate_Fn, &
   Zoltan_Set_Obj_Size_Fn, Zoltan_Set_Pack_Obj_Fn, Zoltan_Set_Unpack_Obj_Fn 

public :: &
   Zoltan_Get_Child_Order ! TEMP child_order

!--------------------------------------------------------------------------
! user defined types corresponding to the C structs

type LB_Struct
   private
   sequence
   type(LB_PTR) :: addr
#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
   integer :: dummy
#endif
end type LB_Struct

!--------------------------------------------------------------------------
! defined constants corresponding to Zoltan enumerated types

! Enumerated type used to indicate which function is to be set by LB_Set_Fn.
! These values must agree with those used in the LB_Set_Fn wrapper in cwrap.c

type ZOLTAN_FN_TYPEF
   private
   integer(LB_INT) :: choice
end type ZOLTAN_FN_TYPEF

type ZOLTAN_FN_TYPES
   private
   integer(LB_INT) :: choice
end type ZOLTAN_FN_TYPES

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPEF) :: &
#else
type(ZOLTAN_FN_TYPEF), parameter :: &
#endif
   ZOLTAN_NUM_EDGES_FN_TYPE        = ZOLTAN_FN_TYPEF(0_LB_INT), &
   ZOLTAN_NUM_GEOM_FN_TYPE         = ZOLTAN_FN_TYPEF(2_LB_INT), &
   ZOLTAN_NUM_OBJ_FN_TYPE          = ZOLTAN_FN_TYPEF(4_LB_INT), &
   ZOLTAN_FIRST_OBJ_FN_TYPE        = ZOLTAN_FN_TYPEF(6_LB_INT), &
   ZOLTAN_NEXT_OBJ_FN_TYPE         = ZOLTAN_FN_TYPEF(7_LB_INT), &
   ZOLTAN_NUM_BORDER_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(8_LB_INT), &
   ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(10_LB_INT), &
   ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(11_LB_INT), &
   ZOLTAN_OBJ_SIZE_FN_TYPE         = ZOLTAN_FN_TYPEF(15_LB_INT), &
   ZOLTAN_NUM_COARSE_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(18_LB_INT), &
   ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(20_LB_INT), &
   ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(21_LB_INT), &
   ZOLTAN_NUM_CHILD_FN_TYPE        = ZOLTAN_FN_TYPEF(22_LB_INT)

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPES) :: &
#else
type(ZOLTAN_FN_TYPES), parameter :: &
#endif
   ZOLTAN_EDGE_LIST_FN_TYPE        = ZOLTAN_FN_TYPES(1_LB_INT), &
   ZOLTAN_GEOM_FN_TYPE             = ZOLTAN_FN_TYPES(3_LB_INT), &
   ZOLTAN_OBJ_LIST_FN_TYPE         = ZOLTAN_FN_TYPES(5_LB_INT), &
   ZOLTAN_BORDER_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(9_LB_INT), &
   ZOLTAN_PRE_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(12_LB_INT), &
   ZOLTAN_MID_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(13_LB_INT), &
   ZOLTAN_POST_MIGRATE_FN_TYPE     = ZOLTAN_FN_TYPES(14_LB_INT), &
   ZOLTAN_PACK_OBJ_FN_TYPE         = ZOLTAN_FN_TYPES(16_LB_INT), &
   ZOLTAN_UNPACK_OBJ_FN_TYPE       = ZOLTAN_FN_TYPES(17_LB_INT), &
   ZOLTAN_COARSE_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(19_LB_INT), &
   ZOLTAN_CHILD_LIST_FN_TYPE       = ZOLTAN_FN_TYPES(23_LB_INT), &
   ZOLTAN_CHILD_WEIGHT_FN_TYPE     = ZOLTAN_FN_TYPES(24_LB_INT)

! Type of refinement used when building a refinement tree
! These values must agree with the values in zoltan.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(LB_INT) :: &
#else
integer(LB_INT), parameter :: &
#endif
  ZOLTAN_OTHER_REF     = 0_LB_INT, &
  ZOLTAN_IN_ORDER      = 1_LB_INT, &
  ZOLTAN_TRI_BISECT    = 2_LB_INT, &
  ZOLTAN_QUAD_QUAD     = 3_LB_INT, &
  ZOLTAN_HEX3D_OCT     = 4_LB_INT

! Error codes for LB library
! These values must agree with the values in zoltan.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(LB_INT) :: &
#else
integer(LB_INT), parameter :: &
#endif
   ZOLTAN_OK     =  0_LB_INT, &
   ZOLTAN_WARN   =  1_LB_INT, &
   ZOLTAN_FATAL  = -1_LB_INT, &
   ZOLTAN_MEMERR = -2_LB_INT

!--------------------------------------------------------------------------
! defined constants for internal use

integer, parameter :: stderr = 6

!--------------------------------------------------------------------------
! interface blocks for the C wrapper functions

interface
!NAS$ ALIEN "F77 lb_fw_get_address_int"
subroutine LB_fw_Get_Address_int(arg,ret_addr)
use zoltan_types
use lb_user_const
integer(LB_INT) :: arg
integer(LB_INT_PTR), intent(out) :: ret_addr
end subroutine LB_fw_Get_Address_int
end interface

interface
!NAS$ ALIEN "F77 lb_fw_initialize"
function LB_fw_Initialize(ver)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Initialize
real(LB_FLOAT), intent(out) :: ver
end function LB_fw_Initialize
end interface

interface
!NAS$ ALIEN "F77 lb_fw_initialize1"
function LB_fw_Initialize1(argc,argv,starts,ver)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Initialize1
integer(LB_INT) INTENT_IN argc
integer(LB_INT), dimension(*) INTENT_IN argv, starts
real(LB_FLOAT), intent(out) :: ver
end function LB_fw_Initialize1
end interface

interface
!NAS$ ALIEN "F77 lb_fw_create"
subroutine LB_fw_Create(communicator,lb,nbytes)
use zoltan_types
use lb_user_const
implicit none
integer INTENT_IN communicator
integer(LB_INT), dimension(*), intent(out) :: lb
integer(LB_INT) INTENT_IN nbytes
end subroutine LB_fw_Create
end interface

interface
!NAS$ ALIEN "F77 lb_fw_destroy"
subroutine LB_fw_Destroy(lb,nbytes)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
end subroutine LB_fw_Destroy
end interface

interface
!NAS$ ALIEN "F77 lb_fw_memory_stats"
subroutine LB_fw_Memory_Stats()
use zoltan_types
use lb_user_const
implicit none
end subroutine LB_fw_Memory_Stats
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn0f"
function LB_fw_Set_Fn0f(lb,nbytes,fn_type,fn_ptr)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn0f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
end function LB_fw_Set_Fn0f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn0s"
function LB_fw_Set_Fn0s(lb,nbytes,fn_type,fn_ptr)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn0s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
end function LB_fw_Set_Fn0s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn1f"
function LB_fw_Set_Fn1f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn1f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
integer(LB_INT), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn1f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn1s"
function LB_fw_Set_Fn1s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn1s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
integer(LB_INT), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn1s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn2f"
function LB_fw_Set_Fn2f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn2f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
real(LB_FLOAT), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn2f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn2s"
function LB_fw_Set_Fn2s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn2s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
real(LB_FLOAT), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn2s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn3f"
function LB_fw_Set_Fn3f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn3f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
real(LB_DOUBLE), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn3f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn3s"
function LB_fw_Set_Fn3s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn3s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
real(LB_DOUBLE), dimension(*) INTENT_IN data
end function LB_fw_Set_Fn3s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn4f"
function LB_fw_Set_Fn4f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn4f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
type(LB_User_Data_1) INTENT_IN data
end function LB_fw_Set_Fn4f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn4s"
function LB_fw_Set_Fn4s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn4s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(LB_User_Data_1) INTENT_IN data
end function LB_fw_Set_Fn4s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn5f"
function LB_fw_Set_Fn5f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn5f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
type(LB_User_Data_2) INTENT_IN data
end function LB_fw_Set_Fn5f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn5s"
function LB_fw_Set_Fn5s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn5s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(LB_User_Data_2) INTENT_IN data
end function LB_fw_Set_Fn5s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn6f"
function LB_fw_Set_Fn6f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn6f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
type(LB_User_Data_3) INTENT_IN data
end function LB_fw_Set_Fn6f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn6s"
function LB_fw_Set_Fn6s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn6s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(LB_User_Data_3) INTENT_IN data
end function LB_fw_Set_Fn6s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn7f"
function LB_fw_Set_Fn7f(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn7f
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(LB_INT), external :: fn_ptr
#endif
type(LB_User_Data_4) INTENT_IN data
end function LB_fw_Set_Fn7f
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_fn7s"
function LB_fw_Set_Fn7s(lb,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(LB_INT) :: LB_fw_Set_Fn7s
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(LB_User_Data_4) INTENT_IN data
end function LB_fw_Set_Fn7s
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_method"
function LB_fw_Set_Method(lb,nbytes,string,length)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Set_Method
integer(LB_INT), dimension(*) INTENT_IN lb, string
integer(LB_INT) INTENT_IN nbytes, length
end function LB_fw_Set_Method
end interface

interface
!NAS$ ALIEN "F77 lb_fw_set_param"
function LB_fw_Set_Param(lb,nbytes,param_name,param_name_len, &
                          new_value,new_value_len)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Set_Param
integer(LB_INT), dimension(*) INTENT_IN lb, param_name, new_value
integer(LB_INT) INTENT_IN nbytes, param_name_len, new_value_len
end function LB_fw_Set_Param
end interface

interface
!NAS$ ALIEN "F77 lb_fw_balance"
function LB_fw_Balance(lb,nbytes,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Balance
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(out) :: changes
integer(LB_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(LB_INT), intent(out) :: num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
end function LB_fw_Balance
end interface


interface
!NAS$ ALIEN "F77 lb_fw_eval"
function LB_fw_Eval(lb,nbytes,print_stats,nobj,obj_wgt, &
                      ncuts,cut_wgt,nboundary,nadj,is_nobj, &
                      is_obj_wgt,is_ncuts,is_cut_wgt,is_nboundary,is_nadj)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Eval
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, print_stats
integer(LB_INT), intent(out) :: nobj, ncuts, nboundary, nadj
real(LB_FLOAT), intent(out) :: obj_wgt(*), cut_wgt(*)
integer(LB_INT), intent(in) :: is_nobj, is_ncuts, is_cut_wgt, is_nboundary, &
                               is_nadj, is_obj_wgt
end function LB_fw_Eval
end interface

interface
!NAS$ ALIEN "F77 lb_fw_point_assign"
function LB_fw_Point_Assign(lb,nbytes,coords,proc)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Point_Assign
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
real(LB_DOUBLE), dimension(*) INTENT_IN coords
integer(LB_INT), intent(out) :: proc
end function LB_fw_Point_Assign
end interface

interface
!NAS$ ALIEN "F77 lb_fw_box_assign"
function LB_fw_Box_Assign(lb,nbytes,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Box_Assign
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
real(LB_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
integer(LB_INT), dimension(*), intent(out) :: procs
integer(LB_INT), intent(out) :: numprocs
end function LB_fw_Box_Assign
end interface

interface
!NAS$ ALIEN "F77 lb_fw_compute_destinations"
function LB_fw_Compute_Destinations(lb,nbytes, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Compute_Destinations
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
integer(LB_INT), dimension(*) INTENT_IN import_global_ids
integer(LB_INT), pointer, dimension(:) :: export_global_ids
integer(LB_INT), dimension(*) INTENT_IN import_local_ids
integer(LB_INT), pointer, dimension(:) :: export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs
integer(LB_INT), pointer, dimension(:) :: export_procs
end function LB_fw_Compute_Destinations
end interface

interface
!NAS$ ALIEN "F77 lb_fw_help_migrate"
function LB_fw_Help_Migrate(lb,nbytes, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Help_Migrate
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), dimension(*) INTENT_IN import_global_ids, export_global_ids
integer(LB_INT), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs, export_procs
end function LB_fw_Help_Migrate
end interface

interface
!NAS$ ALIEN "F77 lb_fw_register_fort_malloc"
subroutine LB_fw_Register_Fort_Malloc(malloc_int,free_int)
use zoltan_types
use lb_user_const
implicit none
#ifdef NASOFTWARE
type(address), intent(in) :: malloc_int, free_int
#else
external malloc_int,free_int
#endif
end subroutine LB_fw_Register_Fort_Malloc
end interface

interface
!NAS$ ALIEN "F77 lb_fw_get_wgt_dim"
function LB_fw_Get_Wgt_Dim(lb,nbytes)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Get_Wgt_Dim
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
end function LB_fw_Get_Wgt_Dim
end interface

interface
!NAS$ ALIEN "F77 lb_fw_get_comm_dim"
function LB_fw_Get_Comm_Dim(lb,nbytes)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Get_Comm_Dim
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
end function LB_fw_Get_Comm_Dim
end interface

! TEMP child_order
interface
!NAS$ ALIEN "F77 lb_fw_get_child_order"
subroutine LB_fw_Get_Child_Order(lb,nbytes,order,ierr)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(inout), dimension(*) :: order
integer(LB_INT), intent(out) :: ierr
end subroutine LB_fw_Get_Child_Order
end interface
! end TEMP child_order

!--------------------------------------------------------------------------
! generic names for the Fortran wrapper procedures

interface Zoltan_Initialize
   module procedure f90LB_Initialize
   module procedure f90LB_Initialize1
end interface

interface Zoltan_Create
   module procedure f90LB_Create
end interface

interface Zoltan_Destroy
   module procedure f90LB_Destroy
end interface

interface Zoltan_Memory_Stats
   module procedure f90LB_Memory_Stats
end interface

interface Zoltan_Set_Fn
   module procedure f90LB_Set_Fn0f
   module procedure f90LB_Set_Fn1f
   module procedure f90LB_Set_Fn2f
   module procedure f90LB_Set_Fn3f
   module procedure f90LB_Set_Fn4f
   module procedure f90LB_Set_Fn5f
   module procedure f90LB_Set_Fn6f
   module procedure f90LB_Set_Fn7f
   module procedure f90LB_Set_Fn0s
   module procedure f90LB_Set_Fn1s
   module procedure f90LB_Set_Fn2s
   module procedure f90LB_Set_Fn3s
   module procedure f90LB_Set_Fn4s
   module procedure f90LB_Set_Fn5s
   module procedure f90LB_Set_Fn6s
   module procedure f90LB_Set_Fn7s
end interface

interface Zoltan_LB_Set_Method
   module procedure f90LB_Set_Method
end interface

interface Zoltan_Set_Param
   module procedure f90LB_Set_Param
end interface

interface Zoltan_LB_Balance
   module procedure f90LB_Balance
end interface

interface Zoltan_LB_Eval
   module procedure f90LB_Eval
end interface

interface Zoltan_LB_Free_Data
   module procedure f90LB_Free_Data
end interface

interface Zoltan_LB_Point_Assign
   module procedure f90LB_Point_Assign
end interface

interface Zoltan_LB_Box_Assign
   module procedure f90LB_Box_Assign
end interface

interface Zoltan_Compute_Destinations
   module procedure f90LB_Compute_Destinations
end interface

interface Zoltan_Help_Migrate
   module procedure f90LB_Help_Migrate
end interface

! TEMP child_order
interface Zoltan_Get_Child_Order
   module procedure f90LB_Get_Child_Order
end interface

#include "set_numgeom.if"
#include "set_geom.if"
#include "set_numedges.if"
#include "set_edgelist.if"
#include "set_numobj.if"
#include "set_objlist.if"
#include "set_firstobj.if"
#include "set_nextobj.if"
#include "set_numborderobj.if"
#include "set_borderobjlist.if"
#include "set_firstborderobj.if"
#include "set_nextborderobj.if"
#include "set_premigrate.if"
#include "set_midmigrate.if"
#include "set_postmigrate.if"
#include "set_objsize.if"
#include "set_packobj.if"
#include "set_unpackobj.if"
#include "set_numcoarseobj.if"
#include "set_coarseobjlist.if"
#include "set_firstcoarseobj.if"
#include "set_nextcoarseobj.if"
#include "set_numchild.if"
#include "set_childlist.if"
#include "set_childweight.if"

!-------------------------------------------------------------------------
! Include LB_* interface for backward compatibility.

#include "lbif.fpp"
#include "set_numgeom.if.lbif"
#include "set_geom.if.lbif"
#include "set_numedges.if.lbif"
#include "set_edgelist.if.lbif"
#include "set_numobj.if.lbif"
#include "set_objlist.if.lbif"
#include "set_firstobj.if.lbif"
#include "set_nextobj.if.lbif"
#include "set_numborderobj.if.lbif"
#include "set_borderobjlist.if.lbif"
#include "set_firstborderobj.if.lbif"
#include "set_nextborderobj.if.lbif"
#include "set_premigrate.if.lbif"
#include "set_midmigrate.if.lbif"
#include "set_postmigrate.if.lbif"
#include "set_objsize.if.lbif"
#include "set_packobj.if.lbif"
#include "set_unpackobj.if.lbif"
#include "set_numcoarseobj.if.lbif"
#include "set_coarseobjlist.if.lbif"
#include "set_firstcoarseobj.if.lbif"
#include "set_nextcoarseobj.if.lbif"
#include "set_numchild.if.lbif"
#include "set_childlist.if.lbif"
#include "set_childweight.if.lbif"


contains

!--------------------------------------------------------------------------
! Utilities

subroutine fort_malloc_int(array,n,ret_addr)
! This gets called by the C special_malloc to do the allocation
integer(LB_INT), pointer :: array(:)
integer(LB_INT), intent(in) :: n
integer(LB_INT_PTR), intent(out) :: ret_addr
integer :: stat
! Allocate the space
allocate(array(n),stat=stat)
if (stat==0) then
! Send the address of the allocated space to C
   call LB_fw_Get_Address_int(array(1),ret_addr)
else
   write(stderr,*) "Error: out of memory during allocation from Fortran"
   ret_addr = 0
endif
end subroutine fort_malloc_int

subroutine fort_free_int(array)
! This gets called by the C special_free to do the deallocation
integer(LB_INT), pointer :: array(:)
integer :: stat
deallocate(array,stat=stat)
if (stat /= 0) then
   write(stderr,*) "Warning: failed to deallocate memory from Fortran"
endif
end subroutine fort_free_int

!--------------------------------------------------------------------------
! Fortran wrapper procedures

function f90LB_Initialize(ver)
integer(LB_INT) :: f90LB_Initialize
real(LB_FLOAT), intent(out) :: ver
#ifdef NASOFTWARE
call LB_fw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_free_int))
#else
call LB_fw_Register_Fort_Malloc(fort_malloc_int,fort_free_int)
#endif
f90LB_Initialize = LB_fw_Initialize(ver)
end function f90LB_Initialize

function f90LB_Initialize1(argc,argv,ver)
integer(LB_INT) :: f90LB_Initialize1
integer(LB_INT) INTENT_IN argc
character(len=*), dimension(*) INTENT_IN argv
real(LB_FLOAT), intent(out) :: ver
integer(LB_INT), allocatable, dimension(:) :: int_argv,starts
integer(LB_INT) :: i, j, leng
#ifdef NASOFTWARE
call LB_fw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_free_int))
#else
call LB_fw_Register_Fort_Malloc(fort_malloc_int,fort_free_int)
#endif
allocate(starts(argc+1), int_argv(len(argv)*argc))
starts(1) = 1
do i=1,argc
   leng = len_trim(argv(i))
   do j=1,leng
      int_argv(j+starts(i)-1) = ichar(argv(i)(j:j))
   end do
   starts(i+1) = starts(i) + leng
end do
f90LB_Initialize1 = LB_fw_Initialize1(argc,int_argv,starts,ver)
deallocate(starts,int_argv)
end function f90LB_Initialize1

function f90LB_Create(communicator)
type(LB_Struct), pointer :: f90LB_Create
integer INTENT_IN communicator
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb
integer(LB_INT) :: nbytes
integer :: i
logical :: isnull
allocate(f90LB_Create)
nbytes = LB_PTR_LENGTH
call LB_fw_Create(communicator,lb,nbytes)
do i=1,LB_PTR_LENGTH
   f90LB_Create%addr%addr(i:i) = char(lb(i))
end do
isnull = (f90LB_Create%addr == LB_NULL_PTR)
if (isnull) then
   deallocate(f90LB_Create)
   nullify(f90LB_Create)
endif
end function f90LB_Create

subroutine f90LB_Destroy(lb)
type(LB_Struct), pointer :: lb
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
call LB_fw_Destroy(lb_addr,nbytes)
deallocate(lb)
nullify(lb)
end subroutine f90LB_Destroy

subroutine f90LB_Memory_Stats()
call LB_fw_Memory_Stats()
end subroutine f90LB_Memory_Stats

function f90LB_Set_Fn0f(lb,fn_type,fn_ptr)
integer(LB_INT) :: f90LB_Set_Fn0f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn0f = LB_fw_Set_Fn0f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
f90LB_Set_Fn0f = LB_fw_Set_Fn0f(lb_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function f90LB_Set_Fn0f

function f90LB_Set_Fn0s(lb,fn_type,fn_ptr)
integer(LB_INT) :: f90LB_Set_Fn0s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn0s = LB_fw_Set_Fn0s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
f90LB_Set_Fn0s = LB_fw_Set_Fn0s(lb_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function f90LB_Set_Fn0s

function f90LB_Set_Fn1f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn1f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
integer(LB_INT) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn1f = LB_fw_Set_Fn1f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn1f = LB_fw_Set_Fn1f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn1f

function f90LB_Set_Fn1s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn1s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(LB_INT) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn1s = LB_fw_Set_Fn1s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn1s = LB_fw_Set_Fn1s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn1s

function f90LB_Set_Fn2f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn2f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
real(LB_FLOAT) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn2f = LB_fw_Set_Fn2f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn2f = LB_fw_Set_Fn2f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn2f

function f90LB_Set_Fn2s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn2s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(LB_FLOAT) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn2s = LB_fw_Set_Fn2s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn2s = LB_fw_Set_Fn2s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn2s

function f90LB_Set_Fn3f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn3f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
real(LB_DOUBLE) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn3f = LB_fw_Set_Fn3f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn3f = LB_fw_Set_Fn3f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn3f

function f90LB_Set_Fn3s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn3s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(LB_DOUBLE) INTENT_IN data(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn3s = LB_fw_Set_Fn3s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn3s = LB_fw_Set_Fn3s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn3s

function f90LB_Set_Fn4f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn4f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
type(LB_User_Data_1) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn4f = LB_fw_Set_Fn4f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn4f = LB_fw_Set_Fn4f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn4f

function f90LB_Set_Fn4s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn4s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_1) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn4s = LB_fw_Set_Fn4s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn4s = LB_fw_Set_Fn4s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn4s

function f90LB_Set_Fn5f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn5f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
type(LB_User_Data_2) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn5f = LB_fw_Set_Fn5f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn5f = LB_fw_Set_Fn5f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn5f

function f90LB_Set_Fn5s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn5s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_2) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn5s = LB_fw_Set_Fn5s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn5s = LB_fw_Set_Fn5s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn5s

function f90LB_Set_Fn6f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn6f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
type(LB_User_Data_3) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn6f = LB_fw_Set_Fn6f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn6f = LB_fw_Set_Fn6f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn6f

function f90LB_Set_Fn6s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn6s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_3) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn6s = LB_fw_Set_Fn6s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn6s = LB_fw_Set_Fn6s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn6s

function f90LB_Set_Fn7f(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn7f
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(LB_INT), external :: fn_ptr
type(LB_User_Data_4) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn7f = LB_fw_Set_Fn7f(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn7f = LB_fw_Set_Fn7f(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn7f

function f90LB_Set_Fn7s(lb,fn_type,fn_ptr,data)
integer(LB_INT) :: f90LB_Set_Fn7s
type(LB_Struct) INTENT_IN lb
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(LB_User_Data_4) INTENT_IN data
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
#ifdef NASOFTWARE
f90LB_Set_Fn7s = LB_fw_Set_Fn7s(lb_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
f90LB_Set_Fn7s = LB_fw_Set_Fn7s(lb_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function f90LB_Set_Fn7s

function f90LB_Set_Method(lb,string)
integer(LB_INT) :: f90LB_Set_Method
type(LB_Struct) INTENT_IN lb
character(len=*) INTENT_IN string
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT), dimension(len_trim(string)) :: int_string
integer(LB_INT) :: nbytes, length, i
nbytes = LB_PTR_LENGTH
length = len_trim(string)
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
do i=1,length
   int_string(i) = ichar(string(i:i))
end do
f90LB_Set_Method = LB_fw_Set_Method(lb_addr,nbytes,int_string,length)
end function f90LB_Set_Method

function f90LB_Set_Param(lb,param_name,new_value)
integer(LB_INT) :: f90LB_Set_Param
type(LB_Struct) INTENT_IN lb
character(len=*) INTENT_IN param_name, new_value
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT), dimension(len_trim(param_name)) :: int_param_name
integer(LB_INT), dimension(len_trim(new_value)) :: int_new_value
integer(LB_INT) :: nbytes, param_name_len, new_value_len, i
nbytes = LB_PTR_LENGTH
param_name_len = len_trim(param_name)
new_value_len = len_trim(new_value)
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
do i=1,param_name_len
   int_param_name(i) = ichar(param_name(i:i))
end do
do i=1,new_value_len
   int_new_value(i) = ichar(new_value(i:i))
end do
f90LB_Set_Param = LB_fw_Set_Param(lb_addr,nbytes,int_param_name, &
                                    param_name_len,int_new_value,new_value_len)
end function f90LB_Set_Param

function f90LB_Balance(lb,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Balance
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
integer(LB_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(LB_INT), intent(out) :: num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_changes
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Balance = LB_fw_Balance(lb_addr,nbytes,int_changes, &
                             num_gid_entries, num_lid_entries, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function f90LB_Balance

function f90LB_Eval(lb,print_stats,nobj,obj_wgt, &
                    ncuts,cut_wgt,nboundary,nadj)
integer(LB_INT) :: f90LB_Eval
type(LB_Struct) INTENT_IN lb
logical INTENT_IN print_stats
integer(LB_INT), intent(out), optional :: nobj, ncuts, nboundary, nadj
real(LB_FLOAT), intent(out), optional :: obj_wgt(*), cut_wgt(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_print_stats, dim, edim
integer(LB_INT) :: loc_nobj, loc_ncuts, loc_nboundary, loc_nadj
real(LB_FLOAT), allocatable :: loc_obj_wgt(:), loc_cut_wgt(:)
integer(LB_INT) :: is_nobj, is_ncuts, is_cut_wgt, is_nboundary, is_nadj, is_obj_wgt
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
if (print_stats) then
   int_print_stats = 1
else
   int_print_stats = 0
endif
if (present(nobj)) then
   is_nobj = 1
else
   is_nobj = 0
endif
if (present(ncuts)) then
   is_ncuts = 1
else
   is_ncuts = 0
endif
if (present(nboundary)) then
   is_nboundary = 1
else
   is_nboundary = 0
endif
if (present(nadj)) then
   is_nadj = 1
else
   is_nadj = 0
endif
if (present(obj_wgt)) then
   is_obj_wgt = 1
   dim = LB_fw_Get_Wgt_Dim(lb_addr,nbytes)
   allocate(loc_obj_wgt(dim))
else
   is_obj_wgt = 0
   allocate(loc_obj_wgt(1))
endif
if (present(cut_wgt)) then
   is_cut_wgt = 1
   edim = LB_fw_Get_Comm_Dim(lb_addr,nbytes)
   allocate(loc_cut_wgt(edim))
else
   is_cut_wgt = 0
   allocate(loc_cut_wgt(1))
endif
f90LB_Eval = LB_fw_Eval(lb_addr,nbytes,int_print_stats,loc_nobj,loc_obj_wgt, &
                loc_ncuts,loc_cut_wgt,loc_nboundary,loc_nadj,is_nobj, &
                is_obj_wgt,is_ncuts,is_cut_wgt,is_nboundary,is_nadj)
if (present(nobj)) nobj = loc_nobj
if (present(obj_wgt)) then
   do i = 1,dim
      obj_wgt(i) = loc_obj_wgt(i)
   end do
endif
if (present(cut_wgt)) then
   do i = 1,edim
      cut_wgt(i) = loc_cut_wgt(i)
   end do
endif
if (present(ncuts)) ncuts = loc_ncuts
if (present(nboundary)) nboundary = loc_nboundary
if (present(nadj)) nadj = loc_nadj
deallocate(loc_obj_wgt)
end function f90LB_Eval

function f90LB_Free_Data(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Free_Data
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
f90LB_Free_Data = ZOLTAN_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data = ZOLTAN_WARN
nullify(export_procs)
end function f90LB_Free_Data

function f90LB_Point_Assign(lb,coords,proc)
integer(LB_INT) :: f90LB_Point_Assign
type(LB_Struct) INTENT_IN lb
real(LB_DOUBLE), dimension(*) INTENT_IN coords
integer(LB_INT), intent(out) :: proc
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Point_Assign = LB_fw_Point_Assign(lb_addr,nbytes,coords,proc)
end function f90LB_Point_Assign

function f90LB_Box_Assign(lb,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
integer(LB_INT) :: f90LB_Box_Assign
type(LB_Struct) INTENT_IN lb
real(LB_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
integer(LB_INT), intent(out), dimension(*) :: procs
integer(LB_INT), intent(out) :: numprocs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Box_Assign = LB_fw_Box_Assign(lb_addr,nbytes,xmin,ymin,zmin,xmax,ymax, &
                                    zmax,procs,numprocs)
end function f90LB_Box_Assign

function f90LB_Compute_Destinations(lb, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Compute_Destinations
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs)) then
   write(stderr,*) "Error from LB_Compute_Destinations: import pointers are not associated"
   f90LB_Compute_Destinations = ZOLTAN_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Compute_Destinations = LB_fw_Compute_Destinations(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Compute_Destinations


function f90LB_Help_Migrate(lb, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Help_Migrate
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
logical :: free_import_global_ids, free_import_local_ids, free_import_procs
logical :: free_export_global_ids, free_export_local_ids, free_export_procs

if ((num_import.gt.0).and.(.not.associated(import_global_ids) .or. &
                           .not.associated(import_local_ids)  .or. &
                           .not.associated(import_procs))) then
   write(stderr,*) "Error from LB_Help_Migrate: import pointers are not associated"
   f90LB_Help_Migrate = ZOLTAN_WARN
   return
endif
if ((num_export.gt.0).and.(.not.associated(export_procs) .or. &
                           .not.associated(export_global_ids) .or. &
                           .not.associated(export_local_ids))) then
   write(stderr,*) "Error from LB_Help_Migrate: export pointers are not associated"

   f90LB_Help_Migrate = ZOLTAN_WARN
   return
endif

! generate place-holders to make call to LB_fw_Help_Migrate valid;
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

nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Help_Migrate = LB_fw_Help_Migrate(lb_addr,nbytes, &
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

end function f90LB_Help_Migrate

! TEMP child_order
subroutine f90LB_Get_Child_Order(lb,order,ierr)
type(LB_Struct), pointer :: lb
integer(LB_INT), intent(inout), dimension(*) :: order
integer(LB_INT), intent(out) :: ierr
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
call LB_fw_Get_Child_Order(lb_addr,nbytes,order,ierr)
end subroutine f90LB_Get_Child_Order
! end TEMP child_order


#include "set_numgeom.fn"
#include "set_geom.fn"
#include "set_numedges.fn"
#include "set_edgelist.fn"
#include "set_numobj.fn"
#include "set_objlist.fn"
#include "set_firstobj.fn"
#include "set_nextobj.fn"
#include "set_numborderobj.fn"
#include "set_borderobjlist.fn"
#include "set_firstborderobj.fn"
#include "set_nextborderobj.fn"
#include "set_premigrate.fn"
#include "set_midmigrate.fn"
#include "set_postmigrate.fn"
#include "set_objsize.fn"
#include "set_packobj.fn"
#include "set_unpackobj.fn"
#include "set_numcoarseobj.fn"
#include "set_coarseobjlist.fn"
#include "set_firstcoarseobj.fn"
#include "set_nextcoarseobj.fn"
#include "set_numchild.fn"
#include "set_childlist.fn"
#include "set_childweight.fn"

end module zoltan
