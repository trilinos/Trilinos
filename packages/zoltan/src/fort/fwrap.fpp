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
   LB_User_Data_4, &
   LB_GID, &
   LB_LID

public :: &
   LB_Struct, &
   LB_FN_TYPEF, &
   LB_FN_TYPES

public :: &
   LB_NUM_EDGES_FN_TYPE, &
   LB_EDGE_LIST_FN_TYPE, &
   LB_NUM_GEOM_FN_TYPE, &
   LB_GEOM_FN_TYPE, &
   LB_NUM_OBJ_FN_TYPE, &
   LB_OBJ_LIST_FN_TYPE, &
   LB_FIRST_OBJ_FN_TYPE, &
   LB_NEXT_OBJ_FN_TYPE, &
   LB_NUM_BORDER_OBJ_FN_TYPE, &
   LB_BORDER_OBJ_LIST_FN_TYPE, &
   LB_FIRST_BORDER_OBJ_FN_TYPE, &
   LB_NEXT_BORDER_OBJ_FN_TYPE, &
   LB_PRE_MIGRATE_FN_TYPE, &
   LB_POST_MIGRATE_FN_TYPE, &
   LB_OBJ_SIZE_FN_TYPE, &
   LB_PACK_OBJ_FN_TYPE, &
   LB_UNPACK_OBJ_FN_TYPE, &
   LB_NUM_COARSE_OBJ_FN_TYPE, &
   LB_COARSE_OBJ_LIST_FN_TYPE, &
   LB_FIRST_COARSE_OBJ_FN_TYPE, &
   LB_NEXT_COARSE_OBJ_FN_TYPE, &
   LB_NUM_CHILD_FN_TYPE, &
   LB_CHILD_LIST_FN_TYPE, &
   LB_CHILD_WEIGHT_FN_TYPE

public :: &
   LB_OTHER_REF, &
   LB_IN_ORDER, &
   LB_TRI_BISECT, &
   LB_QUAD_QUAD, &
   LB_HEX3D_OCT

public :: &
   LB_OK, &
   LB_WARN, &
   LB_FATAL, &
   LB_MEMERR

public :: &
   LB_Initialize, &
   LB_Create, &
   LB_Destroy, &
   LB_Set_Fn, &
   LB_Set_Method, &
   LB_Set_Param, &
   LB_Balance, &
   LB_Eval, &
   LB_Free_Data, &
   LB_Point_Assign, &
   LB_Box_Assign, &
   LB_Compute_Destinations, &
   LB_Help_Migrate

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

type LB_FN_TYPEF
   private
   integer(LB_INT) :: choice
end type LB_FN_TYPEF

type LB_FN_TYPES
   private
   integer(LB_INT) :: choice
end type LB_FN_TYPES

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(LB_FN_TYPEF) :: &
#else
type(LB_FN_TYPEF), parameter :: &
#endif
   LB_NUM_EDGES_FN_TYPE        = LB_FN_TYPEF(0_LB_INT), &
   LB_NUM_GEOM_FN_TYPE         = LB_FN_TYPEF(2_LB_INT), &
   LB_NUM_OBJ_FN_TYPE          = LB_FN_TYPEF(4_LB_INT), &
   LB_FIRST_OBJ_FN_TYPE        = LB_FN_TYPEF(6_LB_INT), &
   LB_NEXT_OBJ_FN_TYPE         = LB_FN_TYPEF(7_LB_INT), &
   LB_NUM_BORDER_OBJ_FN_TYPE   = LB_FN_TYPEF(8_LB_INT), &
   LB_FIRST_BORDER_OBJ_FN_TYPE = LB_FN_TYPEF(10_LB_INT), &
   LB_NEXT_BORDER_OBJ_FN_TYPE  = LB_FN_TYPEF(11_LB_INT), &
   LB_OBJ_SIZE_FN_TYPE         = LB_FN_TYPEF(14_LB_INT), &
   LB_NUM_COARSE_OBJ_FN_TYPE   = LB_FN_TYPEF(17_LB_INT), &
   LB_FIRST_COARSE_OBJ_FN_TYPE = LB_FN_TYPEF(19_LB_INT), &
   LB_NEXT_COARSE_OBJ_FN_TYPE  = LB_FN_TYPEF(20_LB_INT), &
   LB_NUM_CHILD_FN_TYPE        = LB_FN_TYPEF(21_LB_INT)

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(LB_FN_TYPES) :: &
#else
type(LB_FN_TYPES), parameter :: &
#endif
   LB_EDGE_LIST_FN_TYPE        = LB_FN_TYPES(1_LB_INT), &
   LB_GEOM_FN_TYPE             = LB_FN_TYPES(3_LB_INT), &
   LB_OBJ_LIST_FN_TYPE         = LB_FN_TYPES(5_LB_INT), &
   LB_BORDER_OBJ_LIST_FN_TYPE  = LB_FN_TYPES(9_LB_INT), &
   LB_PRE_MIGRATE_FN_TYPE      = LB_FN_TYPES(12_LB_INT), &
   LB_POST_MIGRATE_FN_TYPE     = LB_FN_TYPES(13_LB_INT), &
   LB_PACK_OBJ_FN_TYPE         = LB_FN_TYPES(15_LB_INT), &
   LB_UNPACK_OBJ_FN_TYPE       = LB_FN_TYPES(16_LB_INT), &
   LB_COARSE_OBJ_LIST_FN_TYPE  = LB_FN_TYPES(18_LB_INT), &
   LB_CHILD_LIST_FN_TYPE       = LB_FN_TYPES(22_LB_INT), &
   LB_CHILD_WEIGHT_FN_TYPE     = LB_FN_TYPES(23_LB_INT)

! Type of refinement used when building a refinement tree
! These values must agree with the values in lb/lbi_const.h

integer(LB_INT), parameter :: &
  LB_OTHER_REF     = 0_LB_INT, &
  LB_IN_ORDER      = 1_LB_INT, &
  LB_TRI_BISECT    = 2_LB_INT, &
  LB_QUAD_QUAD     = 3_LB_INT, &
  LB_HEX3D_OCT     = 4_LB_INT

! Error codes for LB library
! These values must agree with the values in lb/lbi_const.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(LB_INT) :: &
#else
integer(LB_INT), parameter :: &
#endif
   LB_OK     =  0_LB_INT, &
   LB_WARN   =  1_LB_INT, &
   LB_FATAL  = -1_LB_INT, &
   LB_MEMERR = -2_LB_INT

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
integer(LB_INT), intent(out) :: ret_addr
end subroutine LB_fw_Get_Address_int
end interface

interface
!NAS$ ALIEN "F77 lb_fw_get_address_gid"
subroutine LB_fw_Get_Address_GID(arg,ret_addr)
use zoltan_types
use lb_user_const
type(LB_GID) :: arg
integer(LB_INT), intent(out) :: ret_addr
end subroutine LB_fw_Get_Address_GID
end interface

interface
!NAS$ ALIEN "F77 lb_fw_get_address_lid"
subroutine LB_fw_Get_Address_LID(arg,ret_addr)
use zoltan_types
use lb_user_const
type(LB_LID) :: arg
integer(LB_INT), intent(out) :: ret_addr
end subroutine LB_fw_Get_Address_LID
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
!NAS$ ALIEN "F77 lb_fw_balance11"
function LB_fw_Balance11(lb,nbytes,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Balance11
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
end function LB_fw_Balance11
end interface

interface
!NAS$ ALIEN "F77 lb_fw_balance12"
function LB_fw_Balance12(lb,nbytes,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Balance12
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
end function LB_fw_Balance12
end interface

interface
!NAS$ ALIEN "F77 lb_fw_balance21"
function LB_fw_Balance21(lb,nbytes,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Balance21
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
end function LB_fw_Balance21
end interface

interface
!NAS$ ALIEN "F77 lb_fw_balance22"
function LB_fw_Balance22(lb,nbytes,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Balance22
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT), intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
end function LB_fw_Balance22
end interface

interface
!NAS$ ALIEN "F77 lb_fw_eval"
subroutine LB_fw_Eval(lb,nbytes,print_stats,nobj,obj_wgt, &
                      cut_wgt,nboundary,nadj,ierr,is_nobj, &
                      is_obj_wgt,is_cut_wgt,is_nboundary,is_nadj)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes, print_stats
integer(LB_INT), intent(out) :: nobj, cut_wgt, nboundary, nadj, ierr
real(LB_FLOAT), intent(out) :: obj_wgt(*)
integer(LB_INT), intent(in) :: is_nobj, is_cut_wgt, is_nboundary, is_nadj, &
                               is_obj_wgt
end subroutine LB_fw_Eval
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
!NAS$ ALIEN "F77 lb_fw_compute_destinations11"
function LB_fw_Compute_Destinations11(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Compute_Destinations11
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
type(LB_GID), dimension(*) INTENT_IN import_global_ids
type(LB_GID), pointer, dimension(:) :: export_global_ids
type(LB_LID), dimension(*) INTENT_IN import_local_ids
type(LB_LID), pointer, dimension(:) :: export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs
integer(LB_INT), pointer, dimension(:) :: export_procs
end function LB_fw_Compute_Destinations11
end interface

interface
!NAS$ ALIEN "F77 lb_fw_compute_destinations12"
function LB_fw_Compute_Destinations12(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Compute_Destinations12
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
type(LB_GID), dimension(*) INTENT_IN import_global_ids
type(LB_GID), pointer, dimension(:) :: export_global_ids
integer(LB_INT), dimension(*) INTENT_IN import_local_ids
integer(LB_INT), pointer, dimension(:) :: export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs
integer(LB_INT), pointer, dimension(:) :: export_procs
end function LB_fw_Compute_Destinations12
end interface

interface
!NAS$ ALIEN "F77 lb_fw_compute_destinations21"
function LB_fw_Compute_Destinations21(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Compute_Destinations21
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
integer(LB_INT), dimension(*) INTENT_IN import_global_ids
integer(LB_INT), pointer, dimension(:) :: export_global_ids
type(LB_LID), dimension(*) INTENT_IN import_local_ids
type(LB_LID), pointer, dimension(:) :: export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs
integer(LB_INT), pointer, dimension(:) :: export_procs
end function LB_fw_Compute_Destinations21
end interface

interface
!NAS$ ALIEN "F77 lb_fw_compute_destinations22"
function LB_fw_Compute_Destinations22(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Compute_Destinations22
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
end function LB_fw_Compute_Destinations22
end interface


interface
!NAS$ ALIEN "F77 lb_fw_help_migrate11"
function LB_fw_Help_Migrate11(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Help_Migrate11
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import, num_export
type(LB_GID), dimension(*) INTENT_IN import_global_ids, export_global_ids
type(LB_LID), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs, export_procs
end function LB_fw_Help_Migrate11
end interface

interface
!NAS$ ALIEN "F77 lb_fw_help_migrate12"
function LB_fw_Help_Migrate12(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Help_Migrate12
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import, num_export
type(LB_GID), dimension(*) INTENT_IN import_global_ids, export_global_ids
integer(LB_INT), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs, export_procs
end function LB_fw_Help_Migrate12
end interface

interface
!NAS$ ALIEN "F77 lb_fw_help_migrate21"
function LB_fw_Help_Migrate21(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Help_Migrate21
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), dimension(*) INTENT_IN import_global_ids, export_global_ids
type(LB_LID), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs, export_procs
end function LB_fw_Help_Migrate21
end interface

interface
!NAS$ ALIEN "F77 lb_fw_help_migrate22"
function LB_fw_Help_Migrate22(lb,nbytes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
implicit none
integer(LB_INT) :: LB_fw_Help_Migrate22
integer(LB_INT), dimension(*) INTENT_IN lb
integer(LB_INT) INTENT_IN nbytes
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), dimension(*) INTENT_IN import_global_ids, export_global_ids
integer(LB_INT), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(LB_INT), dimension(*) INTENT_IN import_procs, export_procs
end function LB_fw_Help_Migrate22
end interface

interface
!NAS$ ALIEN "F77 lb_fw_register_fort_malloc"
subroutine LB_fw_Register_Fort_Malloc(malloc_int,malloc_gid,malloc_lid, &
                                      free_int,free_gid,free_lid)
use zoltan_types
use lb_user_const
implicit none
#ifdef NASOFTWARE
type(address), intent(in) :: malloc_int,malloc_gid,malloc_lid, &
                             free_int,free_gid,free_lid
#else
external malloc_int,malloc_gid,malloc_lid,free_int,free_gid,free_lid
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

!--------------------------------------------------------------------------
! generic names for the Fortran wrapper procedures

interface LB_Initialize
   module procedure f90LB_Initialize
   module procedure f90LB_Initialize1
end interface

interface LB_Create
   module procedure f90LB_Create
end interface

interface LB_Destroy
   module procedure f90LB_Destroy
end interface

interface LB_Set_Fn
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

interface LB_Set_Method
   module procedure f90LB_Set_Method
end interface

interface LB_Set_Param
   module procedure f90LB_Set_Param
end interface

interface LB_Balance
   module procedure f90LB_Balance11
   module procedure f90LB_Balance12
   module procedure f90LB_Balance21
   module procedure f90LB_Balance22
end interface

interface LB_Eval
   module procedure f90LB_Eval
end interface

interface LB_Free_Data
   module procedure f90LB_Free_Data11
   module procedure f90LB_Free_Data12
   module procedure f90LB_Free_Data21
   module procedure f90LB_Free_Data22
end interface

interface LB_Point_Assign
   module procedure f90LB_Point_Assign
end interface

interface LB_Box_Assign
   module procedure f90LB_Box_Assign
end interface

interface LB_Compute_Destinations
   module procedure f90LB_Compute_Destinations11
   module procedure f90LB_Compute_Destinations12
   module procedure f90LB_Compute_Destinations21
   module procedure f90LB_Compute_Destinations22
end interface

interface LB_Help_Migrate
   module procedure f90LB_Help_Migrate11
   module procedure f90LB_Help_Migrate12
   module procedure f90LB_Help_Migrate21
   module procedure f90LB_Help_Migrate22
end interface

contains

!--------------------------------------------------------------------------
! Utilities

subroutine fort_malloc_int(array,n,ret_addr)
! This gets called by the C special_malloc to do the allocation
integer(LB_INT), pointer :: array(:)
integer(LB_INT), intent(in) :: n
integer(LB_INT), intent(out) :: ret_addr
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

subroutine fort_malloc_gid(array,n,ret_addr)
! This gets called by the C special_malloc to do the allocation
type(LB_GID), pointer :: array(:)
integer(LB_INT), intent(in) :: n
integer(LB_INT), intent(out) :: ret_addr
integer :: stat
! Allocate the space
allocate(array(n),stat=stat)
if (stat==0) then
! Send the address of the allocated space to C
   call LB_fw_Get_Address_GID(array(1),ret_addr)
else
   write(stderr,*) "Error: out of memory during allocation from Fortran"
   ret_addr = 0
endif
end subroutine fort_malloc_gid

subroutine fort_malloc_lid(array,n,ret_addr)
! This gets called by the C special_malloc to do the allocation
type(LB_LID), pointer :: array(:)
integer(LB_INT), intent(in) :: n
integer(LB_INT), intent(out) :: ret_addr
integer :: stat
! Allocate the space
allocate(array(n),stat=stat)
if (stat==0) then
! Send the address of the allocated space to C
   call LB_fw_Get_Address_LID(array(1),ret_addr)
else
   write(stderr,*) "Error: out of memory during allocation from Fortran"
   ret_addr = 0
endif
end subroutine fort_malloc_lid

subroutine fort_free_int(array)
! This gets called by the C special_free to do the deallocation
integer(LB_INT), pointer :: array(:)
integer :: stat
deallocate(array,stat=stat)
if (stat /= 0) then
   write(stderr,*) "Warning: failed to deallocate memory from Fortran"
endif
end subroutine fort_free_int

subroutine fort_free_gid(array)
! This gets called by the C special_free to do the deallocation
type(LB_GID), pointer :: array(:)
integer :: stat
deallocate(array,stat=stat)
if (stat /= 0) then
   write(stderr,*) "Warning: failed to deallocate memory from Fortran"
endif
end subroutine fort_free_gid

subroutine fort_free_lid(array)
! This gets called by the C special_free to do the deallocation
type(LB_LID), pointer :: array(:)
integer :: stat
deallocate(array,stat=stat)
if (stat /= 0) then
   write(stderr,*) "Warning: failed to deallocate memory from Fortran"
endif
end subroutine fort_free_lid

!--------------------------------------------------------------------------
! Fortran wrapper procedures

function f90LB_Initialize(ver)
integer(LB_INT) :: f90LB_Initialize
real(LB_FLOAT), intent(out) :: ver
#ifdef NASOFTWARE
call LB_fw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_malloc_gid), &
                                loc(fort_malloc_lid),loc(fort_free_int), &
                                loc(fort_free_gid),loc(fort_free_lid))
#else
call LB_fw_Register_Fort_Malloc(fort_malloc_int,fort_malloc_gid,fort_malloc_lid,&
                                fort_free_int,fort_free_gid,fort_free_lid)
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
call LB_fw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_malloc_gid), &
                                loc(fort_malloc_lid),loc(fort_free_int), &
                                loc(fort_free_gid),loc(fort_free_lid))
#else
call LB_fw_Register_Fort_Malloc(fort_malloc_int,fort_malloc_gid,fort_malloc_lid,&
                                fort_free_int,fort_free_gid,fort_free_lid)
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

function f90LB_Set_Fn0f(lb,fn_type,fn_ptr)
integer(LB_INT) :: f90LB_Set_Fn0f
type(LB_Struct) INTENT_IN lb
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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
type(LB_FN_TYPEF) INTENT_IN fn_type
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
type(LB_FN_TYPES) INTENT_IN fn_type
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

function f90LB_Balance11(lb,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Balance11
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_changes
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Balance11 = LB_fw_Balance11(lb_addr,nbytes,int_changes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function f90LB_Balance11

function f90LB_Balance12(lb,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Balance12
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_changes
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Balance12 = LB_fw_Balance12(lb_addr,nbytes,int_changes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function f90LB_Balance12

function f90LB_Balance21(lb,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Balance21
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
integer(LB_INT), intent(out) :: num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_changes
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Balance21 = LB_fw_Balance21(lb_addr,nbytes,int_changes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function f90LB_Balance21

function f90LB_Balance22(lb,changes,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Balance22
type(LB_Struct) INTENT_IN lb
logical, intent(out) :: changes
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
f90LB_Balance22 = LB_fw_Balance22(lb_addr,nbytes,int_changes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function f90LB_Balance22

subroutine f90LB_Eval(lb,print_stats,nobj,obj_wgt, &
                      cut_wgt,nboundary,nadj,ierr)
type(LB_Struct) INTENT_IN lb
logical INTENT_IN print_stats
integer(LB_INT), intent(out), optional :: nobj, cut_wgt, nboundary, nadj
integer(LB_INT), intent(out) :: ierr
real(LB_FLOAT), intent(out), optional :: obj_wgt(*)
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i, int_print_stats, dim
integer(LB_INT) :: loc_nobj, loc_cut_wgt, loc_nboundary, loc_nadj
real(LB_FLOAT), allocatable :: loc_obj_wgt(:)
integer(LB_INT) :: is_nobj, is_cut_wgt, is_nboundary, is_nadj, is_obj_wgt
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
if (present(cut_wgt)) then
   is_cut_wgt = 1
else
   is_cut_wgt = 0
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
call LB_fw_Eval(lb_addr,nbytes,int_print_stats,loc_nobj,loc_obj_wgt, &
                loc_cut_wgt,loc_nboundary,loc_nadj,ierr,is_nobj,is_obj_wgt, &
                is_cut_wgt,is_nboundary,is_nadj)
if (present(nobj)) nobj = loc_nobj
if (present(obj_wgt)) then
   do i = 1,dim
      obj_wgt(i) = loc_obj_wgt(i)
   end do
endif
if (present(cut_wgt)) cut_wgt = loc_cut_wgt
if (present(nboundary)) nboundary = loc_nboundary
if (present(nadj)) nadj = loc_nadj
deallocate(loc_obj_wgt)
end subroutine f90LB_Eval

function f90LB_Free_Data11(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Free_Data11
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
f90LB_Free_Data11 = LB_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data11 = LB_WARN
nullify(export_procs)
end function f90LB_Free_Data11

function f90LB_Free_Data12(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Free_Data12
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
f90LB_Free_Data12 = LB_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data12 = LB_WARN
nullify(export_procs)
end function f90LB_Free_Data12

function f90LB_Free_Data21(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Free_Data21
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
f90LB_Free_Data21 = LB_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data21 = LB_WARN
nullify(export_procs)
end function f90LB_Free_Data21

function f90LB_Free_Data22(import_global_ids, import_local_ids,import_procs, &
                         export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Free_Data22
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer :: stat
stat = 0
f90LB_Free_Data22 = LB_OK
if (associated(import_global_ids)) deallocate(import_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(import_global_ids)
if (associated(import_local_ids)) deallocate(import_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(import_local_ids)
if (associated(import_procs)) deallocate(import_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(import_procs)
if (associated(export_global_ids)) deallocate(export_global_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(export_global_ids)
if (associated(export_local_ids)) deallocate(export_local_ids,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(export_local_ids)
if (associated(export_procs)) deallocate(export_procs,stat=stat)
if (stat /= 0) f90LB_Free_Data22 = LB_WARN
nullify(export_procs)
end function f90LB_Free_Data22

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

function f90LB_Compute_Destinations11(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Compute_Destinations11
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs)) then
   write(stderr,*) "Error from LB_Compute_Destinations: import pointers are not associated"
   f90LB_Compute_Destinations11 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Compute_Destinations11 = LB_fw_Compute_Destinations11(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Compute_Destinations11

function f90LB_Compute_Destinations12(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Compute_Destinations12
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs)) then
   write(stderr,*) "Error from LB_Compute_Destinations: import pointers are not associated"
   f90LB_Compute_Destinations12 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Compute_Destinations12 = LB_fw_Compute_Destinations12(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Compute_Destinations12

function f90LB_Compute_Destinations21(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Compute_Destinations21
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import
integer(LB_INT), intent(out) :: num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs)) then
   write(stderr,*) "Error from LB_Compute_Destinations: import pointers are not associated"
   f90LB_Compute_Destinations21 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Compute_Destinations21 = LB_fw_Compute_Destinations21(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Compute_Destinations21

function f90LB_Compute_Destinations22(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Compute_Destinations22
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
   f90LB_Compute_Destinations22 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Compute_Destinations22 = LB_fw_Compute_Destinations22(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Compute_Destinations22

function f90LB_Help_Migrate11(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Help_Migrate11
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs) .or. .not.associated(export_procs) &
    .or. .not.associated(export_global_ids) &
    .or. .not.associated(export_local_ids)) then
   write(stderr,*) "Error from LB_Help_Migrate: import or export pointers are not associated"
   f90LB_Help_Migrate11 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Help_Migrate11 = LB_fw_Help_Migrate11(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Help_Migrate11

function f90LB_Help_Migrate12(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Help_Migrate12
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import, num_export
type(LB_GID), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs) .or. .not.associated(export_procs) &
    .or. .not.associated(export_global_ids) &
    .or. .not.associated(export_local_ids)) then
   write(stderr,*) "Error from LB_Help_Migrate: import or export pointers are not associated"
   f90LB_Help_Migrate12 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Help_Migrate12 = LB_fw_Help_Migrate12(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Help_Migrate12

function f90LB_Help_Migrate21(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Help_Migrate21
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
type(LB_LID), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs) .or. .not.associated(export_procs) &
    .or. .not.associated(export_global_ids) &
    .or. .not.associated(export_local_ids)) then
   write(stderr,*) "Error from LB_Help_Migrate: import or export pointers are not associated"
   f90LB_Help_Migrate21 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Help_Migrate21 = LB_fw_Help_Migrate21(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Help_Migrate21

function f90LB_Help_Migrate22(lb,num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(LB_INT) :: f90LB_Help_Migrate22
type(LB_Struct) INTENT_IN lb
integer(LB_INT) INTENT_IN num_import, num_export
integer(LB_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(LB_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(LB_INT), pointer, dimension(:) :: import_procs, export_procs
integer(LB_INT), dimension(LB_PTR_LENGTH) :: lb_addr
integer(LB_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs) .or. .not.associated(export_procs) &
    .or. .not.associated(export_global_ids) &
    .or. .not.associated(export_local_ids)) then
   write(stderr,*) "Error from LB_Help_Migrate: import or export pointers are not associated"
   f90LB_Help_Migrate22 = LB_WARN
   return
endif
nbytes = LB_PTR_LENGTH
do i=1,nbytes
   lb_addr(i) = ichar(lb%addr%addr(i:i))
end do
f90LB_Help_Migrate22 = LB_fw_Help_Migrate22(lb_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function f90LB_Help_Migrate22

end module zoltan
