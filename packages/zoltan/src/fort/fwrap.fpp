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

#include "sppr_header"
!--------------------------------------------------------------------------
! preprocessor directives to handle special case compilers

#ifdef NASOFTWARE
#define INTENT_IN ::
#define INTENT_OUT ::
#define C_BINDING bind(c)
#define PASS_BY_REF ,pass_by("*")
#else
#define INTENT_IN ,intent(in) ::
#define INTENT_OUT ,intent(out) ::
#define C_BINDING
#define PASS_BY_REF
#endif

module zoltan
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
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
   ZOLTAN_CHILD_WEIGHT_FN_TYPE, &
   ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE, &
   ZOLTAN_PACK_OBJ_MULTI_FN_TYPE, &
   ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE

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
   Zoltan_Set_Obj_Size_Fn, Zoltan_Set_Pack_Obj_Fn, Zoltan_Set_Unpack_Obj_Fn, &
   Zoltan_Set_Obj_Size_Multi_Fn, &
   Zoltan_Set_Pack_Obj_Multi_Fn, Zoltan_Set_Unpack_Obj_Multi_Fn 

public :: &
   Zoltan_Get_Child_Order ! TEMP child_order

!--------------------------------------------------------------------------
! user defined types corresponding to the C structs

type Zoltan_Struct
   private
   sequence
   type(Zoltan_PTR) :: addr
#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
   integer :: dummy
#endif
end type Zoltan_Struct

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

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPEF) :: &
#else
type(ZOLTAN_FN_TYPEF), parameter :: &
#endif
   ZOLTAN_NUM_EDGES_FN_TYPE        = ZOLTAN_FN_TYPEF(0_Zoltan_INT), &
   ZOLTAN_NUM_GEOM_FN_TYPE         = ZOLTAN_FN_TYPEF(2_Zoltan_INT), &
   ZOLTAN_NUM_OBJ_FN_TYPE          = ZOLTAN_FN_TYPEF(4_Zoltan_INT), &
   ZOLTAN_FIRST_OBJ_FN_TYPE        = ZOLTAN_FN_TYPEF(6_Zoltan_INT), &
   ZOLTAN_NEXT_OBJ_FN_TYPE         = ZOLTAN_FN_TYPEF(7_Zoltan_INT), &
   ZOLTAN_NUM_BORDER_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(8_Zoltan_INT), &
   ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(10_Zoltan_INT), &
   ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(11_Zoltan_INT), &
   ZOLTAN_OBJ_SIZE_FN_TYPE         = ZOLTAN_FN_TYPEF(15_Zoltan_INT), &
   ZOLTAN_NUM_COARSE_OBJ_FN_TYPE   = ZOLTAN_FN_TYPEF(18_Zoltan_INT), &
   ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE = ZOLTAN_FN_TYPEF(20_Zoltan_INT), &
   ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE  = ZOLTAN_FN_TYPEF(21_Zoltan_INT), &
   ZOLTAN_NUM_CHILD_FN_TYPE        = ZOLTAN_FN_TYPEF(22_Zoltan_INT)

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPES) :: &
#else
type(ZOLTAN_FN_TYPES), parameter :: &
#endif
   ZOLTAN_EDGE_LIST_FN_TYPE        = ZOLTAN_FN_TYPES(1_Zoltan_INT), &
   ZOLTAN_GEOM_FN_TYPE             = ZOLTAN_FN_TYPES(3_Zoltan_INT), &
   ZOLTAN_OBJ_LIST_FN_TYPE         = ZOLTAN_FN_TYPES(5_Zoltan_INT), &
   ZOLTAN_BORDER_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(9_Zoltan_INT), &
   ZOLTAN_PRE_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(12_Zoltan_INT), &
   ZOLTAN_MID_MIGRATE_FN_TYPE      = ZOLTAN_FN_TYPES(13_Zoltan_INT), &
   ZOLTAN_POST_MIGRATE_FN_TYPE     = ZOLTAN_FN_TYPES(14_Zoltan_INT), &
   ZOLTAN_PACK_OBJ_FN_TYPE         = ZOLTAN_FN_TYPES(16_Zoltan_INT), &
   ZOLTAN_UNPACK_OBJ_FN_TYPE       = ZOLTAN_FN_TYPES(17_Zoltan_INT), &
   ZOLTAN_COARSE_OBJ_LIST_FN_TYPE  = ZOLTAN_FN_TYPES(19_Zoltan_INT), &
   ZOLTAN_CHILD_LIST_FN_TYPE       = ZOLTAN_FN_TYPES(23_Zoltan_INT), &
   ZOLTAN_CHILD_WEIGHT_FN_TYPE     = ZOLTAN_FN_TYPES(24_Zoltan_INT), &
   ZOLTAN_OBJ_SIZE_MULTI_FN_TYPE   = ZOLTAN_FN_TYPES(25_Zoltan_INT), &
   ZOLTAN_PACK_OBJ_MULTI_FN_TYPE   = ZOLTAN_FN_TYPES(26_Zoltan_INT), &
   ZOLTAN_UNPACK_OBJ_MULTI_FN_TYPE = ZOLTAN_FN_TYPES(27_Zoltan_INT)

! Type of refinement used when building a refinement tree
! These values must agree with the values in zoltan.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(Zoltan_INT) :: &
#else
integer(Zoltan_INT), parameter :: &
#endif
  ZOLTAN_OTHER_REF     = 0_Zoltan_INT, &
  ZOLTAN_IN_ORDER      = 1_Zoltan_INT, &
  ZOLTAN_TRI_BISECT    = 2_Zoltan_INT, &
  ZOLTAN_QUAD_QUAD     = 3_Zoltan_INT, &
  ZOLTAN_HEX3D_OCT     = 4_Zoltan_INT

! Error codes for LB library
! These values must agree with the values in zoltan.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(Zoltan_INT) :: &
#else
integer(Zoltan_INT), parameter :: &
#endif
   ZOLTAN_OK     =  0_Zoltan_INT, &
   ZOLTAN_WARN   =  1_Zoltan_INT, &
   ZOLTAN_FATAL  = -1_Zoltan_INT, &
   ZOLTAN_MEMERR = -2_Zoltan_INT

!--------------------------------------------------------------------------
! defined constants for internal use

integer, parameter :: stderr = 6

!--------------------------------------------------------------------------
! interface blocks for the C wrapper functions

interface
!NAS$ ALIEN "F77 zfw_get_address_int"
subroutine Zfw_Get_Address_int(arg,ret_addr)
use zoltan_types
use lb_user_const
use zoltan_user_data
integer(Zoltan_INT) :: arg
integer(Zoltan_INT_PTR), intent(out) :: ret_addr
end subroutine Zfw_Get_Address_int
end interface

interface
!NAS$ ALIEN "F77 zfw_initialize"
function Zfw_Initialize(ver)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Initialize
real(Zoltan_FLOAT), intent(out) :: ver
end function Zfw_Initialize
end interface

interface
!NAS$ ALIEN "F77 zfw_initialize1"
function Zfw_Initialize1(argc,argv,starts,ver)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Initialize1
integer(Zoltan_INT) INTENT_IN argc
integer(Zoltan_INT), dimension(*) INTENT_IN argv, starts
real(Zoltan_FLOAT), intent(out) :: ver
end function Zfw_Initialize1
end interface

interface
!NAS$ ALIEN "F77 zfw_create"
subroutine Zfw_Create(communicator,zz,nbytes)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer INTENT_IN communicator
integer(Zoltan_INT), dimension(*), intent(out) :: zz
integer(Zoltan_INT) INTENT_IN nbytes
end subroutine Zfw_Create
end interface

interface
!NAS$ ALIEN "F77 zfw_destroy"
subroutine Zfw_Destroy(zz,nbytes)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
end subroutine Zfw_Destroy
end interface

interface
!NAS$ ALIEN "F77 zfw_memory_stats"
subroutine Zfw_Memory_Stats()
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
end subroutine Zfw_Memory_Stats
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn0f"
function Zfw_Set_Fn0f(zz,nbytes,fn_type,fn_ptr)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn0f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
end function Zfw_Set_Fn0f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn0s"
function Zfw_Set_Fn0s(zz,nbytes,fn_type,fn_ptr)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn0s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
end function Zfw_Set_Fn0s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn1f"
function Zfw_Set_Fn1f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn1f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
integer(Zoltan_INT), dimension(*) INTENT_IN data
end function Zfw_Set_Fn1f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn1s"
function Zfw_Set_Fn1s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn1s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
integer(Zoltan_INT), dimension(*) INTENT_IN data
end function Zfw_Set_Fn1s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn2f"
function Zfw_Set_Fn2f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn2f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
real(Zoltan_FLOAT), dimension(*) INTENT_IN data
end function Zfw_Set_Fn2f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn2s"
function Zfw_Set_Fn2s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn2s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
real(Zoltan_FLOAT), dimension(*) INTENT_IN data
end function Zfw_Set_Fn2s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn3f"
function Zfw_Set_Fn3f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn3f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
real(Zoltan_DOUBLE), dimension(*) INTENT_IN data
end function Zfw_Set_Fn3f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn3s"
function Zfw_Set_Fn3s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn3s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
real(Zoltan_DOUBLE), dimension(*) INTENT_IN data
end function Zfw_Set_Fn3s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn4f"
function Zfw_Set_Fn4f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn4f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
type(Zoltan_User_Data_1) INTENT_IN data
end function Zfw_Set_Fn4f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn4s"
function Zfw_Set_Fn4s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn4s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(Zoltan_User_Data_1) INTENT_IN data
end function Zfw_Set_Fn4s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn5f"
function Zfw_Set_Fn5f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn5f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
type(Zoltan_User_Data_2) INTENT_IN data
end function Zfw_Set_Fn5f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn5s"
function Zfw_Set_Fn5s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn5s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(Zoltan_User_Data_2) INTENT_IN data
end function Zfw_Set_Fn5s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn6f"
function Zfw_Set_Fn6f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn6f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
type(Zoltan_User_Data_3) INTENT_IN data
end function Zfw_Set_Fn6f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn6s"
function Zfw_Set_Fn6s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn6s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(Zoltan_User_Data_3) INTENT_IN data
end function Zfw_Set_Fn6s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn7f"
function Zfw_Set_Fn7f(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn7f
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
integer(Zoltan_INT), external :: fn_ptr
#endif
type(Zoltan_User_Data_4) INTENT_IN data
end function Zfw_Set_Fn7f
end interface

interface
!NAS$ ALIEN "F77 zfw_set_fn7s"
function Zfw_Set_Fn7s(zz,nbytes,fn_type,fn_ptr,data)
use zoltan_types
use lb_user_const
use zoltan_user_data
#ifdef NASOFTWARE
use nas_system
#endif
implicit none
integer(Zoltan_INT) :: Zfw_Set_Fn7s
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, fn_type
#ifdef NASOFTWARE
type(address), intent(in) :: fn_ptr
#else
external fn_ptr
#endif
type(Zoltan_User_Data_4) INTENT_IN data
end function Zfw_Set_Fn7s
end interface

interface
!NAS$ ALIEN "F77 zfw_set_param"
function Zfw_Set_Param(zz,nbytes,param_name,param_name_len, &
                          new_value,new_value_len)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Param
integer(Zoltan_INT), dimension(*) INTENT_IN zz, param_name, new_value
integer(Zoltan_INT) INTENT_IN nbytes, param_name_len, new_value_len
end function Zfw_Set_Param
end interface

interface
!NAS$ ALIEN "F77 zfw_balance"
function Zfw_LB_Balance(zz,nbytes,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Balance
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT), intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
end function Zfw_LB_Balance
end interface


interface
!NAS$ ALIEN "F77 zfw_eval"
function Zfw_LB_Eval(zz,nbytes,print_stats,nobj,obj_wgt, &
                      ncuts,cut_wgt,nboundary,nadj,is_nobj, &
                      is_obj_wgt,is_ncuts,is_cut_wgt,is_nboundary,is_nadj)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Eval
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes, print_stats
integer(Zoltan_INT), intent(out) :: nobj, ncuts, nboundary, nadj
real(Zoltan_FLOAT), intent(out) :: obj_wgt(*), cut_wgt(*)
integer(Zoltan_INT), intent(in) :: is_nobj, is_ncuts, is_cut_wgt, &
                               is_nboundary, is_nadj, is_obj_wgt
end function Zfw_LB_Eval
end interface

interface
!NAS$ ALIEN "F77 zfw_point_assign"
function Zfw_LB_Point_Assign(zz,nbytes,coords,proc)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Point_Assign
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
real(Zoltan_DOUBLE), dimension(*) INTENT_IN coords
integer(Zoltan_INT), intent(out) :: proc
end function Zfw_LB_Point_Assign
end interface

interface
!NAS$ ALIEN "F77 zfw_box_assign"
function Zfw_LB_Box_Assign(zz,nbytes,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Box_Assign
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
real(Zoltan_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), dimension(*), intent(out) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
end function Zfw_LB_Box_Assign
end interface

interface
!NAS$ ALIEN "F77 zfw_compute_destinations"
function Zfw_Compute_Destinations(zz,nbytes, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Compute_Destinations
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN num_import
integer(Zoltan_INT), intent(out) :: num_export
integer(Zoltan_INT), dimension(*) INTENT_IN import_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: export_global_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: export_local_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_procs
integer(Zoltan_INT), pointer, dimension(:) :: export_procs
end function Zfw_Compute_Destinations
end interface

interface
!NAS$ ALIEN "F77 zfw_help_migrate"
function Zfw_Help_Migrate(zz,nbytes, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Help_Migrate
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN num_import, num_export
integer(Zoltan_INT), dimension(*) INTENT_IN import_global_ids, export_global_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_procs, export_procs
end function Zfw_Help_Migrate
end interface

interface
!NAS$ ALIEN "F77 zfw_register_fort_malloc"
subroutine Zfw_Register_Fort_Malloc(malloc_int,free_int)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
#ifdef NASOFTWARE
type(address), intent(in) :: malloc_int, free_int
#else
external malloc_int,free_int
#endif
end subroutine Zfw_Register_Fort_Malloc
end interface

interface
!NAS$ ALIEN "F77 zfw_get_wgt_dim"
function Zfw_Get_Wgt_Dim(zz,nbytes)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Get_Wgt_Dim
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
end function Zfw_Get_Wgt_Dim
end interface

interface
!NAS$ ALIEN "F77 zfw_get_comm_dim"
function Zfw_Get_Comm_Dim(zz,nbytes)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Get_Comm_Dim
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
end function Zfw_Get_Comm_Dim
end interface

! TEMP child_order
interface
!NAS$ ALIEN "F77 zfw_get_child_order"
subroutine Zfw_Reftree_Get_Child_Order(zz,nbytes,order,ierr)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT), intent(inout), dimension(*) :: order
integer(Zoltan_INT), intent(out) :: ierr
end subroutine Zfw_Reftree_Get_Child_Order
end interface
! end TEMP child_order

!--------------------------------------------------------------------------
! generic names for the Fortran wrapper procedures

interface Zoltan_Initialize
   module procedure Zf90_Initialize
   module procedure Zf90_Initialize1
end interface

interface Zoltan_Create
   module procedure Zf90_Create
end interface

interface Zoltan_Destroy
   module procedure Zf90_Destroy
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

interface Zoltan_LB_Balance
   module procedure Zf90_LB_Balance
end interface

interface Zoltan_LB_Eval
   module procedure Zf90_LB_Eval
end interface

interface Zoltan_LB_Free_Data
   module procedure Zf90_LB_Free_Data
end interface

interface Zoltan_LB_Point_Assign
   module procedure Zf90_LB_Point_Assign
end interface

interface Zoltan_LB_Box_Assign
   module procedure Zf90_LB_Box_Assign
end interface

interface Zoltan_Compute_Destinations
   module procedure Zf90_Compute_Destinations
end interface

interface Zoltan_Help_Migrate
   module procedure Zf90_Help_Migrate
end interface

! TEMP child_order
interface Zoltan_Get_Child_Order
   module procedure Zf90_Reftree_Get_Child_Order
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
#include "set_objsizemulti.if"
#include "set_packobjmulti.if"
#include "set_unpackobjmulti.if"
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

!--------------------------------------------------------------------------
! Fortran wrapper procedures

function Zf90_Initialize(ver)
integer(Zoltan_INT) :: Zf90_Initialize
real(Zoltan_FLOAT), intent(out) :: ver
#ifdef NASOFTWARE
call Zfw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_free_int))
#else
call Zfw_Register_Fort_Malloc(fort_malloc_int,fort_free_int)
#endif
Zf90_Initialize = Zfw_Initialize(ver)
end function Zf90_Initialize

function Zf90_Initialize1(argc,argv,ver)
integer(Zoltan_INT) :: Zf90_Initialize1
integer(Zoltan_INT) INTENT_IN argc
character(len=*), dimension(*) INTENT_IN argv
real(Zoltan_FLOAT), intent(out) :: ver
integer(Zoltan_INT), allocatable, dimension(:) :: int_argv,starts
integer(Zoltan_INT) :: i, j, leng
#ifdef NASOFTWARE
call Zfw_Register_Fort_Malloc(loc(fort_malloc_int),loc(fort_free_int))
#else
call Zfw_Register_Fort_Malloc(fort_malloc_int,fort_free_int)
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
Zf90_Initialize1 = Zfw_Initialize1(argc,int_argv,starts,ver)
deallocate(starts,int_argv)
end function Zf90_Initialize1

function Zf90_Create(communicator)
type(Zoltan_Struct), pointer :: Zf90_Create
integer INTENT_IN communicator
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

subroutine Zf90_Memory_Stats()
call Zfw_Memory_Stats()
end subroutine Zf90_Memory_Stats

function Zf90_Set_Fn0f(zz,fn_type,fn_ptr)
integer(Zoltan_INT) :: Zf90_Set_Fn0f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn0f = Zfw_Set_Fn0f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
Zf90_Set_Fn0f = Zfw_Set_Fn0f(zz_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function Zf90_Set_Fn0f

function Zf90_Set_Fn0s(zz,fn_type,fn_ptr)
integer(Zoltan_INT) :: Zf90_Set_Fn0s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn0s = Zfw_Set_Fn0s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr))
#else
Zf90_Set_Fn0s = Zfw_Set_Fn0s(zz_addr,nbytes,fn_type%choice,fn_ptr)
#endif
end function Zf90_Set_Fn0s

function Zf90_Set_Fn1f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn1f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
integer(Zoltan_INT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn1f = Zfw_Set_Fn1f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn1f = Zfw_Set_Fn1f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn1f

function Zf90_Set_Fn1s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn1s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
integer(Zoltan_INT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn1s = Zfw_Set_Fn1s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn1s = Zfw_Set_Fn1s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn1s

function Zf90_Set_Fn2f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn2f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_FLOAT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn2f = Zfw_Set_Fn2f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn2f = Zfw_Set_Fn2f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn2f

function Zf90_Set_Fn2s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn2s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(Zoltan_FLOAT) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn2s = Zfw_Set_Fn2s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn2s = Zfw_Set_Fn2s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn2s

function Zf90_Set_Fn3f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn3f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
real(Zoltan_DOUBLE) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn3f = Zfw_Set_Fn3f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn3f = Zfw_Set_Fn3f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn3f

function Zf90_Set_Fn3s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn3s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
real(Zoltan_DOUBLE) INTENT_IN data(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn3s = Zfw_Set_Fn3s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn3s = Zfw_Set_Fn3s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn3s

function Zf90_Set_Fn4f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn4f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_1) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn4f = Zfw_Set_Fn4f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn4f = Zfw_Set_Fn4f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn4f

function Zf90_Set_Fn4s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn4s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(Zoltan_User_Data_1) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn4s = Zfw_Set_Fn4s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn4s = Zfw_Set_Fn4s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn4s

function Zf90_Set_Fn5f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn5f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_2) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn5f = Zfw_Set_Fn5f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn5f = Zfw_Set_Fn5f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn5f

function Zf90_Set_Fn5s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn5s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(Zoltan_User_Data_2) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn5s = Zfw_Set_Fn5s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn5s = Zfw_Set_Fn5s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn5s

function Zf90_Set_Fn6f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn6f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_3) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn6f = Zfw_Set_Fn6f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn6f = Zfw_Set_Fn6f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn6f

function Zf90_Set_Fn6s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn6s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(Zoltan_User_Data_3) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn6s = Zfw_Set_Fn6s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn6s = Zfw_Set_Fn6s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn6s

function Zf90_Set_Fn7f(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn7f
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPEF) INTENT_IN fn_type
integer(Zoltan_INT), external :: fn_ptr
type(Zoltan_User_Data_4) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn7f = Zfw_Set_Fn7f(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn7f = Zfw_Set_Fn7f(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn7f

function Zf90_Set_Fn7s(zz,fn_type,fn_ptr,data)
integer(Zoltan_INT) :: Zf90_Set_Fn7s
type(Zoltan_Struct) INTENT_IN zz
type(ZOLTAN_FN_TYPES) INTENT_IN fn_type
external fn_ptr
type(Zoltan_User_Data_4) INTENT_IN data
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
#ifdef NASOFTWARE
Zf90_Set_Fn7s = Zfw_Set_Fn7s(zz_addr,nbytes,fn_type%choice,loc(fn_ptr),data)
#else
Zf90_Set_Fn7s = Zfw_Set_Fn7s(zz_addr,nbytes,fn_type%choice,fn_ptr,data)
#endif
end function Zf90_Set_Fn7s

function Zf90_Set_Param(zz,param_name,new_value)
integer(Zoltan_INT) :: Zf90_Set_Param
type(Zoltan_Struct) INTENT_IN zz
character(len=*) INTENT_IN param_name, new_value
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

function Zf90_LB_Balance(zz,changes,num_gid_entries,num_lid_entries, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_LB_Balance
type(Zoltan_Struct) INTENT_IN zz
logical, intent(out) :: changes
integer(Zoltan_INT), intent(out) :: num_gid_entries, num_lid_entries
integer(Zoltan_INT), intent(out) :: num_import, num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i, int_changes
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Balance = Zfw_LB_Balance(zz_addr,nbytes,int_changes, &
                             num_gid_entries, num_lid_entries, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
changes = .not.(int_changes==0)
end function Zf90_LB_Balance

function Zf90_LB_Eval(zz,print_stats,nobj,obj_wgt, &
                    ncuts,cut_wgt,nboundary,nadj)
integer(Zoltan_INT) :: Zf90_LB_Eval
type(Zoltan_Struct) INTENT_IN zz
logical INTENT_IN print_stats
integer(Zoltan_INT), intent(out), optional :: nobj, ncuts, nboundary, nadj
real(Zoltan_FLOAT), intent(out), optional :: obj_wgt(*), cut_wgt(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i, int_print_stats, dim, edim
integer(Zoltan_INT) :: loc_nobj, loc_ncuts, loc_nboundary, loc_nadj
real(Zoltan_FLOAT), allocatable :: loc_obj_wgt(:), loc_cut_wgt(:)
integer(Zoltan_INT) :: is_nobj, is_ncuts, is_cut_wgt, is_nboundary, is_nadj, is_obj_wgt
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
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
   dim = Zfw_Get_Wgt_Dim(zz_addr,nbytes)
   allocate(loc_obj_wgt(dim))
else
   is_obj_wgt = 0
   allocate(loc_obj_wgt(1))
endif
if (present(cut_wgt)) then
   is_cut_wgt = 1
   edim = Zfw_Get_Comm_Dim(zz_addr,nbytes)
   allocate(loc_cut_wgt(edim))
else
   is_cut_wgt = 0
   allocate(loc_cut_wgt(1))
endif
Zf90_LB_Eval = Zfw_LB_Eval(zz_addr,nbytes,int_print_stats,loc_nobj,loc_obj_wgt, &
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
end function Zf90_LB_Eval

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

function Zf90_LB_Point_Assign(zz,coords,proc)
integer(Zoltan_INT) :: Zf90_LB_Point_Assign
type(Zoltan_Struct) INTENT_IN zz
real(Zoltan_DOUBLE), dimension(*) INTENT_IN coords
integer(Zoltan_INT), intent(out) :: proc
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_LB_Point_Assign = Zfw_LB_Point_Assign(zz_addr,nbytes,coords,proc)
end function Zf90_LB_Point_Assign

function Zf90_LB_Box_Assign(zz,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs)
integer(Zoltan_INT) :: Zf90_LB_Box_Assign
type(Zoltan_Struct) INTENT_IN zz
real(Zoltan_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
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

function Zf90_Compute_Destinations(zz, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_Compute_Destinations
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN num_import
integer(Zoltan_INT), intent(out) :: num_export
integer(Zoltan_INT), pointer, dimension(:) :: import_global_ids, export_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_local_ids, export_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: import_procs, export_procs
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
if (.not.associated(import_global_ids) .or. .not.associated(import_local_ids) &
    .or. .not.associated(import_procs)) then
   write(stderr,*) "Error from Zoltan_Compute_Destinations: import pointers are not associated"
   Zf90_Compute_Destinations = ZOLTAN_WARN
   return
endif
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Compute_Destinations = Zfw_Compute_Destinations(zz_addr,nbytes, &
                             num_import,import_global_ids,import_local_ids, &
                             import_procs,num_export,export_global_ids, &
                             export_local_ids,export_procs)
end function Zf90_Compute_Destinations


function Zf90_Help_Migrate(zz, &
                       num_import,import_global_ids, &
                       import_local_ids,import_procs,num_export, &
                       export_global_ids,export_local_ids,export_procs)
integer(Zoltan_INT) :: Zf90_Help_Migrate
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN num_import, num_export
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

! TEMP child_order
subroutine Zf90_Reftree_Get_Child_Order(zz,order,ierr)
type(Zoltan_Struct), pointer :: zz
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
#include "set_objsizemulti.fn"
#include "set_packobjmulti.fn"
#include "set_unpackobjmulti.fn"
#include "set_numcoarseobj.fn"
#include "set_coarseobjlist.fn"
#include "set_firstcoarseobj.fn"
#include "set_nextcoarseobj.fn"
#include "set_numchild.fn"
#include "set_childlist.fn"
#include "set_childweight.fn"

!-------------------------------------------------------------------------
! Include LB_* interface for backward compatibility.

#include "lbfn.fpp"
#include "set_numgeom.fn.lbfn"
#include "set_geom.fn.lbfn"
#include "set_numedges.fn.lbfn"
#include "set_edgelist.fn.lbfn"
#include "set_numobj.fn.lbfn"
#include "set_objlist.fn.lbfn"
#include "set_firstobj.fn.lbfn"
#include "set_nextobj.fn.lbfn"
#include "set_numborderobj.fn.lbfn"
#include "set_borderobjlist.fn.lbfn"
#include "set_firstborderobj.fn.lbfn"
#include "set_nextborderobj.fn.lbfn"
#include "set_premigrate.fn.lbfn"
#include "set_midmigrate.fn.lbfn"
#include "set_postmigrate.fn.lbfn"
#include "set_objsize.fn.lbfn"
#include "set_packobj.fn.lbfn"
#include "set_unpackobj.fn.lbfn"
#include "set_numcoarseobj.fn.lbfn"
#include "set_coarseobjlist.fn.lbfn"
#include "set_firstcoarseobj.fn.lbfn"
#include "set_nextcoarseobj.fn.lbfn"
#include "set_numchild.fn.lbfn"
#include "set_childlist.fn.lbfn"
#include "set_childweight.fn.lbfn"

end module zoltan
