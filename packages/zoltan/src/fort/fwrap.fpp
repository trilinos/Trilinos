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
   ZOLTAN_PARTITION_FN_TYPE, &
   ZOLTAN_PARTITION_MULTI_FN_TYPE, &
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
   Zoltan_Generate_Files, &
   Zoltan_RCB_Box

! Registration functions with strict type checking.
public :: &
   Zoltan_Set_Num_Obj_Fn, Zoltan_Set_Obj_List_Fn, &
   Zoltan_Set_First_Obj_Fn, Zoltan_Set_Next_Obj_Fn, &
   Zoltan_Set_Num_Border_Obj_Fn, Zoltan_Set_Border_Obj_List_Fn, &
   Zoltan_Set_First_Border_Obj_Fn, Zoltan_Set_Next_Border_Obj_Fn, &
   Zoltan_Set_Num_Geom_Fn, Zoltan_Set_Geom_Multi_Fn, Zoltan_Set_Geom_Fn, &
   Zoltan_Set_Partition_Fn, Zoltan_Set_Partition_Multi_Fn, &
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
   ZOLTAN_PARTITION_FN_TYPE        = ZOLTAN_FN_TYPEF(34_Zoltan_INT)

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPES) :: &
#else
type(ZOLTAN_FN_TYPES), parameter :: &
#endif
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
   ZOLTAN_PARTITION_MULTI_FN_TYPE  = ZOLTAN_FN_TYPES(35_Zoltan_INT)

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
!NAS$ ALIEN "F77 zfw_set_param_vec"
function Zfw_Set_Param_Vec(zz,nbytes,param_name,param_name_len, &
                          new_value,new_value_len,index)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Set_Param_Vec
integer(Zoltan_INT), dimension(*) INTENT_IN zz, param_name, new_value
integer(Zoltan_INT) INTENT_IN nbytes, param_name_len, new_value_len
integer(Zoltan_INT) INTENT_IN index
end function Zfw_Set_Param_Vec
end interface

interface
!NAS$ ALIEN "F77 zfw_partition"
function Zfw_LB_Partition(zz,nbytes,changes,num_gid_entries,num_lid_entries, &
                num_import,import_global_ids, &
                import_local_ids,import_procs,import_to_part,num_export, &
                export_global_ids,export_local_ids,export_procs,export_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Partition
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
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
!NAS$ ALIEN "F77 zfw_set_part_sizes"
function Zfw_LB_Set_Part_Sizes(zz,nbytes,global_part,len,partids,&
                               wgtidx,partsizes)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Set_Part_Sizes
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes,global_part,len,partids(*),wgtidx(*)
real(Zoltan_FLOAT) INTENT_IN partsizes(*)
end function Zfw_LB_Set_Part_Sizes
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
!NAS$ ALIEN "F77 zfw_point_pp_assign"
function Zfw_LB_Point_PP_Assign(zz,nbytes,coords,proc,part)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Point_PP_Assign
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
real(Zoltan_DOUBLE), dimension(*) INTENT_IN coords
integer(Zoltan_INT), intent(out) :: proc
integer(Zoltan_INT), intent(out) :: part
end function Zfw_LB_Point_PP_Assign
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
!NAS$ ALIEN "F77 zfw_box_pp_assign"
function Zfw_LB_Box_PP_Assign(zz,nbytes,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs,parts,numparts)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_LB_Box_PP_Assign
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
real(Zoltan_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
integer(Zoltan_INT), dimension(*), intent(out) :: procs
integer(Zoltan_INT), intent(out) :: numprocs
integer(Zoltan_INT), dimension(*), intent(out) :: parts
integer(Zoltan_INT), intent(out) :: numparts
end function Zfw_LB_Box_PP_Assign
end interface

interface
!NAS$ ALIEN "F77 zfw_invert_lists"
function Zfw_Invert_Lists(zz,nbytes, &
                       num_input,input_global_ids,input_local_ids, &
                       input_procs,input_to_part, &
                       num_output,output_global_ids,output_local_ids, &
                       output_procs,output_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Invert_Lists
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), dimension(*) INTENT_IN input_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_global_ids
integer(Zoltan_INT), dimension(*) INTENT_IN input_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_local_ids
integer(Zoltan_INT), dimension(*) INTENT_IN input_procs
integer(Zoltan_INT), pointer, dimension(:) :: output_procs
integer(Zoltan_INT), dimension(*) INTENT_IN input_to_part
integer(Zoltan_INT), pointer, dimension(:) :: output_to_part
end function Zfw_Invert_Lists
end interface

interface
!NAS$ ALIEN "F77 zfw_compute_destinations"
function Zfw_Compute_Destinations(zz,nbytes, &
                       num_input,input_global_ids, &
                       input_local_ids,input_procs,num_output, &
                       output_global_ids,output_local_ids,output_procs)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Compute_Destinations
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN num_input
integer(Zoltan_INT), intent(out) :: num_output
integer(Zoltan_INT), dimension(*) INTENT_IN input_global_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_global_ids
integer(Zoltan_INT), dimension(*) INTENT_IN input_local_ids
integer(Zoltan_INT), pointer, dimension(:) :: output_local_ids
integer(Zoltan_INT), dimension(*) INTENT_IN input_procs
integer(Zoltan_INT), pointer, dimension(:) :: output_procs
end function Zfw_Compute_Destinations
end interface

interface
!NAS$ ALIEN "F77 zfw_migrate"
function Zfw_Migrate(zz,nbytes, &
                     num_import,import_global_ids,import_local_ids, &
                     import_procs,import_to_part, &
                     num_export,export_global_ids,export_local_ids, &
                     export_procs,export_to_part)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Migrate
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN num_import, num_export
integer(Zoltan_INT), dimension(*) INTENT_IN import_global_ids, export_global_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_local_ids, export_local_ids
integer(Zoltan_INT), dimension(*) INTENT_IN import_procs, export_procs
integer(Zoltan_INT), dimension(*) INTENT_IN import_to_part, export_to_part
end function Zfw_Migrate
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
!NAS$ ALIEN "F77 zfw_order"
function Zfw_Order(zz,nbytes,num_gid_entries,num_lid_entries,num_obj, &
                   gids,lids,rank,iperm)
use zoltan_types
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Order
INTEGER(Zoltan_INT), dimension(*), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(IN) :: nbytes
INTEGER(Zoltan_INT), INTENT(OUT) :: num_gid_entries, num_lid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*), lids(*)
INTEGER(Zoltan_INT) :: rank(*), iperm(*)
end function Zfw_Order
end interface

interface
!NAS$ ALIEN "F77 zfw_generate_files"
function Zfw_Generate_Files(zz,nbytes,filename,filename_len, &
                          base_index, gen_geom, gen_graph, gen_hg)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_Generate_Files
integer(Zoltan_INT), dimension(*) INTENT_IN zz, filename
integer(Zoltan_INT) INTENT_IN nbytes, filename_len, base_index
integer(Zoltan_INT) INTENT_IN gen_geom, gen_graph, gen_hg
end function Zfw_Generate_Files
end interface

interface
!NAS$ ALIEN "F77 zfw_rcb_box"
function Zfw_RCB_Box(zz,nbytes,part,ndim,xmin,ymin,zmin,xmax,ymax,zmax)
use zoltan_types
use lb_user_const
use zoltan_user_data
implicit none
integer(Zoltan_INT) :: Zfw_RCB_Box
integer(Zoltan_INT), dimension(*) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN nbytes
integer(Zoltan_INT) INTENT_IN part
integer(Zoltan_INT), intent(out) :: ndim
real(Zoltan_DOUBLE), intent(out) :: xmin,ymin,zmin,xmax,ymax,zmax
end function Zfw_RCB_Box
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

interface Zoltan_Generate_Files
   module procedure Zf90_Generate_Files
end interface

! TEMP child_order
interface Zoltan_Get_Child_Order
   module procedure Zf90_Reftree_Get_Child_Order
end interface

#include "set_numgeom.if"
#include "set_geommulti.if"
#include "set_geom.if"
#include "set_partition.if"
#include "set_partitionmulti.if"
#include "set_numedges.if"
#include "set_numedgesmulti.if"
#include "set_edgelist.if"
#include "set_edgelistmulti.if"
#include "set_numobj.if"
#include "set_objlist.if"
#include "set_firstobj.if"
#include "set_nextobj.if"
#include "set_numborderobj.if"
#include "set_borderobjlist.if"
#include "set_firstborderobj.if"
#include "set_nextborderobj.if"
#include "set_premigratepp.if"
#include "set_midmigratepp.if"
#include "set_postmigratepp.if"
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

!--------------------------------------------------------------------------
! Fortran wrapper procedures
!--------------------------------------------------------------------------

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
subroutine Zf90_Memory_Stats()
call Zfw_Memory_Stats()
end subroutine Zf90_Memory_Stats

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Set_Param_Vec(zz,param_name,new_value,index)
integer(Zoltan_INT) :: Zf90_Set_Param_Vec
type(Zoltan_Struct) INTENT_IN zz
character(len=*) INTENT_IN param_name, new_value
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
type(Zoltan_Struct) INTENT_IN zz
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
type(Zoltan_Struct) INTENT_IN zz
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
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN global_part,len,partids(*),wgtidx(*)
real(Zoltan_FLOAT) INTENT_IN partsizes(*)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Point_PP_Assign(zz,coords,proc,part)
integer(Zoltan_INT) :: Zf90_LB_Point_PP_Assign
type(Zoltan_Struct) INTENT_IN zz
real(Zoltan_DOUBLE), dimension(*) INTENT_IN coords
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_LB_Box_PP_Assign(zz,xmin,ymin,zmin,xmax,ymax,zmax,procs,numprocs,parts,numparts)
integer(Zoltan_INT) :: Zf90_LB_Box_PP_Assign
type(Zoltan_Struct) INTENT_IN zz
real(Zoltan_DOUBLE) INTENT_IN xmin,ymin,zmin,xmax,ymax,zmax
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
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN num_input
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
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN num_input
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
type(Zoltan_Struct) INTENT_IN zz
integer(Zoltan_INT) INTENT_IN num_import, num_export
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Order(zz,num_gid_entries,num_lid_entries,num_obj,gids,lids,rank,iperm)
integer(Zoltan_INT) :: Zf90_Order
TYPE(Zoltan_Struct), INTENT(IN) :: zz 
INTEGER(Zoltan_INT), INTENT(OUT) :: num_gid_entries, num_lid_entries
INTEGER(Zoltan_INT), INTENT(IN) :: num_obj
INTEGER(Zoltan_INT) :: gids(*), lids(*)
INTEGER(Zoltan_INT) :: rank(*), iperm(*)
integer(Zoltan_INT), dimension(Zoltan_PTR_LENGTH) :: zz_addr
integer(Zoltan_INT) :: nbytes, i
nbytes = Zoltan_PTR_LENGTH
do i=1,nbytes
   zz_addr(i) = ichar(zz%addr%addr(i:i))
end do
Zf90_Order = Zfw_Order(zz_addr,nbytes,num_gid_entries,num_lid_entries,num_obj,&
                       gids,lids,rank,iperm)
end function Zf90_Order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Zf90_Generate_Files(zz,filename,base_index,gen_geom,gen_graph,gen_hg)
integer(Zoltan_INT) :: Zf90_Generate_Files
type(Zoltan_Struct) INTENT_IN zz
character(len=*) INTENT_IN filename
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
type(Zoltan_Struct) INTENT_IN zz
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
! TEMP child_order
subroutine Zf90_Reftree_Get_Child_Order(zz,order,ierr)
type(Zoltan_Struct) INTENT_IN zz
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#include "set_numgeom.fn"
#include "set_geommulti.fn"
#include "set_geom.fn"
#include "set_partition.fn"
#include "set_partitionmulti.fn"
#include "set_numedges.fn"
#include "set_numedgesmulti.fn"
#include "set_edgelist.fn"
#include "set_edgelistmulti.fn"
#include "set_numobj.fn"
#include "set_objlist.fn"
#include "set_firstobj.fn"
#include "set_nextobj.fn"
#include "set_numborderobj.fn"
#include "set_borderobjlist.fn"
#include "set_firstborderobj.fn"
#include "set_nextborderobj.fn"
#include "set_premigratepp.fn"
#include "set_midmigratepp.fn"
#include "set_postmigratepp.fn"
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
!-------------------------------------------------------------------------

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
