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

public :: &
   LB_INT, &
   LB_FLOAT, &
   LB_DOUBLE, &
   LB_User_Data_1, &
   LB_User_Data_2, &
   LB_User_Data_3, &
   LB_User_Data_4

public :: &
   LB_Struct

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
   LB_MID_MIGRATE_FN_TYPE, &
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
   LB_Memory_Stats, &
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

! Registration functions with strict type checking.
public :: &
   LB_Set_Num_Obj_Fn, LB_Set_Obj_List_Fn, &
   LB_Set_First_Obj_Fn, LB_Set_Next_Obj_Fn, &
   LB_Set_Num_Border_Obj_Fn, LB_Set_Border_Obj_List_Fn, &
   LB_Set_First_Border_Obj_Fn, LB_Set_Next_Border_Obj_Fn, &
   LB_Set_Num_Geom_Fn, LB_Set_Geom_Fn, &
   LB_Set_Num_Edges_Fn, LB_Set_Edge_List_Fn, &
   LB_Set_Num_Coarse_Obj_Fn, LB_Set_Coarse_Obj_List_Fn, &
   LB_Set_First_Coarse_Obj_Fn, LB_Set_Next_Coarse_Obj_Fn, &
   LB_Set_Num_Child_Fn, LB_Set_Child_List_Fn, LB_Set_Child_Weight_Fn, &
   LB_Set_Pre_Migrate_Fn, LB_Set_Mid_Migrate_Fn, LB_Set_Post_Migrate_Fn, &
   LB_Set_Obj_Size_Fn, LB_Set_Pack_Obj_Fn, LB_Set_Unpack_Obj_Fn 

public :: &
   LB_Get_Child_Order ! TEMP child_order

!--------------------------------------------------------------------------
! user defined types corresponding to the C structs

type LB_Struct
   private
   sequence
   type(Zoltan_PTR) :: addr
#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
   integer :: dummy
#endif
end type LB_Struct

!--------------------------------------------------------------------------
! defined constants corresponding to Zoltan enumerated types

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPEF) :: &
#else
type(ZOLTAN_FN_TYPEF), parameter :: &
#endif
! KDD:  The following is not the most elegant way to write these assignments.
! KDD:  I'd prefer, e.g.,  LB_NUM_EDGES_FN_TYPE = ZOLTAN_NUM_EDGES_FN_TYPE,
! KDD:  Such assignment works for integer, but for user-defined types
! KDD:  the pgf90 compiler for tflops either reports an invalid assignment
! KDD:  or core dumps when we use the above assignments.
   LB_NUM_EDGES_FN_TYPE        = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_EDGES_FN_TYPE%choice), &
   LB_NUM_GEOM_FN_TYPE         = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_GEOM_FN_TYPE%choice), &
   LB_NUM_OBJ_FN_TYPE          = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_OBJ_FN_TYPE%choice), &
   LB_FIRST_OBJ_FN_TYPE        = &
      ZOLTAN_FN_TYPEF(ZOLTAN_FIRST_OBJ_FN_TYPE%choice), &
   LB_NEXT_OBJ_FN_TYPE         = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NEXT_OBJ_FN_TYPE%choice), &
   LB_NUM_BORDER_OBJ_FN_TYPE   = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_BORDER_OBJ_FN_TYPE%choice), &
   LB_FIRST_BORDER_OBJ_FN_TYPE = &
      ZOLTAN_FN_TYPEF(ZOLTAN_FIRST_BORDER_OBJ_FN_TYPE%choice), &
   LB_NEXT_BORDER_OBJ_FN_TYPE  = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NEXT_BORDER_OBJ_FN_TYPE%choice), &
   LB_OBJ_SIZE_FN_TYPE         = &
      ZOLTAN_FN_TYPEF(ZOLTAN_OBJ_SIZE_FN_TYPE%choice), &
   LB_NUM_COARSE_OBJ_FN_TYPE   = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_COARSE_OBJ_FN_TYPE%choice), &
   LB_FIRST_COARSE_OBJ_FN_TYPE = &
      ZOLTAN_FN_TYPEF(ZOLTAN_FIRST_COARSE_OBJ_FN_TYPE%choice), &
   LB_NEXT_COARSE_OBJ_FN_TYPE  = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NEXT_COARSE_OBJ_FN_TYPE%choice), &
   LB_NUM_CHILD_FN_TYPE        = &
      ZOLTAN_FN_TYPEF(ZOLTAN_NUM_CHILD_FN_TYPE%choice)

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
type(ZOLTAN_FN_TYPES) :: &
#else
type(ZOLTAN_FN_TYPES), parameter :: &
#endif
! KDD:  The following is not the most elegant way to write these assignments.
! KDD:  I'd prefer, e.g.,  LB_EDGE_LIST_FN_TYPE = ZOLTAN_EDGE_LIST_FN_TYPE,
! KDD:  Such assignment works for integer, but for user-defined types
! KDD:  the pgf90 compiler for tflops either reports an invalid assignment
! KDD:  or core dumps when we use the above assignments.
   LB_EDGE_LIST_FN_TYPE        = &
      ZOLTAN_FN_TYPES(ZOLTAN_EDGE_LIST_FN_TYPE%choice), &
   LB_GEOM_FN_TYPE             = &
      ZOLTAN_FN_TYPES(ZOLTAN_GEOM_FN_TYPE%choice), &
   LB_OBJ_LIST_FN_TYPE         = &
      ZOLTAN_FN_TYPES(ZOLTAN_OBJ_LIST_FN_TYPE%choice), &
   LB_BORDER_OBJ_LIST_FN_TYPE  = &
      ZOLTAN_FN_TYPES(ZOLTAN_BORDER_OBJ_LIST_FN_TYPE%choice), &
   LB_PRE_MIGRATE_FN_TYPE      = &
      ZOLTAN_FN_TYPES(ZOLTAN_PRE_MIGRATE_FN_TYPE%choice), &
   LB_MID_MIGRATE_FN_TYPE      = &
      ZOLTAN_FN_TYPES(ZOLTAN_MID_MIGRATE_FN_TYPE%choice), &
   LB_POST_MIGRATE_FN_TYPE     = &
      ZOLTAN_FN_TYPES(ZOLTAN_POST_MIGRATE_FN_TYPE%choice), &
   LB_PACK_OBJ_FN_TYPE         = &
      ZOLTAN_FN_TYPES(ZOLTAN_PACK_OBJ_FN_TYPE%choice), &
   LB_UNPACK_OBJ_FN_TYPE       = &
      ZOLTAN_FN_TYPES(ZOLTAN_UNPACK_OBJ_FN_TYPE%choice), &
   LB_COARSE_OBJ_LIST_FN_TYPE  = &
      ZOLTAN_FN_TYPES(ZOLTAN_COARSE_OBJ_LIST_FN_TYPE%choice), &
   LB_CHILD_LIST_FN_TYPE       = &
      ZOLTAN_FN_TYPES(ZOLTAN_CHILD_LIST_FN_TYPE%choice), &
   LB_CHILD_WEIGHT_FN_TYPE     = &
      ZOLTAN_FN_TYPES(ZOLTAN_CHILD_WEIGHT_FN_TYPE%choice)

! Type of refinement used when building a refinement tree
! These values must agree with the values in lb/lbi_const.h

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(Zoltan_INT) :: &
#else
integer(Zoltan_INT), parameter :: &
#endif
  LB_OTHER_REF     = ZOLTAN_OTHER_REF, &
  LB_IN_ORDER      = ZOLTAN_IN_ORDER, &
  LB_TRI_BISECT    = ZOLTAN_TRI_BISECT, &
  LB_QUAD_QUAD     = ZOLTAN_QUAD_QUAD, &
  LB_HEX3D_OCT     = ZOLTAN_HEX3D_OCT

#ifdef SUNSOFT
! bug in SunSoft compiler prevents using parameter
integer(Zoltan_INT) :: &
#else
integer(Zoltan_INT), parameter :: &
#endif
   LB_OK     = ZOLTAN_OK, &
   LB_WARN   = ZOLTAN_WARN, &
   LB_FATAL  = ZOLTAN_FATAL, &
   LB_MEMERR = ZOLTAN_MEMERR

!--------------------------------------------------------------------------
! generic names for the Fortran wrapper procedures

interface LB_Initialize
   module procedure Zf90_Initialize
   module procedure Zf90_Initialize1
end interface

interface LB_Create
   module procedure Zf90_Create
end interface

interface LB_Destroy
   module procedure Zf90_Destroy
end interface

interface LB_Memory_Stats
   module procedure Zf90_Memory_Stats
end interface

interface LB_Set_Fn
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

interface LB_Set_Method
   module procedure Zf90_LB_Set_Method
end interface

interface LB_Set_Param
   module procedure Zf90_Set_Param
end interface

interface LB_Balance
   module procedure Zf90_LB_Balance
end interface

interface LB_Eval
   module procedure Zf90_LB_Eval
end interface

interface LB_Free_Data
   module procedure Zf90_LB_Free_Data
end interface

interface LB_Point_Assign
   module procedure Zf90_LB_Point_Assign
end interface

interface LB_Box_Assign
   module procedure Zf90_LB_Box_Assign
end interface

interface LB_Compute_Destinations
   module procedure Zf90_Compute_Destinations
end interface

interface LB_Help_Migrate
   module procedure Zf90_Help_Migrate
end interface

! TEMP child_order
interface LB_Get_Child_Order
   module procedure Zf90_LB_Get_Child_Order
end interface

