module lb_user_const
use zoltan_types
implicit none

!  This file contains the user-defined data types and comparison functions
!  for global and local IDs used by the application and the
!  load-balancing library.  Application developers should modify these
!  data types to match those of identifiers used in their applications.
!
!  In this example, LB_LID is INTEGER(KIND=LB_INT), and
!  LB_GID is a structure containing a pair of integers, with the second
!  integer acting as a high order word, and knowledge of how the compiler
!  implements structures is used to perform operations in C.
!
!  LB_GID are the unique global identifiers for objects in the application.  
!  The LB_GID are used as identifiers within the load-balancing routines
!  as well.  Thus, functions defining methods to compare global identifiers 
!  must be provided.
!
!  LB_LID are local identifiers that are not used by the load-balancing
!  routines.  They are stored with objects in the load-balancing routines,
!  however, and are passed to the application query functions.  An 
!  application can provide any values it wants for local identifiers, and
!  can use them to make access of object information in the query routines
!  more efficient.

!  LB_GID and LB_LID data type definitions.
!  For this example, local IDs (LB_LID) are integers, so this is just a
!  dummy definition to satisfy the linker, with LIDs actually being
!  handled in lb_user_const.h.  LB_GID is defined as a structure and
!  used in the Fortran code, but the copy and comparison operators are
!  defined in lb_user_const.h using knowledge of how the compiler
!  implements structures (this is highly compiler specific).

type LB_GID
   integer(LB_INT) :: int1, int2
end type LB_GID

type LB_LID
   integer(LB_INT) :: id
end type LB_LID

! User defined data types for passing data to the query functions.  These can
! be used any way you want, but one suggestion is to use them as "wrapper"
! types for your own user defined types, e.g.
! type LB_User_Data_1 
!    type(your_type), pointer :: ptr
! end type
! Exactly four data types must be defined, but you don't have to use them.

type LB_User_Data_1
   integer :: dummy
end type LB_User_Data_1

type LB_User_Data_2
   integer :: dummy
end type LB_User_Data_2

type LB_User_Data_3
   integer :: dummy
end type LB_User_Data_3

type LB_User_Data_4
   integer :: dummy
end type LB_User_Data_4

!contains

!  Subroutines to copy LB_GIDs and LB_LIDs.
!  These subroutines are used by the load-balancing routines to copy LB_GID and
!  LB_LID values to new storage locations.
!  In this example, these operations are performed by C macros in
!  lb_user_const_int.h, so the routines are omitted.

!  Functions to compare LB_GIDs.
!  Functions must be provided to test whether two LB_GIDs are equal (EQ),
!  not equal (NE), less than (LT), less than or equal (LE), 
!  greater than (GT), and greater than or equal (GE).
!  The function must return the value 1 if the comparison yields .true. and
!  0 if the comparison yeilds .false.
!  Comparison functions are not needed for LB_LIDs as LB_LIDs are not used
!  within the load-balancing routines.
!  In this example, these operations are performed by C macros in
!  lb_user_const_int.h, so the routines are omitted.

end module lb_user_const
