!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zoltan Dynamic Load-Balancing Library for Parallel Applications            !
! Copyright (c) 2000, Sandia National Laboratories.                          !
! Zoltan is distributed under the GNU Lesser General Public License 2.1.     !
! For more info, see the README file in the top-level Zoltan directory.      ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CVS File Information :
!     $RCSfile$
!     $Author$
!     $Date$
!     $Revision$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!  integer acting as a high order word.
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
!  For this example, local IDs (LB_LID) are integers, so these are just
!  dummy definitions to satisfy the linker, with LIDs actually being
!  handled in lb_user_const.h

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

contains

!  Subroutines to copy LB_GIDs and LB_LIDs.
!  These subroutines are used by the load-balancing routines to copy LB_GID and
!  LB_LID values to new storage locations.
!  In this example, these operations are performed by C macros in
!  lb_user_const_int.h, so these are just dummy routines that will never
!  be called.

subroutine LB_SET_GID(a,b)
type(LB_GID), intent(out) :: a
type(LB_GID), intent(in) :: b
a = b
end subroutine LB_SET_GID

subroutine LB_SET_LID(a,b)
integer(LB_INT), intent(out) :: a
integer(LB_INT), intent(in) :: b
a = b
end subroutine LB_SET_LID

!  Functions to compare LB_GIDs.
!  Functions must be provided to test whether two LB_GIDs are equal (EQ),
!  not equal (NE), less than (LT), less than or equal (LE), 
!  greater than (GT), and greater than or equal (GE).
!  The function must return the value 1 if the comparison yields .true. and
!  0 if the comparison yeilds .false.
!  Comparison functions are not needed for LB_LIDs as LB_LIDs are not used
!  within the load-balancing routines.
!  In this example, these operations are performed by C macros in
!  lb_user_const_int.h, so these are just dummy routines that will never
!  be called.

function LB_EQ_GID(a,b)
integer(LB_INT) :: LB_EQ_GID
type(LB_GID), intent(in) :: a,b
if (a%int1 == b%int1 .and. a%int2 == b%int2) then
   LB_EQ_GID = 1_LB_INT
else
   LB_EQ_GID = 0_LB_INT
endif
end function LB_EQ_GID

function LB_NE_GID(a,b)
integer(LB_INT) :: LB_NE_GID
type(LB_GID), intent(in) :: a,b
if (a%int1 /= b%int1 .or. a%int2 /= b%int2) then
   LB_NE_GID = 1_LB_INT
else
   LB_NE_GID = 0_LB_INT
endif
end function LB_NE_GID

function LB_LT_GID(a,b)
integer(LB_INT) :: LB_LT_GID
type(LB_GID), intent(in) :: a,b
if (a%int2 == b%int2) then
   if (a%int1 < b%int1) then
      LB_LT_GID = 1_LB_INT
   else
      LB_LT_GID = 0_LB_INT
   endif
else if (a%int2 < b%int2) then
   LB_LT_GID = 1_LB_INT
else
   LB_LT_GID = 0_LB_INT
endif
end function LB_LT_GID

function LB_LE_GID(a,b)
integer(LB_INT) :: LB_LE_GID
type(LB_GID), intent(in) :: a,b
if (a%int2 == b%int2) then
   if (a%int1 <= b%int1) then
      LB_LE_GID = 1_LB_INT
   else
      LB_LE_GID = 0_LB_INT
   endif
else if (a%int2 <= b%int2) then
   LB_LE_GID = 1_LB_INT
else
   LB_LE_GID = 0_LB_INT
endif
end function LB_LE_GID

function LB_GT_GID(a,b)
integer(LB_INT) :: LB_GT_GID
type(LB_GID), intent(in) :: a,b
if (a%int2 == b%int2) then
   if (a%int1 > b%int1) then
      LB_GT_GID = 1_LB_INT
   else
      LB_GT_GID = 0_LB_INT
   endif
else if (a%int2 > b%int2) then
   LB_GT_GID = 1_LB_INT
else
   LB_GT_GID = 0_LB_INT
endif
end function LB_GT_GID

function LB_GE_GID(a,b)
integer(LB_INT) :: LB_GE_GID
type(LB_GID), intent(in) :: a,b
if (a%int2 == b%int2) then
   if (a%int1 >= b%int1) then
      LB_GE_GID = 1_LB_INT
   else
      LB_GE_GID = 0_LB_INT
   endif
else if (a%int2 >= b%int2) then
   LB_GE_GID = 1_LB_INT
else
   LB_GE_GID = 0_LB_INT
endif
end function LB_GE_GID

end module lb_user_const
