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

module lb_user_const
use zoltan_types
implicit none

! Note:  this file is provided for backward compatibility with prior versions
! of Zoltan (using the LB_* naming convention).
! New zoltan users should define their data types in zoltan_user_data.f90.
!
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

end module lb_user_const
