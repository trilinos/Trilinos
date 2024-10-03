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

module zoltan_user_data
use zoltan_types
implicit none


! User defined data types for passing data to the query functions.  These can
! be used any way you want, but one suggestion is to use them as "wrapper"
! types for your own user defined types, e.g.
! type Zoltan_User_Data_1 
!    type(your_type), pointer :: ptr
! end type
! Exactly four data types must be defined, but you don't have to use them.

type Zoltan_User_Data_1
   integer :: dummy
end type Zoltan_User_Data_1

type Zoltan_User_Data_2
   integer :: dummy
end type Zoltan_User_Data_2

type Zoltan_User_Data_3
   integer :: dummy
end type Zoltan_User_Data_3

type Zoltan_User_Data_4
   integer :: dummy
end type Zoltan_User_Data_4

end module zoltan_user_data
