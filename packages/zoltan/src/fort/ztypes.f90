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

module zoltan_types

!--------------------------------------------------------------------------
! kind parameters to match Fortran types with C types.
! This might have to be adjusted on some machines.

integer, parameter :: &
   Zoltan_INT = selected_int_kind(9), &
!  Always assume 64-bit pointers.  If sizeof(void*) < 8, extra bytes are
!  padded to zero in Zfw_Create.
   Zoltan_INT_PTR = selected_int_kind(17), &
   Zoltan_FLOAT = selected_real_kind(6), &
   Zoltan_DOUBLE = selected_real_kind(15)


! Always assume 64-bit pointers.  If sizeof(void*) < 8, extra bytes are
! padded to zero in Zfw_Create.
integer, parameter :: Zoltan_PTR_LENGTH = 8

type Zoltan_PTR
   sequence
   character(len=Zoltan_PTR_LENGTH) :: addr
end type Zoltan_PTR

type(Zoltan_PTR), parameter :: &
   Zoltan_NULL_PTR = Zoltan_PTR( &
                 char(0)//char(0)//char(0)//char(0)// &
                 char(0)//char(0)//char(0)//char(0))

!--------------------------------------------------------------------------
! user defined types corresponding to the C structs

type Zoltan_Struct
   sequence
   type(Zoltan_PTR) :: addr
!#ifdef ABSOFT
! workaround for a bug in the Absoft compiler
!   integer :: dummy
!#endif
end type Zoltan_Struct

interface operator(==)
   module procedure ptrcompare
end interface

!  Include this section for backward compatibility with old Zoltan interface.
integer, parameter :: &
   LB_INT = Zoltan_INT, &
   LB_INT_PTR = Zoltan_INT_PTR, &
   LB_FLOAT = Zoltan_FLOAT, &
   LB_DOUBLE = Zoltan_DOUBLE

private :: ptrcompare

contains

! comparison operator for C pointers

function ptrcompare(p1,p2)
logical :: ptrcompare
type(Zoltan_PTR), intent(in) :: p1, p2
integer :: i
! should be able to compare the strings with a single statement
!   ptrcompare = (p1%addr == p2%addr)
! but bugs in PGI 3.1-2 and NAG 4.0 require comparing one character
! at a time
ptrcompare = .true.
do i=1,Zoltan_PTR_LENGTH
   ptrcompare = ptrcompare .and. (p1%addr(i:i) == p2%addr(i:i))
end do
end function ptrcompare

end module zoltan_types
