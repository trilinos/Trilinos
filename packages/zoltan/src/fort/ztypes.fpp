#include 'sppr_header'
module zoltan_types

!--------------------------------------------------------------------------
! kind parameters to match Fortran types with C types.
! This might have to be adjusted on some machines.

integer, parameter :: &
   LB_INT = selected_int_kind(9), &
   LB_FLOAT = selected_real_kind(6), &
   LB_DOUBLE = selected_real_kind(15)

integer, parameter :: LB_PTR_LENGTH = 4

type LB_PTR
   character(len=LB_PTR_LENGTH) :: addr
end type LB_PTR

type(LB_PTR), parameter :: &
   LB_NULL_PTR = LB_PTR(char(0)//char(0)//char(0)//char(0))

interface operator(==)
   module procedure ptrcompare
end interface

private :: ptrcompare

contains

! comparison operator for C pointers

function ptrcompare(p1,p2)
logical :: ptrcompare
type(LB_PTR), intent(in) :: p1, p2
ptrcompare = (p1%addr == p2%addr)
end function ptrcompare

end module zoltan_types
