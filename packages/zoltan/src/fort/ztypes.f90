!! 
!! @HEADER
!!
!!!!**********************************************************************
!!
!!  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
!!                  Copyright 2012 Sandia Corporation
!!
!! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
!! the U.S. Government retains certain rights in this software.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are
!! met:
!!
!! 1. Redistributions of source code must retain the above copyright
!! notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright
!! notice, this list of conditions and the following disclaimer in the
!! documentation and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the Corporation nor the names of the
!! contributors may be used to endorse or promote products derived from
!! this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
!! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
!! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! Questions? Contact Karen Devine	kddevin@sandia.gov
!!                    Erik Boman	egboman@sandia.gov
!!
!!!!**********************************************************************
!!
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
