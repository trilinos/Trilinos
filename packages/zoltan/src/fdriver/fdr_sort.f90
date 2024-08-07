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

!--------------------------------------------------------------------------
! Purpose: Driver for dynamic load-balance library, ZOLTAN.                
!                                                                          
!--------------------------------------------------------------------------
! Author(s):  Matthew M. St.John (9226)                                    
!   Translated to Fortran by William F. Mitchell
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Revision History:                                                        
!                                                                          
!    30 March 1999:    Date of creation                                    
!       1 September 1999: Fortran translation
!--------------------------------------------------------------------------

!**************************************************************************
!**************************************************************************
!**************************************************************************


module dr_sort
use zoltan
implicit none
public :: dr_sort_index
public :: dr_sort2_index

contains

subroutine dr_sort_index_sub(sorted, val1, starti, endi, equal, larger)
use zoltan
integer(Zoltan_INT) :: starti, endi, equal, larger
integer(Zoltan_INT) :: sorted(0:)
integer(Zoltan_INT) :: val1(0:)
integer(Zoltan_INT) :: i, key, next, key_next

  i = (endi + starti) / 2
  key = val1(sorted(i))

  equal = starti
  larger = starti
  do i = starti, endi
     next = sorted(i)
     key_next = val1(next)
     if (key_next < key) then
        sorted(i) = sorted(larger)
        sorted(larger) = sorted(equal)
        larger = larger + 1
        sorted(equal)  = next
        equal = equal + 1
     else 
       if (key_next == key) then
        sorted(i) = sorted(larger)
        sorted(larger) = next
        larger = larger + 1
       endif
     endif
  end do
end subroutine dr_sort_index_sub

recursive subroutine dr_sort_index(starti, endi, ra, indx)
use zoltan
integer(Zoltan_INT) :: starti, endi
integer(Zoltan_INT) :: ra(0:)
integer(Zoltan_INT) :: indx(0:)

integer(Zoltan_INT) :: equal, larger

  if (starti < endi) then
     call dr_sort_index_sub(indx,ra,starti,endi,equal,larger)
     call dr_sort_index(starti, equal-1, ra, indx)
     call dr_sort_index(larger, endi, ra, indx)
  endif

end subroutine dr_sort_index

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dr_sort2_index_sub(sorted, val1, val2, starti, endi, equal, larger)
use zoltan
integer(Zoltan_INT) :: starti, endi, equal, larger
integer(Zoltan_INT) :: sorted(0:)
integer(Zoltan_INT) :: val1(0:)
integer(Zoltan_INT) :: val2(0:)
integer(Zoltan_INT) :: i, key1, next, key1_next, key2, key2_next

  i = (endi + starti) / 2
  key1 = val1(sorted(i))
  key2 = val2(sorted(i))

  equal = starti
  larger = starti
  do i = starti, endi
     next = sorted(i)
     key1_next = val1(next)
     key2_next = val2(next)
     if ((key1_next < key1) .or. ((key1_next == key1) .and. (key2_next < key2))) then
        sorted(i) = sorted(larger)
        sorted(larger) = sorted(equal)
        larger = larger + 1
        sorted(equal)  = next
        equal = equal + 1
     else 
       if ((key1_next == key1) .and. (key2_next == key2)) then
        sorted(i) = sorted(larger)
        sorted(larger) = next
        larger = larger + 1
       endif
     endif
  end do
end subroutine dr_sort2_index_sub

recursive subroutine dr_sort2_index(starti, endi, val1, val2, indx)
use zoltan
integer(Zoltan_INT) :: starti, endi
integer(Zoltan_INT) :: val1(0:)
integer(Zoltan_INT) :: val2(0:)
integer(Zoltan_INT) :: indx(0:)

integer(Zoltan_INT) :: equal, larger

  if (starti < endi) then
     call dr_sort2_index_sub(indx,val1,val2,starti,endi,equal,larger)
     call dr_sort2_index(starti, equal-1, val1, val2, indx)
     call dr_sort2_index(larger, endi, val1, val2, indx)
  endif

end subroutine dr_sort2_index

end module dr_sort
