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

!/*--------------------------------------------------------------------------*/
!/* Purpose: Driver for dynamic load-balance library, ZOLTAN.                */
!/*                                                                          */
!/*--------------------------------------------------------------------------*/
!/* Author(s):  Matthew M. St.John (9226)                                    */
!   Translated to Fortran by William F. Mitchell
!/*--------------------------------------------------------------------------*/
!/*--------------------------------------------------------------------------*/
!/* Revision History:                                                        */
!/*                                                                          */
!/*    30 March 1999:    Date of creation                                    */
!       1 September 1999: Fortran translation
!/*--------------------------------------------------------------------------*/

!/****************************************************************************/
!/****************************************************************************/
!/****************************************************************************/


module dr_sort
use zoltan
implicit none
public :: dr_sort_index

contains

subroutine dr_sort_index(n, ra, indx)
use zoltan
integer(Zoltan_INT) :: n
integer(Zoltan_INT) :: ra(0:)
integer(Zoltan_INT) :: indx(0:)

!/*
!*       Numerical Recipies in C source code
!*       modified to have first argument an integer array
!*
!*       Sorts the array ra[0,..,(n-1)] in ascending numerical order using
!*       heapsort algorithm.
!*
!*/

  integer(Zoltan_INT) :: l, j, ir, i
  integer(Zoltan_INT) :: rra, irra
!  /*
!   *  No need to sort if one or fewer items.
!   */
  if (n <= 1) return

  l=n/2
  ir=n-1
  do
    if (l > 0) then
      l = l-1
      rra=ra(indx(l))
      irra=indx(l)
    else
      rra=ra(indx(ir))
      irra=indx(ir)

      indx(ir)=indx(0)
      ir = ir-1
      if (ir == 0) then
        indx(0)=irra
        return
      endif
    endif
    i=l
    j=2*l+1
    do while (j <= ir)
      if (j < ir .and. ra(indx(j)) < ra(indx(j+1))) j = j+1
      if (rra < ra(indx(j))) then
        indx(i)=indx(j)
        i = j
        j = j+i+1
      else
        j=ir+1
      endif
    end do
    indx(i)=irra
  end do
end subroutine dr_sort_index

end module dr_sort
