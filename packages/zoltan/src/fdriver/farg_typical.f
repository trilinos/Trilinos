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

! Command line argument functions for typical iargc, getarg implementations

      integer function mpir_iargc()
      mpir_iargc = iargc()
      return
      end

      subroutine mpir_getarg( i, s )
      integer       i
      character*(*) s
      call getarg(i,s)
      return
      end
