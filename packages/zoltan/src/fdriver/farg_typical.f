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
