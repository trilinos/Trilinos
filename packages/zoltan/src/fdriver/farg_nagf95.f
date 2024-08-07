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

! Command line argument functions for NAGWare f95 4.0

      integer function mpir_iargc()
      use f90_unix_env
      mpir_iargc = iargc()
      return
      end

      subroutine mpir_getarg( i, s )
      use f90_unix_env
      integer       i
      character*(*) s
      integer lenarg, ierr
      call getarg(i,s,lenarg,ierr)
      return
      end

