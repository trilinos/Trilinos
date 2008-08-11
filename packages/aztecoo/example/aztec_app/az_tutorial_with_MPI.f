C @HEADER
C ************************************************************************
C 
C        AztecOO: An Object-Oriented Aztec Linear Solver Package 
C                 Copyright (2002) Sandia Corporation
C 
C Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
C license for use of this work by or on behalf of the U.S. Government.
C 
C This library is free software; you can redistribute it and/or modify
C it under the terms of the GNU Lesser General Public License as
C published by the Free Software Foundation; either version 2.1 of the
C License, or (at your option) any later version.
C  
C This library is distributed in the hope that it will be useful, but
C WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
C Lesser General Public License for more details.
C  
C You should have received a copy of the GNU Lesser General Public
C License along with this library; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
C USA
C Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
C 
C ************************************************************************
C @HEADER
C====================================================================
C ------------------------
C | CVS File Information |
C ------------------------
C
C $RCSfile$
C
C $Author$
C
C $Date$
C
C $Revision$
C
C $Name$
C====================================================================*/

C
C
       program main
C
C---------------------------------------------------------------
C      Set up a 2D Poisson test problem and solve it with AZTEC.
C      Author:   Ray Tuminaro, Div 1422, Sandia National Labs
C      date:     11/10/94
C
       implicit none
       include "az_aztecf.h"
       include "mpif.h"
C
       integer   n, nrow
       common    /global/ n
C           POISSON EQUATION WILL BE SOLVED ON an n x n GRID.
C	    NOTE: n should be odd for rhs to be properly set.
C
       double precision b(0:1024),x(0:1024)
C           rhs and approximate solution
       integer    i
C
C             See Aztec User's Guide for the variables that follow:
C
       integer proc_config(0:AZ_PROC_SIZE), options(0:AZ_OPTIONS_SIZE)
       double precision params(0:AZ_PARAMS_SIZE)
       integer data_org(0:1024)
       double precision status(0:AZ_STATUS_SIZE)
       integer update(0:1024), external(0:1024)
       integer update_index(0:1024), extern_index(0:1024)
       integer bindx(0:1024)
       double  precision val(0:1024)
       integer N_update,ierror
C
C
       call MPI_INIT(ierror)
C
C           # of unknowns updated on this node
	n = 6
C
C      get number of processors and the name of this processor
C
       call AZ_set_proc_config(proc_config, MPI_COMM_WORLD)
C
C      Define paritioning:matrix rows in ascending order assigned
C      to this node
C
       nrow = n*n
       call AZ_read_update(N_update,update,proc_config,nrow,1,0)
C
C      create the matrix: each processor creates only rows
C      appearing in update[] (using global col. numbers).
C
       bindx(0) = N_update+1
       do 250 i = 0, N_update-1
          call create_matrix_row_5pt(update(i),i,val,bindx)
250    continue
C
C      convert matrix to a local distributed matrix */
C
       call AZ_transform(proc_config,external,bindx,val,update,
     $                   update_index,extern_index,data_org,
     $                   N_update,0,0,0,0,AZ_MSR_MATRIX)
C
C      initialize AZTEC options
C
       call AZ_defaults(options, params)
C
C      Set rhs (delta function at grid center) and initialize guess
C
       do 350 i = 0, N_update-1
          x(update_index(i)) = 0.0
          b(update_index(i)) = 0.0
          if (update(i) .eq. 0) b(update_index(i)) = 1.0
350    continue
C
C      solve the system of equations using b  as the right hand side
C
       call AZ_solve(x,b, options, params, 0,bindx,0,0,
     $               0,val, data_org, status, proc_config)

C
       call MPI_FINALIZE(ierror)
C
       stop
       end

C*********************************************************************
C*********************************************************************
C
       subroutine create_matrix_row_5pt(row,location,val,bindx)
C
       integer row,location,bindx(0:*)
       double precision val(0:*)
       integer   n
       common    /global/ n
C
C Add one row to an MSR matrix corresponding to a 5pt discrete
C approximation to the 2D Poisson operator on an n x n square.
C
C  Parameters:
C     row          == global row number of the new row to be added.
C     location     == local row where diagonal of the new row will be stored.
C     val,bindx    == (see user's guide). On output, val[] and bindx[]
C                     are appended such that the new row has been added.
C
       integer k
C
C      check neighbors in each direction and add nonzero if neighbor exits
C
       k = bindx(location)
       bindx(k)  = row + 1
       if (mod(row,n) .ne. n-1) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row - 1
       if (mod(row,n) .ne.   0) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row + n
       if (mod(row/n,n) .ne. n-1) then
          val(k) = -1.
          k = k + 1
       endif
       bindx(k)  = row - n
       if (mod(row/n,n) .ne.   0) then
          val(k) = -1.
          k = k + 1
       endif

       bindx(location+1) = k
       val(location)   = 4.
       return
       end
