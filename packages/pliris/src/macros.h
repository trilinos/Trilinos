/*
// @HEADER
// ***********************************************************************
// 
//                Pliris: Parallel Dense Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#define grey_c(P)     ((P)^((P)>>1))

#define lrow_to_grow(R) ( (mesh_row(me) + nprocs_col*(R))  )

#define grow_to_lrow(R) ( (R/nprocs_col)  )

/* #define col_owner(C)  (((C)%nprocs_row) + (me - me%nprocs_row)) */
#define col_owner(C)  ( proc_num(mesh_row(me) , (C)%nprocs_row) )

/* #define row_owner(R)  ((((R)%nprocs_col)*nprocs_row) + (me%nprocs_row)) */
#define row_owner(R)  ( proc_num((R)%nprocs_col , mesh_col(me)) )

#define owner(R, C)   ((((R)%nprocs_col)*nprocs_row) + ((C)%nprocs_row))

#define mesh_col(P)   ((P)%nprocs_row)

#define mesh_row(P)   ((P)/nprocs_row)

#define proc_num(R,C) ((R)*nprocs_row + (C))

#define mac_send_msg(D,B,S,T)  MPI_Send(B,S,MPI_CHAR,D,T,MPI_COMM_WORLD)
