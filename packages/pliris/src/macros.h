/*
//@HEADER
// ************************************************************************
//
//               Pliris: Parallel Dense Solver Package
//                 Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
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

#if defined(Pliris_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Pliris package is deprecated"
#endif
#endif
