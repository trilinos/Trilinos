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

Authors:

Brian Driessen
Sandia National Labs
(505)-844-9297
bjdries@sandia.gov

Joseph D. Kotulski
Sandia National Labs
(505)-845-7955
jdkotul@sandia.gov

*/


#include <math.h>
#include <stdio.h>

#include "mpi.h"
#include "defines.h"
#include "x_solve.h"
#include "permute_mat.h"
#include "permute_rhs.h"
#include "forward.h"
#include "solve.h"
#include "macros.h"
#include "exchange_pivots.h"
#include "perm1.h"
#include "mytime.h"




extern int me;	               /* processor id information */

extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */
extern int my_rhs;             /* num of right hand side I own */
extern int matrix_size;        /* Size of the Matrix  */
extern int nrhs;
extern int nprocs_row;

void  X_SOLVE_ (DATA_TYPE *mat, int *permutations,
        DATA_TYPE *rhs,int *num_rhs)
{
  int begin_rhs;

  int my_first_col;

#ifdef TIMING0

   double permtime, f_solve,b_solve, t1;

#endif
    MPI_Comm_rank(MPI_COMM_WORLD, &me);


    my_first_col = mesh_col(me);

  /* Distribution for the rhs on me */

  nrhs = *num_rhs;
  my_rhs = nrhs / nprocs_row;
  if (my_first_col < nrhs % nprocs_row) ++my_rhs;

  begin_rhs = my_cols*my_rows;



  /* Exchange Pivoting Information
       Now all processors have complete pivot information  */

#ifdef TIMING0

   t1 = MPI_Wtime();

#endif

  permute_rhs(rhs,permutations);

#ifdef TIMING0

   permtime  = MPI_Wtime() - t1;

#endif

#ifdef TIMING0

   t1 = MPI_Wtime();

#endif

  /* Forward Solve  */
  forward(mat, rhs);

#ifdef TIMING0

   f_solve = MPI_Wtime() -t1;

#endif


  MPI_Barrier(MPI_COMM_WORLD);

#ifdef TIMING0

   t1 = MPI_Wtime();

#endif

  /* Back Solve  */
  back_solve6(mat,rhs);

#ifdef TIMING0

   b_solve = MPI_Wtime() -t1;

#endif

  MPI_Barrier(MPI_COMM_WORLD);

  /* Permute the answer -- torus map inverse */

  perm1_((rhs),&my_rhs);

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef TIMING0

  showtime("Time to permute matrix",&permtime);
  showtime("Time to do forward solve",&f_solve);
  showtime("Time to do back solve",&b_solve);

#endif

}
