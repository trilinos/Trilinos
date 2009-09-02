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
