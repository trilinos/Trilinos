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
#include "permute.h"
#include "permute_mat.h"
#include "permute_rhs.h"
#include "forward.h"
#include "solve.h"
#include "mytime.h"
#include "exchange_pivots.h"


extern int me;	               /* processor id information */

extern int my_rows;            /* num of rows I own */
extern int my_cols;            /* num of cols I own */
extern int my_rhs;             /* num of right hand side I own */




void  X_PERMUTE_ (DATA_TYPE *mat, int *permutations)
{    

#ifdef TIMING0

   double ex_time,p_mat,t1;

#endif
 

  /* Exchange Pivoting Information 
       Now all processors have complete pivot information  */

#ifdef TIMING0

   t1 = MPI_Wtime();

#endif

  exchange_pivots(permutations);

#ifdef TIMING0

   ex_time  = MPI_Wtime()-t1;

   t1 = MPI_Wtime();
#endif

  permute_mat(mat,permutations);

#ifdef TIMING0

   p_mat  = MPI_Wtime()-t1;

   t1 = MPI_Wtime();

  showtime("Time to permute matrix",&p_mat);
  showtime("Time to do exchange pivots",&ex_time);

#endif

}
