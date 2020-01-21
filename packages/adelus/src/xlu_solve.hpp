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

Author

Joseph D. Kotulski
Sandia National Labs
(505)-845-7955
jdkotul@sandia.gov


*/

#ifndef __XLUSOLVE_HPP__
#define __XLUSOLVE_HPP__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "mpi.h"
#include "vars.hpp"
#include "macros.h"
#include "block.h"
#include "solve.hpp"
#include "factor.hpp"
#include "perm1.hpp"
#include "pcomm.h"
#include "mytime.hpp"

#include "Kokkos_Core.hpp"

template<class ZDView>
void lusolve_(ZDView& ZV, int *matrix_size, int *num_procsr, int *num_rhs, double *secs)
{
  typedef typename ZDView::value_type value_type;
  typedef typename ZDView::device_type::execution_space execution_space;
  typedef typename ZDView::device_type::memory_space memory_space;

  //int begin_rhs;                /* Beginning index for the RHS   */
  double run_secs;              /* time (in secs) during which the prog ran */
  //double seconds();             /* function to generate timings */
  // Comment the line above because /home/projects/ppc64le/gcc/7.2.0/include/c++/7.2.0/chrono:600:37: note:   â€˜std::chrono::secondsâ€™
  // typedef duration<int64_t>       seconds;

  double tsecs;                 /* intermediate storage of timing info */

  int totmem;
/*
   Determine who I am (me ) and the total number of nodes (nprocs_cube)
                                                                        */
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs_cube);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  //mat = matrix;
  //rhs = rhsides;

  nrows_matrix = *matrix_size;
  ncols_matrix = *matrix_size;
  nprocs_row   = *num_procsr;

  totmem=0;                      /* Initialize the total memory used */
  nprocs_col = nprocs_cube/nprocs_row;
  max_procs = (nprocs_row < nprocs_col) ? nprocs_col : nprocs_row;

  /* set up communicators for rows and columns */

  myrow = mesh_row(me);
  mycol = mesh_col(me);

  MPI_Comm_split(MPI_COMM_WORLD,myrow,mycol,&row_comm);

  MPI_Comm_split(MPI_COMM_WORLD,mycol,myrow,&col_comm);

  /* Distribution for the matrix on me */

  my_first_col = mesh_col(me);
  my_first_row = mesh_row(me);

  my_rows = nrows_matrix / nprocs_col;
  if (my_first_row < nrows_matrix % nprocs_col)
    ++my_rows;
  my_cols = ncols_matrix / nprocs_row;
  if (my_first_col < ncols_matrix % nprocs_row)
    ++my_cols;

  /* blksz paramter must be set */

  blksz = DEFBLKSZ;

  /* Distribution for the rhs on me */

  nrhs = *num_rhs;
  my_rhs = nrhs / nprocs_row;
  if (my_first_col < nrhs % nprocs_row) ++my_rhs;

  /* Beginning position in array for the RHS   */

  //begin_rhs = my_cols * my_rows;

  printf("Rank %i -- lusolve_() Begin LU+Solve+Perm, value_type %s, execution_space %s, memory_space %s\n", me, typeid(value_type).name(), typeid(execution_space).name(), typeid(memory_space).name());

  /* allocate arrays for factor/solve */

  typedef Kokkos::View<value_type*,  Kokkos::LayoutLeft, memory_space> ViewType1D;
  typedef Kokkos::View<value_type**, Kokkos::LayoutLeft, memory_space> ViewType2D;
  typedef Kokkos::View<int*, Kokkos::LayoutLeft, memory_space> ViewIntType1D;

  //col1 = (DATA_TYPE *) malloc((blksz+1) * (my_rows + 1) * sizeof(DATA_TYPE));
  //row1 = (DATA_TYPE *) malloc(blksz*(my_cols+blksz+nrhs)* sizeof(DATA_TYPE));
  //row2 = (DATA_TYPE *) malloc((my_cols + blksz + nrhs)  * sizeof(DATA_TYPE));
  //row3 = (DATA_TYPE *) malloc((my_cols + blksz + nrhs)  * sizeof(DATA_TYPE));
  //pivot_vec = (int *) malloc(my_cols * sizeof(int));

  totmem += (blksz) * (my_rows) * sizeof(DATA_TYPE);     //col1_view
  totmem += blksz * (my_cols + blksz + nrhs) * sizeof(DATA_TYPE);//row1_view
  totmem += (my_cols + blksz + nrhs) * sizeof(DATA_TYPE);        //row2_view
  totmem += (my_cols + blksz + nrhs) * sizeof(DATA_TYPE);        //row3_view
  totmem += my_cols * sizeof(int);                               //pivot_vec_view
  
  //ViewType2D    col1_view      ( "col1_view",      my_rows + 1, blksz + 1 );
  ViewType2D    col1_view      ( "col1_view",      my_rows, blksz );
  ViewType2D    row1_view      ( "row1_view",      blksz, my_cols + blksz + nrhs );
  ViewType1D    row2_view      ( "row2_view",      my_cols + blksz + nrhs );
  ViewType1D    row3_view      ( "row3_view",      my_cols + blksz + nrhs );  
  ViewIntType1D pivot_vec_view ( "pivot_vec_view", my_cols );

  //row1_stride = my_cols+blksz+1;//Note: need to ask: why add 1 here? What happens when nrhs>1?
                                  //Note: comment out since row1_stride is not used
  col1_stride = my_rows;          //Note: need to ask: why not (my_rows+1) here?
  mat_stride  = my_rows;

  //if (pivot_vec == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  //
  //if (row3 == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  //
  //if (row2 == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  //
  //if (row1 == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  //
  //if (col2 == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  //
  //if (col1 == NULL) {
  //  fprintf(stderr, "Node %d: Out of memory\n", me);
  //  exit(-1);
  //}
  
  { // OLD Note: To avoid segmentation fault when XLU_SOLVE_ (containing unmanaged view) is called mulyiple times, it's safest to make sure unmanaged Views fall out of scope before freeing their memory.
  /* Wrap matrix and rhs raw pointers in unmanaged Views */
  	
  /* Factor and Solve the system */

  tsecs = get_seconds(0.0);

  initcomm();

  //factor(mat);
  printf("OpenMP or Cuda: Rank %i -- factor() starts ...\n", me);
  factor(ZV,
         col1_view,
         row1_view,
         row2_view, 
         row3_view, 
         pivot_vec_view);

  if (nrhs > 0) {

    /* Delele Unneccesary temp variables  */

    //free(row2);

    /* Perform the backsolve  */

    //back_solve6(mat, rhs);
    printf("OpenMP or Cuda: Rank %i -- back_solve6() starts ...\n", me);
    back_solve6(ZV);

    /* Permute the results -- undo the torus map    */

    //perm1_((mat+begin_rhs),&my_rhs);
    printf("OpenMP or Cuda: Rank %i -- perm1_()(permute the results -- undo the torus map) starts ...\n", me);
    auto sub_ZV = subview(ZV, Kokkos::ALL(), Kokkos::make_pair(my_cols, my_cols + my_rhs + 6));
    perm1_(sub_ZV, &my_rhs);
  }

  tsecs = get_seconds(tsecs);

  run_secs = (double) tsecs;

  /* Solve time secs */

  *secs = run_secs;
  showtime("Total time in Factor and Solve",&run_secs);

  } // Note: To avoid segmentation fault when XLU_SOLVE_ (containing unmanaged view) is called mulyiple times, it's safest to make sure unmanaged Views fall out of scope before freeing their memory

  //free(col1);
  //free(col2);
  //free(row1);

  //free(row3);
  //free(pivot_vec);

}

#endif
