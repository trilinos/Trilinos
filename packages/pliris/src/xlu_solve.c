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

Author

Joseph D. Kotulski
Sandia National Labs
(505)-845-7955
jdkotul@sandia.gov


*/





#include <stdio.h>
#include <math.h>
#include "defines.h"
#include "BLAS_prototypes.h"
double  timing(double secs, int type);
#include "xlu_solve.h"
#include "mpi.h"
#include "vars.h"
#include "macros.h"
#include "block.h"
#include "perm1.h"


#define PERMTYPE ((1 << 5) + (1 << 4))


void XLU_SOLVE_ (DATA_TYPE *matrix, int *matrix_size, int *num_procsr, 
    DATA_TYPE *rhsides, int *num_rhs, double *secs)
{
  int i, j, k, l;
  DATA_TYPE *mat;
  int MSPLIT;
  int begin_rhs;                /* Beginning index for the RHS   */
  double run_secs;              /* time (in secs) during which the prog ran */
  double seconds();             /* function to generate timings */
  double tsecs;                 /* intermediate storage of timing info */

  int totmem;
/*
   Determine who I am (me ) and the total number of nodes (nprocs_cube)
                                                                        */

  MPI_Comm_size(MPI_COMM_WORLD,&nprocs_cube);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);


  mat = matrix;
  rhs = rhsides;
  nrows_matrix = *matrix_size;
  ncols_matrix = *matrix_size;
  nprocs_row = *num_procsr;

  totmem=0;                      /* Initialize the total memory used */
  nprocs_col = nprocs_cube/nprocs_row;
  max_procs = (nprocs_row < nprocs_col) ? nprocs_col : nprocs_row;

    /* set up communicators for rows and columns */

    myrow = mesh_row(me);
    mycol = mesh_col(me);
    MPI_Comm_split(MPI_COMM_WORLD,myrow,mycol,&row_comm);
    MPI_Comm_split(MPI_COMM_WORLD,mycol,myrow,&col_comm);

    {int checkcol,checkrow;
     MPI_Comm_rank(col_comm, &checkrow) ;
     MPI_Comm_rank(row_comm, &checkcol) ;
     if (myrow != checkrow)
       printf("Node %d: my row = %d but rank in col = %d\n",me,myrow,checkrow);     if (mycol != checkcol)
       printf("Node %d: my col = %d but rank in row = %d\n",me,mycol,checkcol);   }

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

  begin_rhs = my_cols * my_rows; 

  /* allocate arrays for factor/solve */



  pivot_vec = (int *) malloc(my_cols * sizeof(int));
  totmem += my_cols * sizeof(int);
  if (pivot_vec == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  row3 = (DATA_TYPE *) malloc((my_cols + blksz + nrhs) * sizeof(DATA_TYPE));
  totmem += (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row3 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }
  rhs_temp = (DATA_TYPE *) malloc((my_rows*my_rhs+1) * sizeof(DATA_TYPE));

  row2 = (DATA_TYPE *) malloc(2*(my_cols + blksz + nrhs) * sizeof(DATA_TYPE));
  totmem += (my_cols + blksz + nrhs) * sizeof(DATA_TYPE);
  if (row2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  row1_stride = my_cols+blksz+1;
  row1 = (DATA_TYPE *) malloc(blksz*(my_cols+blksz+nrhs+10)*sizeof(DATA_TYPE));
  totmem += blksz * (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row1 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  col2 = (DATA_TYPE *) malloc((my_rows + 1) * sizeof(DATA_TYPE));
  totmem += (my_rows + 1) * sizeof(DATA_TYPE);
  if (col2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  col1_stride = my_rows;
  col1 = (DATA_TYPE *) malloc((blksz+1) * (my_rows + 1) * sizeof(DATA_TYPE));
  totmem += blksz * (my_rows + 1) * sizeof(DATA_TYPE);
  if (col1 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  mat_stride = my_rows;

  /* Factor and Solve the system */

  tsecs = seconds(0.0);
  initcomm();
  factor(mat);
  printf( " %d leaving factor\n",me);

  if (nrhs > 0) {
     back_solve6(mat, rhs);
     printf( " %d leaving solve\n",me);
 /* Permute the results -- undo the torus map    */
   perm1_((mat+begin_rhs),&my_rhs);   
      
  }
  tsecs = seconds(tsecs);
  run_secs = (double) tsecs;

  /* Solve time secs */

  *secs = run_secs;
  free(rhs_temp);
  free(col1);
  free(col2);
  free(row1);
  free(row2);
  free(row3);
  free(pivot_vec);
}


/*  Permutes -- unwraps the torus-wrap for the solution  

    using the communication buffer             */


void perm1_(DATA_TYPE *vec, int *num_my_rhs)
{
  int i;
  int j;
  int k;

  int start_col;
  int end_row;
  int num_rows;
  int stride;
  int one=1;
  int my_rhs;

  int root;
  int bytes;
  int dest;
  int type;

  int global_index;
  int local_index;
  int my_index;

  int col_offset, row_offset;
  int ncols_proc1, ncols_proc2, nprocs_row1;
  int nrows_proc1, nrows_proc2, nprocs_col1;

  int change_count;
  int change_nosend;
  int num;

  MPI_Request msgrequest;
  MPI_Status msgstatus;

  DATA_TYPE buf[2];
  DATA_TYPE *ptr1;
  DATA_TYPE *ptr2;


  my_rhs=*num_my_rhs;


  if (my_rhs > 0) {
    ncols_proc1 = ncols_matrix/nprocs_row;
    ncols_proc2 = ncols_proc1;
    if (ncols_matrix%nprocs_row > 0) ncols_proc1++;
    nprocs_row1 = (ncols_matrix%nprocs_row);
    row_offset = ncols_proc1 * nprocs_row1;

    nrows_proc1 = nrows_matrix/nprocs_col;
    nrows_proc2 = nrows_proc1;
    if (nrows_matrix%nprocs_col > 0) nrows_proc1++;
    nprocs_col1 = (nrows_matrix%nprocs_col);
    col_offset = nrows_proc1 * nprocs_col1;

    ptr1 = vec;
    change_count = 0;
    change_nosend = 0;
    num = my_rhs*my_rows;

    XCOPY(num,(vec),one,rhs_temp,one);
    for (i=0; i<my_rows; i++) {
      global_index = my_first_row + i*nprocs_col;
      /* break down global index using torus wrap in row direction */
      dest = global_index%nprocs_row;
      local_index = global_index/nprocs_row;

      /* rebuild global index using block in row direction */
      if (dest < nprocs_row1) {
        global_index = dest*ncols_proc1 + local_index;
      } else {
        global_index = row_offset + (dest-nprocs_row1)*ncols_proc2
                       + local_index;
      }

      /* break down global index using blocks in the column direction */
      if (global_index < col_offset) {
        dest = global_index/nrows_proc1;
        local_index = global_index%nrows_proc1;
      } else {
        dest = (global_index - col_offset)/nrows_proc2 + nprocs_col1;
        local_index = (global_index - col_offset)%nrows_proc2;
      }

      dest = dest*nprocs_row + my_first_col;
      if ((local_index != i) || (dest != me)) {

        XCOPY(my_rhs, ptr1, my_rows, row1, one);
#ifdef COMPLEX
        (*(row1+my_rhs)).r = local_index + 0.1;
#else
        *(row1+my_rhs) = local_index + 0.1;
#endif
        bytes = (my_rhs + 1)*sizeof(DATA_TYPE);
        type = PERMTYPE+change_count;
        MPI_Send((char *) row1,bytes,MPI_CHAR,dest,
                 type,MPI_COMM_WORLD);
        change_count++;

      }
      ptr1++;
    }
    for (i=0; i< change_count; i++) {
      bytes = (my_rhs + 1)*sizeof(DATA_TYPE);
      type = PERMTYPE+i;
      dest = -1;
      MPI_Recv((char *) row1,bytes,MPI_CHAR,MPI_ANY_SOURCE,
                MPI_ANY_TAG,MPI_COMM_WORLD,&msgstatus);
#ifdef COMPLEX
      local_index = (*(row1+my_rhs)).r;
#else
      local_index = *(row1+my_rhs);
#endif
      XCOPY(my_rhs,row1,one,(vec+local_index),my_rows);

    }
  }
} 

