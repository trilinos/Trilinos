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
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "BLAS_prototypes.h"
double  timing(double secs, int type);
#include "xlu_solve.h"
#include <mpi.h>
#include "vars.h"
#include "macros.h"
#include "block.h"
#include "solve.h"
#include "factor.h"
#include "pcomm.h"

#define PERMTYPE ((1 << 5) + (1 << 4))


void XLU_SOLVE_ (DATA_TYPE *matrix, int *matrix_size, int *num_procsr, 
    DATA_TYPE *rhsides, int *num_rhs, double *secs)
{
 
  DATA_TYPE *mat;
 
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


  row2 = (DATA_TYPE *) malloc((my_cols + blksz + nrhs) * sizeof(DATA_TYPE));
  totmem += (my_cols + blksz + 1) * sizeof(DATA_TYPE);
  if (row2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
  }

  row1_stride = my_cols+blksz+1;
  row1 = (DATA_TYPE *) malloc(blksz*(my_cols+blksz+nrhs)*sizeof(DATA_TYPE));
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
  if (nrhs > 0) {
     if(me == 0 ) printf(" Entering backsolve \n");
      free(row2);
      back_solve6(mat, rhs);

 /* Permute the results -- undo the torus map    */
   if(me == 0 ) printf(" Entering perm \n");
    perm1_((mat+begin_rhs),&my_rhs);    
      
  }
  tsecs = seconds(tsecs);
  run_secs = (double) tsecs;

  /* Solve time secs */

  *secs = run_secs;
  
  free(col1);
  free(col2);
  free(row1);
 
  free(row3);
  free(pivot_vec);
}


/*  Permutes -- unwraps the torus-wrap for the solution  

    using the communication buffer             */


void perm1_(DATA_TYPE *vec, int *num_my_rhs)
{
  int i;
  int one=1;
  int my_rhs;

 
  int bytes;
  int dest;
  int type;

  int global_index;
  int local_index;
  

  int col_offset, row_offset;
  int ncols_proc1, ncols_proc2, nprocs_row1;
  int nrows_proc1, nrows_proc2, nprocs_col1;

  int change_count;
  int change_nosend;
  int change_send;
  int next, next_s;
  int num;
  int inc;

  MPI_Request msgrequest;
  MPI_Status msgstatus;

 
  DATA_TYPE *ptr1;
 
  DATA_TYPE *temp_s;
  DATA_TYPE *rhs_temp;

  
   my_rhs=*num_my_rhs;

   temp_s = (DATA_TYPE *) malloc((my_rhs + 1) 
                                    * sizeof(DATA_TYPE));
    if (temp_s == NULL) {
	fprintf(stderr, "Node %d: Out of memory for temp_s vector\n", me);
	exit(-1);
    }


    rhs_temp = (DATA_TYPE *) malloc((my_rows*(my_rhs+1)) * sizeof(DATA_TYPE));

    if (rhs_temp == NULL) {
	fprintf(stderr, "Node %d: Out of memory for rhs_temp vector\n", me);
	exit(-1);
    }



  


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
    next = 0;
    change_send = 0;
    next_s = 0; 
    num = my_rhs*my_rows;

    /* XCOPY(num,(vec),one,rhs_temp,one);  */
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

        /* Check if I need to send the data or just change position */

        if( dest == me ) {

         XCOPY(my_rhs, ptr1, my_rows, (row1 + next), one);
#ifdef COMPLEX
         (*(row1+ next + my_rhs)).r = local_index + 0.1;
#else
         *(row1+ next + my_rhs) = local_index + 0.1;
#endif
          change_nosend++;

          next = change_nosend * (my_rhs + 1);
    
	}

        if( dest !=me ) {
         
          bytes = (my_rhs + 1)*sizeof(DATA_TYPE);

          MPI_Irecv( (char *)(rhs_temp +next_s),bytes,MPI_CHAR,MPI_ANY_SOURCE,
                MPI_ANY_TAG,MPI_COMM_WORLD,&msgrequest);

         XCOPY(my_rhs, ptr1, my_rows, temp_s, one);
#ifdef COMPLEX
         (*(temp_s+my_rhs)).r = local_index + 0.1;
#else
         *(temp_s+my_rhs) = local_index + 0.1;
#endif
       
         type = PERMTYPE+change_send;
         MPI_Send((char *) temp_s,bytes,MPI_CHAR,dest,
                 type,MPI_COMM_WORLD);
         change_send++;
        
         next_s = change_send * (my_rhs+1);

         MPI_Wait(&msgrequest,&msgstatus);

	}

      }
      ptr1++;
    }

    /* Unpack changes from other processors  */

    next_s = 0;
    inc = 0;
    for (i = 0; i < change_send; i++) {
#ifdef COMPLEX
      local_index = (*(rhs_temp+next_s+my_rhs)).r;
#else
      local_index = *(rhs_temp +next_s+my_rhs);
#endif
      XCOPY(my_rhs,(rhs_temp + next_s),one,(vec+local_index),my_rows);
      
       inc++;
       next_s = inc * (my_rhs+1);
       }  

    /* Unpack my changes */
    next = 0;
    inc = 0;
    for (i = 0; i < change_nosend; i++) {
#ifdef COMPLEX
      local_index = (*(row1+next+my_rhs)).r;
#else
      local_index = *(row1+next+my_rhs);
#endif
      XCOPY(my_rhs,(row1 + next),one,(vec+local_index),my_rows);
       inc++;
       next = inc * (my_rhs+1);
    }
    
    
  }
  
  free(rhs_temp);
  free(temp_s);
} 






