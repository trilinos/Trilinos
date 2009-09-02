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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "defines.h"
#include "BLAS_prototypes.h"
#include "solve.h"
#include "macros.h"
#include "pcomm.h"

extern int me;

extern int nrows_matrix;       /* number of rows in the matrix */
extern int ncols_matrix;       /* number of cols in the matrix */

extern int nprocs_col;		/* num of procs to which a col is assigned */
extern int nprocs_row;		/* num of procs to which a row is assigned */

extern int my_first_col;        /* proc position in a col */
extern int my_first_row;	/* proc position in a row */

extern int my_rows;		/* num of rows I own */
extern int my_cols;            /* num of cols I own */

extern DATA_TYPE *col1;		/* ptrs to col used for updating a col */
extern int col1_stride;
extern int mat_stride;


extern DATA_TYPE *row1;		/* ptr to diagonal row */
extern DATA_TYPE *row2;		/* ptr to pivot row */

extern int nrhs;                /* number of right hand sides */
extern int my_rhs;              /* number of right hand sides that I own */

extern int colcnt;	        /* number of columns stored for BLAS 3 ops */
extern int blksz;	       	/* block size for BLAS 3 operations */
extern int rhs_blksz;

extern MPI_Comm col_comm;



#define SOSTATUSINT 32768

#define SOCOLTYPE (1<<26)
#define SOROWTYPE (1<<27)
#define SOHSTYPE (SOCOLTYPE + SOROWTYPE)
#define PERMTYPE ((1 << 25) + (1 << 26))

void
back_solve6(DATA_TYPE *fseg, DATA_TYPE *rhs)
{
  int  j;                      /* loop counters */

  int end_row;                  /* row num to end column operations */
  

  DATA_TYPE *ptr2;             /* running ptrs into mat and update cols */
  DATA_TYPE *ptr3, *ptr4;       /* running ptrs to diagonal and pivot rows */
  DATA_TYPE *ptr5;              /* running ptr into matrix */

  int bytes[16];                /* number of bytes in messages */
  int root;                     /* root processor for fanout */
  int type[16];                 /* mesage type for messages */
  int dest[16];                 /* dest for message sends */
 

  int one = 1;                  /* constant for the daxpy routines */

  DATA_TYPE d_one = CONST_ONE;  /* constant for the daxpy routines */
  DATA_TYPE d_min_one = CONST_MINUS_ONE; /* constant for the daxpy routines */

  int j1,j2;
 
  double diag_mag;              /* magnitude of matrix diagonal */

  int n_rhs_this;               /* num rhs that I currently own */


  char transA = 'N';            /* all dgemm's don't transpose matrix A  */
  char transB = 'N';            /* nor transpose matrix B */

  int col_offset;               /* which processor starts the pipeline */
  int my_pos;                   /* my position in the new linup */
  int extra;                    /* extra loop to realign data after pipeline */
  int act_col;                  /* act this column (that I own) */
  int on_col;                   /* on this collection of rhs's */
  int global_col;               /* global col number for act_col */
  int max_bytes;                /* max number of bytes of rhs I can receive */

  int my_col_id, my_row_id, id_temp;
  int dest_right, dest_left;

  int blas_length;

  MPI_Request msgrequest;
  MPI_Status msgstatus;

  /* set the block size for backsolve */

  my_col_id = mesh_col(me);
  my_row_id = mesh_row(me);

  id_temp = my_col_id + 1;
  if (id_temp >= nprocs_row) id_temp = 0;
  dest_right = proc_num(my_row_id,id_temp);

  id_temp = my_col_id - 1;
  if (id_temp < 0) id_temp = nprocs_row-1;
  dest_left = proc_num(my_row_id,id_temp);


  /* set j2 to be first column in last group of columns */
  rhs_blksz=1;
  max_bytes = nrhs/nprocs_row;
  if (nrhs%nprocs_row > 0) max_bytes++;
  max_bytes = max_bytes*sizeof(DATA_TYPE)*my_rows;

  row2 = (DATA_TYPE *) malloc( max_bytes);
   if (row2 == NULL) {
    fprintf(stderr, "Node %d: Out of memory\n", me);
    exit(-1);
   } 


  n_rhs_this = my_rhs;
  j2 = ncols_matrix-1;
  col_offset = (j2%nprocs_row);
  my_pos = my_first_col - col_offset;
  if (my_pos < 0) my_pos += nprocs_row;
  extra = (nprocs_row - (col_offset-1))%nprocs_row;
  
  act_col = my_cols-1;
  if (my_pos != 0) act_col++;

  on_col = my_pos;

  for (j = j2; j >= 1-nprocs_row-extra; j--) {
    if ((j+nprocs_row-1 >= 0) && (n_rhs_this > 0)) {

      if ((act_col < my_cols) && (act_col >= 0)) {

        global_col = act_col*nprocs_row + my_first_col;

        end_row = global_col/nprocs_col;
        if (my_first_row <= global_col%nprocs_col) ++end_row;

        ptr5 = fseg + act_col*my_rows;

        /* do an elimination step on the rhs that I own */

        ptr2 = ptr5;
        root = row_owner(global_col);

        if (me == root) {
          diag_mag = ABS_VAL(*(ptr2+end_row-1));
          for (j1=0; j1<n_rhs_this; j1++) {
            ptr3 = rhs + j1*my_rows + end_row - 1;
            ptr4 = row1 + j1;
            DIVIDE(*ptr3,*(ptr2+end_row-1),diag_mag,*ptr4);
            *ptr3 = *ptr4;
          }
          end_row--;
        }

        bytes[0] = n_rhs_this*sizeof(DATA_TYPE);
        type[0] = SOCOLTYPE+j;
        MPI_Bcast((char *) row1, bytes[0], MPI_CHAR, mesh_row(root), col_comm);

        /* added this barrier for CPLANT operation  */

        MPI_Barrier(col_comm);
 
       /* Changed XGEMM_ to XGEMM   removed all &   */

        XGEMM_(&transA, &transB, &end_row, &n_rhs_this, &one, &d_min_one, 
               ptr5, &my_rows,
               row1, &one, &d_one, 
               rhs, &my_rows);
      }
    }

    if (j != 1-nprocs_row-extra) {
      dest[0] = dest_right;
      if (me != dest[0]) {
	/*  Removed the follwoeing old code  */
/*        dest[1] = dest_right;
        bytes[1] = 0;
        type[1] = SOHSTYPE+j;

        MPI_Send((char *) rhs, bytes[1], MPI_CHAR, dest[1], type[1],MPI_COMM_WORLD);

        dest[1] = dest_left;
        bytes[1] = 0;
        type[1] = SOHSTYPE+j;

	 MPI_Recv((char *) rhs, bytes[1], MPI_CHAR, dest[1], type[1],MPI_COMM_WORLD,&msgstatus); 

         */      

	/*   dest[1] = dest_left;
        bytes[1] = n_rhs_this * sizeof(DATA_TYPE) * my_rows;
        type[1] = SOROWTYPE+j;
        MPI_Send((char *) rhs, bytes[1], MPI_CHAR, dest[1], type[1],MPI_COMM_WORLD); */

        bytes[0] = max_bytes;
        type[0] = SOROWTYPE+j;

        MPI_Irecv((char *) row2, bytes[0], MPI_CHAR, MPI_ANY_SOURCE,type[0],MPI_COMM_WORLD,&msgrequest);

	/*   MPI_Recv((char *)rhs,bytes[0],MPI_CHAR,MPI_ANY_SOURCE,type[0],MPI_COMM_WORLD,&msgstatus); */
        n_rhs_this = bytes[0]/sizeof(DATA_TYPE)/my_rows;

        dest[1] = dest_left;
        bytes[1] = n_rhs_this * sizeof(DATA_TYPE) * my_rows;
        type[1] = SOROWTYPE+j;
        MPI_Send((char *) rhs, bytes[1], MPI_CHAR, dest[1], type[1],MPI_COMM_WORLD);

        MPI_Wait(&msgrequest,&msgstatus);

        blas_length = n_rhs_this*my_rows;
        XCOPY(blas_length,row2,one,rhs,one);



      }
      on_col++;
      if (on_col >= nprocs_row) {
        on_col = 0;
        act_col--;
      }
    }
  }
  /* free(row2);  */
}

void collect_vector(DATA_TYPE *vec)
{
  int j, k;

  int start_col;
  int end_row;

  int bytes;
  int dest;
  int type;
  
  MPI_Status msgstatus;

  for (j=0; j<ncols_matrix; j++) {
    if (me == col_owner(j)) {
      start_col = (j) / nprocs_row;
      if (my_first_col < (j)%nprocs_row) ++start_col;

      if (j%rhs_blksz == 0) {
        for (k=0; k<rhs_blksz; k++) {
          end_row = (j+k)/nprocs_col;
          if (my_first_row <= (j+k)%nprocs_col) ++end_row;
          if (j+k < ncols_matrix) {
            if (me == row_owner(j+k)) {
              dest = col_owner(0);
              if (me != dest) {
                bytes = sizeof(DATA_TYPE);
                type = PERMTYPE + j + k;
                MPI_Send((vec+end_row-1),bytes,MPI_BYTE,
                   dest,type,MPI_COMM_WORLD);
              }
            }
          }
        }
      }
    }
    if (me == col_owner(0)) {
      if (me == row_owner(j)) {
        end_row = (j)/nprocs_col;
        if (my_first_row <= (j)%nprocs_col) ++end_row;
        dest = col_owner((j/rhs_blksz)*rhs_blksz);
        if (me != dest) {
          bytes = sizeof(DATA_TYPE);
          type = PERMTYPE + j;
          MPI_Recv((vec+end_row-1),bytes,MPI_BYTE,dest,
                   type,MPI_COMM_WORLD,&msgstatus);
        }
      }
    }
  }      
}
