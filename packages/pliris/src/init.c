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
#include <mpi.h>
#include "defines.h"
#include "init.h"
#include "pcomm.h"
#include "BLAS_prototypes.h"
#include "dblassp.h"
#include "macros.h"

#define INITTYPE1 ((1 << 17) + 1)
#define INITTYPE2 ((1 << 17) + 2)
#define INITTYPE3 ((1 << 17) + 3)
#define INITTYPE4 ((1 << 17) + 4)
#define INITTYPE5 ((1 << 17) + 5)
#define INITTYPE6 ((1 << 17) + 6)

#define MATVECTYPE ((1 << 25) + (1 << 26))

#define RANDOM

/*
#define RANDOM
#define TOEPLITZ
#define RANDOM
#define ONEM
#define ONES
*/

extern int my_cols;            /* num of cols I own */
extern int my_cols_seg;        /* num of cols I own in a seg */
extern int my_cols_last;       /* num of cols I own in the rightmost seg */

extern int nsegs_row;          /* num of segs to which each row is assigned */
extern int ncols_seg;          /* num of cols in each segment */
extern int ncols_last;         /* num of cols in last segment */
extern int my_rows;            /* num of rows I own */
extern int nprocs_row;         /* num of procs to which a row is assigned */
extern int my_first_col;       /* proc position in a col */
extern int nprocs_col;         /* num of procs to which a col is assigned */
extern int my_first_row;       /* proc position in a row */
extern int nrows_matrix;       /* number of rows in the matrix */
extern int nrhs;               /* number of right hand sides */
extern int my_rhs;             /* number of right hand sides I own */
extern int me;                 /* who I am */
extern int rhs_blksz;          /* block size for backsolve */
extern DATA_TYPE *col1, *col2;
extern DATA_TYPE *row1, *row2;

extern MPI_Comm row_comm,col_comm;

/*     static int i_one =1;      */

void init_seg(DATA_TYPE *seg, int seg_num)
{

  int k, l;          /* loop counters */
  int my_cols_this;  /* num of cols I have in a seg */
  DATA_TYPE *tp;     /* temp ptr for stepping through the matrix entries */
  double drand48();  /* Rand number generator */ 
  long seedval;
 

  if (me == 0) {
#ifdef TOEPLITZ
    printf("USING TOEPLITZ MATRIX\n");
#endif
#ifdef RANDOM
    printf("Using random matrix \n");
#endif
#ifdef ONEM
    printf("USING ONEM MATRIX\n");
#endif
#ifdef ONES
    printf("USING ALL ONES MATRIX\n");
#endif
  }
  
  my_cols_this = my_cols;
    
  seedval = nrows_matrix + me;
  srand48(seedval);

  for ( l = 0; l < my_cols_this; l++) {
    tp = seg + l*(my_rows);
    for (k = 0 ; k < my_rows; k++) {
#ifdef COMPLEX
#ifdef TOEPLITZ
      (*tp).r = ( seg_num*ncols_seg + l*nprocs_row + my_first_col + 
                k*nprocs_col + my_first_row ) % nrows_matrix;
      (*tp).i = 0.0;
      tp++;
#endif
#ifdef RANDOM
      (*tp).r = 1.0 - 2.0*drand48();
      (*tp).i = 1.0 - 2.0*drand48();
      tp++;
#endif
#ifdef ONEM
      indj = l*nprocs_row + my_first_col + 1;
      indi = k*nprocs_col + my_first_row + 1;
      (*tp).r = (indj < indi) ? indj : indi;
      (*tp++).i = 0.0;
#endif
#ifdef ONES
      (*tp).r = 1.0;
      (*tp++).i = 0.0;
#endif

#else

#ifdef TOEPLITZ
      *tp++ = ( seg_num*ncols_seg + l*nprocs_row + my_first_col + 
                k*nprocs_col + my_first_row ) % nrows_matrix;
#endif
#ifdef RANDOM
      *tp++ = 1.0 - 2.0*drand48();
#endif
#ifdef ONEM
      indj = l*nprocs_row + my_first_col + 1;
      indi = k*nprocs_col + my_first_row + 1;
      *tp++ = (indj < indi) ? indj : indi;
#endif
#ifdef ONES
      *tp++ = 1.0;
#endif
#endif
    }
  }
}

void init_rhs(DATA_TYPE *rhs, DATA_TYPE *seg, int seg_num)
{

  int k, l;          /* loop counters */
  int stride;        /* stride for dasum routine */
  DATA_TYPE *tp,*tp1;   /* temp ptr for stepping through the matrix entries */
  int my_cols_this;  /* num of cols I have in a seg */

  my_cols_this = my_cols;

#ifdef COMPLEX
  stride = my_rows;
  for (k= 0; k < my_rows; k++) {
    tp = rhs + k;
    tp1 = seg + k;
    (*tp).r = 0.0;
    (*tp).i = 0.0;
    for (l=0; l < my_cols_this; l++) {
      (*tp).r += (*tp1).r;
      (*tp).i += (*tp1).i;
      tp1 += stride;
    }
    tp++;
  }
#ifdef DCPLX
  MPI_Allreduce((double *)rhs,(double *)col1,2*my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
#else 
MPI_Allreduce((float *)rhs,(float *)col1,2*my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
#endif
#ifdef ONES
  for (k = 0; k < my_rows; k++) {
    *(rhs+k) = nrows_matrix - (k*nprocs_col + my_first_row);
  }
#endif

#else

  stride = my_rows;
  for (k= 0; k < my_rows; k++) {
    tp = rhs + k;
    tp1 = seg + k;
    *tp = 0.0;
    for (l=0; l < my_cols_this; l++) {
      *tp += *tp1;
      tp1 += stride;
    }
    tp++;
  }
  MPI_Allreduce(rhs,col1,my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
#ifdef ONES
  for (k = 0; k < my_rows; k++) {
    *(rhs+k) = nrows_matrix - (k*nprocs_col + my_first_row);
  }
#endif
#endif
}

double one_norm(DATA_TYPE *seg, int seg_num)
{
  int i;
  int my_cols_this;
  DATA_TYPE *ptr1;
  double *ptr2;
  double local_max;
  int indmax, stride;
  int one=1;               /* constant for BLAS routine  */
 
  my_cols_this = my_cols;

  stride = my_rows;
  ptr1 = seg;
  ptr2 = (double *) row1;
  for (i=0; i<my_cols_this; i++) {
#ifdef CBLAS
	*ptr2 = XASUM(my_rows, ptr1, one);
#else
    	*ptr2 = XASUM(my_rows, ptr1, one);
#endif
    ptr1 += (my_rows);
    ptr2++;
  }

  MPI_Allreduce((double *) row1, (double *) row2,my_cols_this,MPI_DATA_TYPE,MPI_SUM,col_comm);
   indmax = IXAMAX(my_cols_this, row2, one);
#ifdef CBLAS
   /* do nothing index is ok    */

#else
   /*  Fortran to C adjust array index  */

   indmax = indmax -1;
#endif
  ptr2 = (double *) row2;
  local_max = *(ptr2+indmax);
  return max_all(local_max, INITTYPE3);
}

double inf_norm(DATA_TYPE *seg, int seg_num)
{
  int i;
  int my_cols_this;
  DATA_TYPE *ptr1;
  double *ptr2;
  double local_max;
  int indmax, stride;
  int one = 1 ;                  /* Constant for BLAS operation  */
  
  my_cols_this = my_cols;

  stride = my_rows;
  ptr1 = seg;
  ptr2 = (double *) col1;
  for (i=0; i<my_rows; i++) {
#ifdef CBLAS
	*ptr2 = XASUM(my_cols_last, ptr1, stride);
#else
    	*ptr2 = XASUM(my_cols_last, ptr1, stride);
#endif
    ptr1++;
    ptr2++;
  }

  MPI_Allreduce(col1,col2,my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
  indmax = IXAMAX(my_rows, col2, one);
#ifdef CBLAS
   /* do nothing index is ok    */

#else
   /*  Fortran to C adjust array index  */

   indmax = indmax -1;
#endif
  ptr2 = (double *) col2;
  local_max = *(ptr2+indmax);
  return max_all(local_max, INITTYPE6);
}

double init_eps(void)
{
  double tempa, tempb, tempc, eps;

  tempa= 4.0/3.0;
  for (eps = 0.0; eps == 0.0; ) {
    tempb= tempa - 1.0;
    tempc= tempb + tempb + tempb;
    eps = fabs(tempc-1.0);
  }
  return eps;
} 

void mat_vec(DATA_TYPE *seg, int seg_num, DATA_TYPE *vec) 
{
  int j;
  int k;

  int start_col;
  int end_row;
  int num_rows;    
  int stride;

  int root;
  int bytes;
  int dest;
  int type;
  int one = 1;       /*  Defined for BLAS call   */

  DATA_TYPE *ptr1;
  DATA_TYPE *ptr2;


  MPI_Status stat;
  extern MPI_Comm col_comm;

  rhs_blksz=1;
  /* collect vector in a row distribution */

  for (j=0; j<ncols_last; j++) {

    if (me == col_owner(j)) {
    
      start_col = (j) / nprocs_row;
      if (my_first_col < (j)%nprocs_row) ++start_col;

      num_rows = end_row;

      if (j%rhs_blksz == 0) {
        for (k=1; k<rhs_blksz; k++) {
          end_row = (j+k)/nprocs_col;
          if (my_first_row <= (j+k)%nprocs_col) ++end_row;
          if (j+k < ncols_last) {
            if (me == row_owner(j+k)) {
              dest = col_owner(j+k);
              if (me != dest) {
                bytes = sizeof(DATA_TYPE);
                type = MATVECTYPE + j + k;
                MPI_Send((char *) (vec+end_row-1),bytes, MPI_CHAR,dest,type,MPI_COMM_WORLD);
              }
            }
          }
        }
      } else {
        if (me == row_owner(j)) {
          end_row = (j)/nprocs_col;
          if (my_first_row <= (j)%nprocs_col) ++end_row;
          dest = col_owner((j/rhs_blksz)*rhs_blksz);
          if (me != dest) {
            bytes = sizeof(DATA_TYPE);
            type = MATVECTYPE + j;
            MPI_Recv((char *) (vec+end_row-1),bytes, MPI_CHAR,dest,type,MPI_COMM_WORLD,&stat);
          }
        }
      }

      end_row = (j)/nprocs_col;
      if (my_first_row <= (j)%nprocs_col) ++end_row;

      ptr1 = row1 + start_col;
      root = row_owner(j);
      if (me == root) {
        *ptr1 = *(vec + end_row - 1);
      }
      bytes = sizeof(DATA_TYPE);
      type = MATVECTYPE + j;
      MPI_Bcast((char *) ptr1, bytes, MPI_CHAR, mesh_row(root), col_comm);
    }
  }      

  MPI_Barrier(MPI_COMM_WORLD);

  ptr1 = vec;
  ptr2 = seg;
  stride = my_rows;
  for (j=0; j<my_rows; j++) {
#ifdef CBLAS
    *ptr1 = XDOT(my_cols_last, ptr2, stride, row1, one);
#else
     *ptr1 = XDOT(my_cols_last, ptr2, stride, row1, one);
#endif

    ptr1++;
    ptr2++;
  }

  MPI_Barrier(MPI_COMM_WORLD);
#ifdef COMPLEX
  MPI_Allreduce(vec,col1,2*my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
#else
  MPI_Allreduce(vec,col1,my_rows,MPI_DATA_TYPE,MPI_SUM,row_comm);
#endif
}  
