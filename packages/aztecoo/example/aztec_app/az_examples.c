/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char rcsid[] = "$Id$";
#endif

/******************************************************************************
 *
 * These routines correspond to the matrix setup phase for a number of sample
 * PDE problems using AZTEC. All of the examples are discretizations of
 * the Poisson equation (2D structured using 5pt operator, 2D structured using
 * 9pt operator, 2D unstructured using finite elements, 3D structured using
 * 7pt operator). All the examples are for MSR matrices excluding one VBR
 * example.  The examples are discussed in further detail in the paper:
 *
 *    Tuminaro, R.S., Shadid, J.N., and Hutchinson, H.A., Parallel Sparse
 *    Matrix Vector Multiply Software for Matrix with Data Locality",
 *    Sandia Technical Report, 1995.
 *
 *****************************************************************************/

extern int N_grid_pts, num_PDE_eqns;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"

/* external function declarations */

extern void create_msr_matrix(int update[], double **val, int **bindx, 
                              int N_update);
extern void add_row_3D(int row, int i, double val[], int bindx[], int *n);
extern void add_row_5pt(int row, int location, double val[], int bindx[], 
			int *n);
extern void add_row_9pt(int row, int location, int *next_nonzero, double val[],
                 int bindx[], int *n);
extern void create_vbr_matrix(int update_blks[], double **val, int **indx,
                       int num_update_blks, int **rnptr, int **bnptr,
                       int **bindx);
extern void check_memory_space(int bnptr_next, int total_blks,
                               int indx_next, int total_nz);
extern void add_vbr_row(int row, int location, double val[], int indx[],
                        int rnptr[], int bnptr[], int bindx[]);
extern int get_block_size(int );
extern void new_block(int *next , int row_blk_size, int col_blk_size,
                      int blk_col, double number, int bindx[], int indx[],
                      double val[]);

extern void read_triangles(int proc, int N_update);
extern void init_msr(double **val, int **bindx, int N_update);
extern void compress_matrix(double val[], int bindx[], int N_update);
extern void add_to_element(int row, int column, double element, double val[],
                           int bindx[], int diag);
extern void setup_Ke(double Ke[][3], double x1, double , double ,
                     double y2, double x3, double y3);
extern void read_coordinates(double **x, double **y, int update_indx[],int,int);
extern void add_to_element(int row, int column, double element, double val[],
                           int bindx[], int diag);
extern void create_fe_matrix(int update[], int proc, int **bindx, double **val,
                      int N_update);
extern void fill_fe_matrix(double val[], int bindx[], int update[], 
	int update_indx[], int external[], int extern_indx[], int data_org[],
	int proc_config[]);
extern void add_row_matrix_free(int row, int location, int *next_nonzero, 
        double val[], int bindx[]);
extern void matrix_vector_mult(int N_update, int update[], int update_indx[],
                      int N_external, int external[], int extern_indx[],
                      double x[], double y[], int data_org[],int proc_config[]);





/*****************************************************************************
 *                                                                           *
 *                             MSR                                           *
 *                                                                           *
 *****************************************************************************/

void create_msr_matrix(int update[], double **val, int **bindx, int N_update)

{
  int i,total_nz, n = -1;
  int avg_nonzeros_per_row = 9;


  total_nz = N_update*avg_nonzeros_per_row + 1;
  *bindx   = (int *) AZ_allocate(total_nz*sizeof(int));
  *val     = (double *) AZ_allocate(total_nz*sizeof(double));
  if ((*val == NULL) && (total_nz != 0) ) {
    (void) fprintf(stderr, "Error: Not enough space to create matrix\n");
    (void) fprintf(stderr,
                   "      Try reducing the variable 'avg_nonzeros_per_row'\n");
    exit(1);
  }
  for (i = 0 ; i < total_nz ; i++ ) (*bindx)[i] = 0;

  (*bindx)[0] = N_update+1;
  for (i = 0; i < N_update; i++) {
    add_row_3D(update[i], i, *val, *bindx,&n);
    if ( (*bindx)[i+1] > total_nz) {
      (void) fprintf(stderr, "Error:total_nz not large enough to accomodate");
      (void) fprintf(stderr, " nonzeros\n       Try increasing the variable");
      (void) fprintf(stderr, " 'avg_nonzeros_per_row'\n       Finished \n");
      (void) fprintf(stderr, "first %d rows out of %d\n", i, N_update);
      exit(1);
    }
  }

} /* create_msr_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_row_5pt(int row, int location, double val[], int bindx[], int *n)

/*
 * Add one row to an MSR matrix corresponding to a 5pt discrete approximation
 * to the 2D Poisson operator on an n x n square.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    val,bindx    == (see user's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{

/*
  static int n = -1;
*/
  int m;
  int        k;
  int        NP;

  /* determine grid dimensions */

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) sqrt( ((double) m)+0.01);
    if ( m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
    if (m == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }
  }
  *n = m;

  k  = bindx[location];
  NP = num_PDE_eqns;

  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */

  bindx[k] = row + NP;   if ((row/NP)%m !=     m-1) val[k++] = -1.00;
  bindx[k] = row - NP;   if ((row/NP)%m !=       0) val[k++] = -1.00;
  bindx[k] = row + m*NP; if ((row/(NP*m))%m != m-1) val[k++] = -1.00;
  bindx[k] = row - m*NP; if ((row/(NP*m))%m !=   0) val[k++] = -1.00;

  bindx[location+1] = k;
  val[location]     = 4.0;

} /* add_row_5pt */


/***************************************************************
 *                                                             *
 *           Other examples                                    *
 *                                                             *
 **************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_row_3D(int row, int location, double val[], int bindx[], int *n)

/*
 * Add one row to an MSR matrix corresponding to a 7pt discrete approximation
 * to the 3D Poisson operator on an n x n x n square.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    val,bindx    == (see Aztec User's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{

  int        m;
  int        k, NP;

  /* determine grid dimensions */

  m = *n; 
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) pow( ((double) m )+0.01, 0.33334);
    if (m == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }

    if (m*m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^3 where k is an integer\n");
      exit(1);
    }
  }
  *n = m;
  k = bindx[location];
  NP = num_PDE_eqns;

  /*
   * Check neighboring points in each direction and add nonzero entry if
   * neighbor exists.
   */

  bindx[k] = row + NP;     if ((row/NP)%m      != m-1) val[k++] = -1.0;
  bindx[k] = row - NP;     if ((row/NP)%m      != 0  ) val[k++] = -1.0;
  bindx[k] = row + m*NP;   if ((row/(NP*m))%m  != m-1) val[k++] = -1.0;
  bindx[k] = row - m*NP;   if ((row/(NP*m))%m  != 0  ) val[k++] = -1.0;
  bindx[k] = row + m*m*NP; if (bindx[k]    < m*m*m*NP) val[k++] = -1.0;
  bindx[k] = row - m*m*NP; if (bindx[k]          >= 0) val[k++] = -1.0;

  bindx[location+1] = k;
  val[location]     = 6.0;

} /* add_row_3D */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_row_9pt(int row, int location, int *next_nonzero, double val[],
                 int bindx[], int *n)

/*
 * Add one row to an MSR matrix corresponding to a 9pt discrete approximation
 * to the 2D Poisson operator on an n x n square.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    next_nonzero == points to the next free storage location in val[] and
 *                    bindx[].
 *                    On input, next_nonzero points to location where off
 *                    diagonals in new row will be stored.
 *                    On output, next_nonzero points to the next free location
 *                    (i.e. after the new row that was added).
 *    val,bindx    == (see file 'parameters'). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{
  int        i, j, point, old_ptr, stride, m;

  m = *n;
  if (m == -1) {
    m = N_grid_pts;
    m = (int ) sqrt(((double) m ) + 0.01);
    if (m <= 1) {
      (void) fprintf(stderr, "grid size too small\n");
      exit(1);
    }

    if (m*m != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
  }
  *n = m; 

  /* figure out where we are in the global grid */

  point = row/num_PDE_eqns;
  i     = point%m;
  j     = (point-i)/m;

  old_ptr = *next_nonzero;
  val[location] = 0.0;

  for (stride = 1; stride <= m; stride = stride*m) {

    /*
     * If we are on the bottom boundary, create entries for 2nd order
     * discretization. Otherwise create high order discretization.
     */

    if ((i != 0) && (i != m-1)) {
      if (i != 1) {
        val[*next_nonzero]   = 1.0/12.0;
        bindx[*next_nonzero] = row - stride*2*num_PDE_eqns;
        *next_nonzero        = *next_nonzero + 1;
      }
      val[*next_nonzero]   = -16.0/12.0;
      bindx[*next_nonzero] = row - stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
      if (i != m-2) {
        val[*next_nonzero]   =   1.0/12.0;
        bindx[*next_nonzero] = row + stride*2*num_PDE_eqns;
        *next_nonzero        = *next_nonzero + 1;
      }
      val[*next_nonzero]   = -16.0/12.0;
      bindx[*next_nonzero] = row + stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
    }
    else if (i != m-1) {
      val[*next_nonzero]   = -1.0;
      bindx[*next_nonzero] = row + stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
    }
    else if (i != 0) {
      val[*next_nonzero]   = -1.0;
      bindx[*next_nonzero] = row - stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
    }
    if ((i != 0) && (i != m-1)) val[location] += 30.0/12.0;
    else val[location] += 2.0;
    i = j;
  }
  bindx[location+1] = bindx[location] + (*next_nonzero - old_ptr);

} /* add_row_9pt */


/*****************************************************************************
 *                                                                           *
 *                             VBR                                           *
 *                                                                           *
 *****************************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void create_vbr_matrix(int update_blks[], double **val, int **indx,
                       int num_update_blks, int **rnptr, int **bnptr,
                       int **bindx)
/*
 * Create a matrix in VBR format (containing 'num_update_blks' block rows).
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/1/95
 *
 * Parameter list:
 *
 *  update_blks     == list of blocks updated by this processor
 *  val             == On output, contains matrix nonzeros
 *  indx            == On output, contains points to blocks of nonzeros
 *                     in val[]. IN particular, the pointers to the nonzero
 *                     blocks in block row i are
 *                        indx[bnptr[i], ... , bnptr[i+1]-1 ]
 *  num_update_blks == Number of blocks updated by this processor.
 *  rnptr           == On output, rnptr[0] = 0 and rnptr[i+1]-rnptr[i] gives
 *                     the row dimension of global block row 'update_blks[i]'.
 *  bnptr           == On output, bnptr[i] points to
 *                        1) the list of nonzero blocks in bindx[] of the
 *                           global row 'update_blks[i]'.
 *                        2) the list of pointers to nonzero matrix values
 *                           in indx[] of the global row 'update_blks[i]'.
 *  bindx           == On output, bindx[bnptr[i], ..., bnptr[i+1]-1] lists
 *                     the global block column indices of the nonzero
 *                     blocks in global row 'update_blks[i]'.
 */

{

  int avg_blks_per_row = 7;
  int avg_blk_size     = num_PDE_eqns*num_PDE_eqns;
  int i, total_nz, total_blks;

  /* --------------------- execution begins ----------------------------------*/

  total_blks = num_update_blks*avg_blks_per_row;
  total_nz   = total_blks*avg_blk_size + 1;

  *indx  = (int *)    AZ_allocate( (total_blks+1)*sizeof(int));
  *bindx = (int *)    AZ_allocate( (total_blks+1)*sizeof(int));
  *rnptr = (int *)    AZ_allocate( (num_update_blks+1)*sizeof(int));
  *bnptr = (int *)    AZ_allocate( (num_update_blks+1)*sizeof(int));
  *val   = (double *) AZ_allocate( (total_nz)*sizeof(double));

  if (*val == NULL) {
    (void) fprintf(stderr, "Error: Not enough space to create matrix\n"
                   "      Try reducing avg_blks_per_row,avg_blk_size\n");
    exit(1);
  }

  (*bnptr)[0] = 0;
  (*rnptr)[0] = 0;
  (*indx )[0] = 0;
  for (i = 0; i < num_update_blks; i++) {
    add_vbr_row(update_blks[i], i, *val, *indx, *rnptr, *bnptr, *bindx);

    /* do a series of memory checks to make sure we have not exceeded space */

    check_memory_space((*bnptr)[i+1], total_blks, (*indx)[(*bnptr)[i+1]],
                       total_nz);
  }

} /* create_vbr_matrix */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void add_vbr_row(int row, int location, double val[], int indx[], int rnptr[],
                 int bnptr[], int bindx[])

/*
 * Create one block row of a VBR matrix corresponding to block row 'row'
 * of a series of m uncoupled 2D Poisson equations.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/1/95
 *
 * Parameter list:
 *   row           == Global index of new block row to be added.
 *   location      == Place (local index) in rnptr,bnptr where the
 *                    new row is to be added.
 *   val           == On output, contains matrix nonzeros for new row.
 *   indx          == On output, contains new pointers to blocks of nonzeros
 *                    in val[] corresponding to new row. In particular, the
 *                    pointers to the nonzero blocks are in
 *                      indx[bnptr[location], ... , bnptr[location+1]-1]
 *   rnptr         == On output, rnptr[location+1] - rnptr[location] gives the
 *                    row dimension of the new block.
 *   bnptr         == On output, bnptr[location] points to
 *                       1) the global block column list of nonzero blocks in
 *                          bindx[] of the new block row.
 *                       2) the list of pointers to nonzero matrix values
 *                          in indx[] of the new block row.
 *   bindx         == On output, bindx[bnptr[location],...,bnptr[location+1]-1]
 *                    lists the global block column indices of the nonzero
 *                    blocks in the new row.
 */

{

  int    i, j, k, stride;
  static int n = -1;
  int        row_blk_size, col_blk_size;

  /* -------------------------- execution begins -----------------------------*/

  if (n == -1) {
    n = N_grid_pts;
    n = (int ) pow( ((double) n )+0.01,0.33334);
    if (n <= 1) {
      (void) fprintf(stderr, "grid size too small\n");
      exit(1);
    }

    if ( n*n*n != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^3 where k is an integer\n");
      exit(1);
    }
  }

  /* figure out where we are in the global grid */

  i = row%n;
  k = (row-i)/n;
  j = k%n;
  k = (k-j)/n;

 /* put in the diagonal block */

  row_blk_size = get_block_size(row);
  bnptr[location+1] = bnptr[location];
  rnptr[location+1] = rnptr[location] + row_blk_size;
  new_block(&bnptr[location+1], row_blk_size, row_blk_size, row, 6.0,
            bindx, indx, val);

  /* check in each of the directions and add a new block column */

  for (stride = 1; stride <= n*n; stride *= n) {
    if (i != 0) {
      col_blk_size = get_block_size(row-stride);
      new_block(&bnptr[location+1], row_blk_size, col_blk_size, row - stride,
                -1.0, bindx, indx, val);
    }

    if (i != n-1) {
      col_blk_size = get_block_size(row+stride);
      new_block(&bnptr[location+1], row_blk_size, col_blk_size, row + stride,
                -1.0, bindx, indx, val);
    }
    if (stride == 1) i = j;
    else i = k;
  }

} /* add_vbr_row */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void check_memory_space(int indx_next, int total_blks, int val_next,
                        int total_nz)

/*
 * Check that we have not used more space then we allocated in creating
 * the vbr matrix.
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/1/95
 *
 * Parameters:
 *
 *    indx_next        == pointer to next available location in indx[].
 *    total_blks       == total number of blocks for which we have allocated
 *                        space.
 *    val_next         == point to next available location in val[].
 *    total_nz         == total number of nonzeros allocated for matrix.
 */

{

  if (indx_next > total_blks) {
    (void) fprintf(stderr, "Error: total_blks not large enough to accomodate "
                   "nonzeros\n       Try increasing the variable "
                   "'avg_blks_per_row'\n");
    exit(1);
  }

  if (val_next > total_nz) {
    (void) fprintf(stderr,
                   "Error: total_nz not large enough to accomodate nonzeros\n"
                   "      Try increasing the variable 'avg_blk_size'\n");
    exit(1);
  }

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int get_block_size(int location)

/*
 * Return the row block size of the global block row 'location'.
 */

{
  if (location < 0)  {
     (void) fprintf(stderr,"Error: Improper grid location (%d) ",location);
     (void) fprintf(stderr,"inside get_block_size()\n");
     exit(-1);
  }

  return num_PDE_eqns;

}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void new_block(int *next, int row_blk_size, int col_blk_size, int blk_col,
               double number, int bindx[], int indx[], double val[])

/*
 * Add a new block to the VBR matrix, initialize all the offdiagonals to zeros,
 * and initialize all the diagonal elements to 'number'.
 *
 * Author: Ray Tuminaro
 * Date:   3/15/95
 * Parameters:
 *   next          == On input, indx[*next] points to the new block.
 *                    On output, *next points to the next free location
 *                    in indx[] and bindx[].
 *   row_blk_size  == number of rows in the new block to be added.
 *   col_blk_size  == number of columns in the new block to be added.
 *   blk_col       == Block column index of new block to be added.
 *   number        == On output, each element of the 'diagonal' of the new
 *                    block contains 'number'
 *   bindx,indx,val== See file 'parameters'. On output, bindx[],indx[],val[]
 *                    have been appended such that the new block is added
 *                    to the VBR matrix.
 */

{

  int k, kk, start;

  bindx[*next]  = blk_col;
  indx[*next+1] = indx[*next] + row_blk_size*col_blk_size;
  start = indx[*next];
  for (k = 0; k < row_blk_size; k++) {
    for (kk = 0; kk < col_blk_size; kk++) {
      if (k == kk) val[start + kk*row_blk_size+k] = number;
      else val[start + kk*row_blk_size+k] = 0.0;
    }
  }

  (*next)++;

}


/*****************************************************************************
 *                                                                           *
 *                             FINITE ELEMENT MSR                            *
 *                                                                           *
 *****************************************************************************/

/******************************************************************************/
/******************************************************************************/

#define NOT_FOUND -1
FILE *fp;
int  *T[3];
int  N_triangle;

void create_fe_matrix(int update[], int proc, int **bindx, double **val,
                      int N_update)

/*
 * Create a matrix in MSR format corresponding to a finite element
 * discretization of the 2D Poisson operator. The finite element grid
 * is contained in the files 'fe_grid_#' (where each file contains only
 * the grid information for processor #).
 *
 * NOTE: this routine does not actually compute the nonzero entries
 * of the matrix. It only computes the nonzero pattern (i.e. bindx[]).
 * The subroutine fill_fe_matrix() computes the nonzero entries and
 * must be called after the communication pattern has been set.
 *
 *
 *      Author:         Lydie Prevost Div 1422 SNL
 *      Date:           11/25/1994
 *
 * Parameter List:
 *
 *  update  == list of unknowns updated by this processor
 *  bindx   == On output, contains nonzero locations in matrix where
 *
 *                 bindx[i],  0 <=i< N_update: ptr to offdiagonal nonzeros
 *                                             for row 'update[i]'.
 *                 bindx[i],  N_update > i     contains column numbers
 *                                             corresponding to nonzero val[i].
 *  N_update== Number of unknowns updated by this processor.
 *  proc    == Node number of this processor.
 */

{

  int i, j,k;
  int row;

  /**************************** execution begins ******************************/

  /* read data */

  read_triangles(proc, N_update);

  /* initialize msr matrix */

  init_msr(val, bindx, N_update);

  for (k = 0; k < N_triangle; k++) {
    for (i = 0; i < 3; i++) {
      row = AZ_find_index(T[i][k], update, N_update);
      if (row != NOT_FOUND ) {
        for (j = 0; j < 3; j++)
          add_to_element(row, T[j][k], 0.0, *val, *bindx, i==j);
      }
    }
  }

  compress_matrix(*val, *bindx, N_update);

} /* create_fe_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void fill_fe_matrix(double val[], int bindx[], int update[], int update_indx[],
                    int external[], int extern_indx[], int data_org[],
		    int proc_config[])

/*
 * Fill in the nonzero entries of an MSR matrix corresponding to a finite
 * element discretization of the 2D Poisson operator. It is assumed that
 * the routine create_fe_matrix() has already been called.
 *
 *      Author:         Lydie Prevost Div 1422 SNL
 *      Date:           11/25/1994
 *
 * Parameter List:
 *
 *  update   == list of unknowns updated by this processor
 *  val      == On output, contains matrix nonzeros where
 *
 *               val[i],0 <= i < N_update: diagonal for row update[i].
 *               val[bindx[i]-->val[bindx[i+1]-1]: nonzeros for row 'update[i]'.
 *
 *  N_update == Number of unknowns updated by this processor.
 *  proc     == Node number of this processor.
 *
 *      Author:         Lydie Prevost Div 1422 SNL
 *      Date:           11/25/1994
 *
 */

{

  /* local variables */

  int    i_triangle, i, j;
  double Ke[3][3];
  double *x,*y;
  int    row;
  int    N_update,N_external;

  /**************************** execution begins ******************************/

  N_update   = data_org[AZ_N_internal] + data_org[AZ_N_border];
  N_external = data_org[AZ_N_external];

  read_coordinates(&x, &y, update_indx, N_update, N_external);

  AZ_exchange_bdry(x, data_org, proc_config);
  AZ_exchange_bdry(y, data_org, proc_config);

  /* relabel the triangle with the appropriate local index number */

  for (i_triangle = 0; i_triangle < N_triangle; i_triangle++) {
    for (i = 0; i < 3; i++) {
      row = AZ_find_index(T[i][i_triangle], update, N_update);

      if (row == NOT_FOUND) {
        row = AZ_find_index(T[i][i_triangle], external, N_external);
        T[i][i_triangle] = extern_indx[row];
      }
      else T[i][i_triangle] = update_indx[row];
    }
  }

  /* fill in the discretization */

  for (i_triangle = 0; i_triangle < N_triangle; i_triangle++) {

    /* initialize submatrix Ke containing contributions for this triangle */

    setup_Ke(Ke, x[T[0][i_triangle]], y[T[0][i_triangle]], x[T[1][i_triangle]],
             y[T[1][i_triangle]], x[T[2][i_triangle]], y[T[2][i_triangle]]);

    /* store submatrix Ke elements in a[] */

    for (i = 0; i < 3; i++) {
      if (T[i][i_triangle] < N_update) {
        for (j = 0; j < 3; j++) {
          add_to_element(T[i][i_triangle], T[j][i_triangle], Ke[i][j], val,
                         bindx, i==j);
        }
      }
    }
  }

} /* fill_fe_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_to_element(int row, int column, double element, double val[],
                    int bindx[], int diag)

/*
 * Add the value 'element' to the matrix entry A(row,column) and store the
 * result back in the matrix.
 *
 * Note: A is encoded in the MSR arrays val[] and bindx[].
 */

{

  int i, start, end;

  /* (row,column) corresponds to a diagonal element of the matrix */

  if (diag) val[row] += element;
  else {
    start = bindx[row];
    end   = bindx[row+1] - 1;

    /* search for the column */

    for (i = start; i <= end; i++) {
      if ( (bindx[i] == -1) || (bindx[i] == column)) break;
    }

    /* this entry did not previously exist in the matrix */

    if (bindx[i] == -1) {
      val[i]   = element;
      bindx[i] = column;
    }

    else if (bindx[i] == column) {
      val[i] += element;
    }

    else {
      (void) fprintf(stderr, "Error:  Not enough room for element (%d,%d) in "
                     "matrix\n", row, column);
      (void) fprintf(stderr, "      Increase MAX_NZ_ROW\n");
      for (i = start; i <= end; i++)
        (void) fprintf(stderr, "%4d ", bindx[i]);
      exit(1);
    }
  }

} /* add_to_element */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void compress_matrix(double val[], int bindx[], int N_points)

/*
 * Take an existing MSR matrix which contains some columns labeled '-1'
 * and compress the matrix by eliminating all the columns labeled '-1'.
 * Note: this routine is usually used in conjunction with init_msr().
 */

{

  int free_ptr,start,end,i,j;
  int temp;

  free_ptr = bindx[0];
  start    = bindx[0];

  for (i = 0 ; i < N_points ; i++ ) {
    bindx[i] = bindx[i+1] - bindx[i];
  }
  for (i = 0 ; i < N_points ; i++ ) {
    end   = start + bindx[i] -1;

    for (j = start ; j <= end ; j++ ) {
      if (bindx[j] != -1) {
        val[free_ptr]   = val[j];
        bindx[free_ptr] = bindx[j];
        free_ptr++;
      }
      else {
        bindx[i]--;
      }
    }
    start = end+1;
  }

  start = N_points+1;
  for (i = 0 ; i <= N_points ; i++ ) {
    temp = bindx[i];
    bindx[i] = start;
    start += temp;
  }

} /* compress_matrix */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void read_triangles(int proc, int N_update)

/*
 * Open the finite element grid files and read the triangle information into the
 * array T[][].
 */

{

  char string[40];
  int  i, *ttt;

  sprintf(string, "fe_grid_%d", proc);
  fp = (FILE *) fopen(string, "r");
  if (fp == NULL) {
    (void) fprintf(stderr, "Error: finite element grid file (%s) not found\n",
                    string);
    exit(-1);
  }

  fscanf(fp, "%d", &N_triangle);
  T[0] = (int *) AZ_allocate(N_triangle*sizeof(int));
  T[1] = (int *) AZ_allocate(N_triangle*sizeof(int));
  T[2] = (int *) AZ_allocate(N_triangle*sizeof(int));
  if ( (N_triangle != 0) && ( T[2] == NULL) ) {
     printf("Not enough space to read the triangles\n");
     exit(1);
  }
  ttt = T[0];
  for (i = 0 ; i < N_triangle ; i++ ) ttt[i] = 0;
  ttt = T[1];
  for (i = 0 ; i < N_triangle ; i++ ) ttt[i] = 0;
  ttt = T[2];
  for (i = 0 ; i < N_triangle ; i++ ) ttt[i] = 0;

  for (i = 0; i < N_triangle; i++)
    fscanf(fp, "%d%d%d", &(T[0][i]), &(T[1][i]), &(T[2][i]));

  /*
   * Now read in the number of points on this processor and make sure it matches
   * the number of points specified in main() by the user.
   */

  fscanf(fp, "%d", &i);
  if (i != N_update) {
    (void) fprintf(stderr, "Error: The number of points assigned to this\n"
                   "       processor does not match the number of\n"
                   "       points in finite element file: %d vs. %d\n",
                   N_update, i);
    exit(1);
  }

} /* read_triangles */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

#define MAX_NZ_ROW 20

void init_msr(double **val, int **bindx, int N_points)

/*
 * Initialize an MSR matrix containing N_points rows and MAX_NZ_ROW
 * nonzeros per row. The nonzero values of the matrix are set to
 * 0 and the columns are set to -1.
 */

{
  int i;

  /* allocate space for sparse matrix */

  *bindx = (int    *) AZ_allocate((N_points*MAX_NZ_ROW+1)*sizeof(int));
  *val   = (double *) AZ_allocate((N_points*MAX_NZ_ROW+1)*sizeof(double));
  if ( *val == NULL) {
     printf("Not enough room for msr matrix\n");
     exit(1);
  }

  for (i = 0; i < N_points*MAX_NZ_ROW+1; i++) {
    (*val)[i]   = 0.0;
    (*bindx)[i] = -1;
  }

  (*bindx)[0] = N_points+1;
    for (i = 0; i < N_points; i++) {
      (*bindx)[i+1] = (*bindx)[i] + MAX_NZ_ROW - 1;
    }

} /* init_msr */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void setup_Ke(double Ke[][3], double xa, double ya, double xb, double yb,
              double xc, double yc)

/*
 * Compute finite element contributions corresponding to 1 triangle
 * and store them in the 3 x 3 matrix Ke.
 *
 * Parameters:
 *    (xa,ya),(xb,yb),(xc,yc)  == Coordinates of triangle vertices.
 *    Ke                       == On output, contains fe contributions
 *                                corresponding to above triangle.
 */

{

  int    ii, jj;
  double det_J;

  Ke[0][0] = (yc-yb)*(yc-yb) + (xc-xb)*(xc-xb);
  Ke[0][1] = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
  Ke[0][2] = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);

  Ke[1][0] = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
  Ke[1][1] = (yc-ya)*(yc-ya) + (xc-xa)*(xc-xa);
  Ke[1][2] = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
  Ke[2][0] = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);
  Ke[2][1] = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
  Ke[2][2] = (yb-ya)*(yb-ya) + (xb-xa)*(xb-xa);

  det_J = (xb-xa)*(yc-ya)-(xc-xa)*(yb-ya);
  det_J = 2*det_J;

  for (ii = 0; ii < 3; ii++) {
    for (jj = 0; jj < 3; jj++) {
      Ke[ii][jj] = Ke[ii][jj] / det_J;
    }
  }

} /* setup_Ke */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void read_coordinates(double **x, double **y, int update_indx[], int N_update,
                      int N_external)

/*
 * Read in coordinates from previously openned file pointed to by 'fp'.
 */

{
  int i;
  int total;

  total = N_update + N_external;

  *x   = (double *) AZ_allocate((total+1)*sizeof(double));
  *y   = (double *) AZ_allocate((total+1)*sizeof(double));
  if ( *y == NULL) {
     printf("Out of space in read_coordinates\n");
     exit(1);
  }

   for (i = 0 ; i < total+1 ; i++ ) (*x)[i] = 0.0;
   for (i = 0 ; i < total+1 ; i++ ) (*y)[i] = 0.0;

  for (i = 0; i < N_update; i++) {
    fscanf(fp,"%lf%lf", &((*x)[update_indx[i]]), &((*y)[update_indx[i]]));
  }

  fclose(fp);

}


/*****************************************************************************
 *                                                                           *
 *                             MATRIX FREE                                   *
 *                                                                           *
 *****************************************************************************/


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void add_row_matrix_free(int row, int location, int *next_nonzero, double val[],
                         int bindx[])

/*
 * Add the nonzero pointers for one new row to an MSR matrix corresponding to
 * a 5pt discrete approximation to the 2D Poisson operator on an n x n square.
 * NOTE: the actual nonzero values of the matrix are not computed and stored
 * in the routine (only the pointers: bindx).
 *
 * Author: Ray Tuminaro, Div 1422
 * Date:   3/15/95
 *
 * Parameters:
 *    row          == the global row number of the new row to be added.
 *    location     == the local row number where the diagonal of the new row
 *                    will be stored.
 *    next_nonzero == points to the next free storage location in bindx[].
 *                    On input, next_nonzero points to location where off
 *                    diagonals in new row will be stored.
 *                    On output, next_nonzero points to the next free location
 *                    (i.e. after the new row that was added).
 *    bindx        == (see Aztec User's guide). On output, bindx[] is
 *                    appended such that the new row has been added.
 */

{
  int        i, j, point, old_ptr, stride;
  static int n = -1;

  /* determine grid dimensions */

  if (n == -1) {
    n = N_grid_pts;
    n = (int ) sqrt(((double) n )+0.01);

    if (n == 1) {

      /* special case */

      val[0] = 1.0;
      bindx[1] = bindx[0];
      return;
    }

    if (n*n != N_grid_pts) {
      (void) fprintf(stderr, "Error: the total number of points (%d) needs\n",
                     N_grid_pts);
      (void) fprintf(stderr, "       to equal k^2 where k is an integer\n");
      exit(1);
    }
  }

  /* figure out where we are in the global grid */

  point   = row/num_PDE_eqns;
  i       = point%n;
  j       = (point-i)/n;
  old_ptr = *next_nonzero;

  /*
   * Check for neighboring points in each direction and add nonzero entry if the
   * neighbor exits.
   */

  for (stride = 1; stride <= n; stride = stride*n) {
    if (i != 0 ) {
      bindx[*next_nonzero] = row - stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
    }

    if (i != n-1) {
      bindx[*next_nonzero] = row + stride*num_PDE_eqns;
      *next_nonzero        = *next_nonzero + 1;
    }
    i = j;
  }

  bindx[location+1] = bindx[location] + (*next_nonzero - old_ptr);

} /* add_row_matrix_free */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void matrix_vector_mult(int N_update, int update[], int update_indx[],
                      int N_external, int external[], int extern_indx[],
                      double x[], double y[], int data_org[], int proc_config[])

/*
 * Example of a matrix-free matrix-vector multiply corresponding to a 5-pt
 * 2D Poisson operator on a rectangular grid.
 *
 * Author: Ray Tuminaro
 * date: 3/15/96
 *
 * Parameters:
 *    N_update  == Number of points updated on this processor.
 *    update    == List (global indecies) of points updated on this node.
 *    update_indx   == update_index[i] gives the local numbering of global
 *                     point 'update[i]'.
 *    N_external  == Number of external points on this processor.
 *    external    == List (global indecies) of external points on this node.
 *    extern_indx   == extern_index[i] gives the local numbering of global
 *                     point 'external[i]'.
 *    x             == Input vector to be multiplied with matrix.
 *    y             == On output, y = A x (where A is discrete approximation
 *                     to the Poisson operator.
 */

{

  int        i, j, ii;
  int        point, tnew;
  int        stride;
  static int n = -1;

  if (n == -1) {
    n = N_grid_pts;
    n = (int ) sqrt(((double) n )+0.01);
  }

  AZ_exchange_bdry(x, data_org, proc_config);

  for (ii = 0; ii < N_update; ii = ii + 1) {

    /* compute the location */

    point = update[ii];
    i     = point%n;
    j     = (point-i)/n;

    y[ii] = 4.0*x[ii];

    for (stride = 1; stride <= n; stride = stride*n) {
      if (i != 0) {
        tnew = AZ_find_index(point-stride, update, N_update);

        if (tnew == NOT_FOUND) {
          tnew = AZ_find_index(point-stride, external, N_external);
          tnew = extern_indx[tnew];
        }
        else
          tnew = update_indx[tnew];

        y[ii] = y[ii] - x[tnew];
      }

      if (i != n-1) {
        tnew = AZ_find_index(point+stride, update, N_update);
        if (tnew == NOT_FOUND) {
          tnew = AZ_find_index(point+stride, external, N_external);
          tnew = extern_indx[tnew];
        }
        else
          tnew = update_indx[tnew];

        y[ii] -= x[tnew];
      }

      i = j;
    }
  }

} /* matrix_vector_mult */
