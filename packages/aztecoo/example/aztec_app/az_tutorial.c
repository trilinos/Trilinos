/*@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }

int n =6;            /* POISSON EQUATION IS SOLVED ON an n x n GRID. This    */
                     /* corresponds to an x-vector of length n-squared and   */ 
                     /* thus the matrix A is of size n-squared by n-squared. */
#define MAX_NZ_ROW 5 /*  Max number of nonzero elements in any matrix row    */

extern void create_matrix_row(int row,int i,double val[],int bindx[]);

int main(int argc, char *argv[])

/* Set up a Poisson test problem and solve it with AZTEC.                 */
{

  double *b,*x;                    /* rhs and approximate solution          */
  int    i, nrow;

  /* See Aztec User's Guide for the variables that follow:         */

  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
  int    *data_org;                /* Array to specify data layout          */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
  int    *update,                  /* vector elements updated on this node. */
         *external;                /* vector elements needed by this node.  */
  int    *update_index;            /* ordering of update[] and external[]   */
  int    *extern_index;            /* locally on this processor.            */
  int    *bindx;                   /* Sparse matrix to be solved is stored  */
  double *val;                     /* in these MSR arrays.                  */
  int    N_update;                 /* # of unknowns updated on this node    */


  /* get number of processors and the name of this processor */
#ifdef AZTEC_MPI
  MPI_Init(&argc,&argv);
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif



  /* Define partitioning:  matrix rows (ascending order) owned by this node */

  nrow = n*n;
  AZ_read_update(&N_update, &update, proc_config, nrow, 1, AZ_linear);

  /*
   * Create the matrix: each processor creates only rows appearing in update[]
   * (using global col. numbers).
   */

  bindx = (int    *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(int));
  val   = (double *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(double));
  if (val == NULL) perror("Error: Not enough space to create matrix");

  bindx[0] = N_update+1;

  for (i = 0; i < N_update; i++) {
    create_matrix_row(update[i], i, val, bindx);
  }

  /* convert matrix to a local distributed matrix */

  AZ_transform(proc_config, &external, bindx, val, update, &update_index,
               &extern_index, &data_org, N_update, NULL, NULL, NULL, NULL,
               AZ_MSR_MATRIX);

  /* initialize AZTEC options */

  AZ_defaults(options, params);

  /* Set rhs (delta function at lower left corner) and initialize guess */

  b = (double *) malloc(N_update*sizeof(double));
  x = (double *) malloc((N_update + data_org[AZ_N_external])*sizeof(double));
          /* NOTE: SOLUTION VECTOR MUST CONTAIN SPACE FOR EXTERNAL ELEMENTS */
  if ((x == NULL) && (i != 0)) perror("Not enough space in rhs");
  for (i = 0; i < N_update; i++) {
    x[update_index[i]] = 0.0;
    b[update_index[i]] = 0.0;
    if (update[i] == 0) b[update_index[i]] = 1.0;
  }

  /* update[], update_index[], external[], extern_index[] are used to map
   * between Aztec's ordering of the equations and the user's ordering
   * (see the User's guide for more details). If these mapping arrays
   * are not needed by the user, they can be deallocated as they are not
   * used by AZ_solve().
   */
 
  free((void *) update);   free((void *) update_index);
  free((void *) external); free((void *) extern_index);


  /* solve the system of equations using b  as the right hand side */

  AZ_solve(x, b, options, params, NULL, bindx, NULL, NULL, NULL, val, data_org,
           status, proc_config);

  /* Free allocated memory */
 
  free((void *) x);    free((void *) b);       free((void *) bindx);    
  free((void *) val);  free((void *) data_org);
 

#ifdef AZTEC_MPI
  MPI_Finalize();
#endif
  return(1);

} /* main */

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

void create_matrix_row(int row, int location, double val[], int bindx[])

/* Add one row to an MSR matrix corresponding to a discrete approximation to the
 * 2D Poisson operator on an n x n square.
 *
 * Parameters:
 *    row          == global row number of the new row to be added.
 *    location     == local row where diagonal of the new row will be stored.
 *    val,bindx    == (see user's guide). On output, val[] and bindx[]
 *                    are appended such that the new row has been added.
 */

{
  int k;

  /* check neighbors in each direction and add nonzero if neighbor exits */

  k = bindx[location];
  bindx[k]  = row + 1; if ((row  )%n != n-1) val[k++] = -1.;
  bindx[k]  = row - 1; if ((row  )%n !=   0) val[k++] = -1.;
  bindx[k]  = row + n; if ((row/n)%n != n-1) val[k++] = -1.;
  bindx[k]  = row - n; if ((row/n)%n !=   0) val[k++] = -1.;

  bindx[location+1] = k;  val[location]     = 4.; /* matrix diagonal */
}
