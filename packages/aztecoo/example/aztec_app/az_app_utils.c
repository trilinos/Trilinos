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

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"

/* ----------------- external function declarations ------------------- */

extern void create_msr_matrix(int *, double **, int **, int );
extern void AZ_transform(int *, int **, int *, double *, int *, int **,int **,
                         int **, int, int *, int *, int *, int **, int );
extern void create_fe_matrix(int *, int , int **, double **, int );
extern void create_vbr_matrix(int *, double **, int **, int, 
                              int **, int **, int **);
extern void init_options(int *, double *);
extern void init_guess_and_rhs(int *, int *, double **,double **,int *, 
                               double *, int *, int *, int *, int *, 
                               int *, int *);
extern void init_matrix_vector_structures(int *, int **, int **, int  **, 
                               int **, int **, int , double **, int **, 
                               int **, int **, int **, int **);
extern void print_global_element(int ,int *,int *,int *, int *, double *,int *);

extern int num_PDE_eqns;   /* number of PDEs being solved         */
extern int application;    /* problem being solved                */
extern int N_grid_pts;     /* N_grid_pts is the total number of   */
                           /* grid points used in the simulation. */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_options(int options[], double params[])

{

  /*
   * Choose among AZTEC options (see User's Guide).
   */

  AZ_defaults(options, params);

  options[AZ_solver]   = AZ_cgs;
  options[AZ_scaling]  = AZ_none;
  options[AZ_precond]  = AZ_ls;
  options[AZ_conv]     = AZ_r0;
  options[AZ_output]   = 1;
  options[AZ_pre_calc] = AZ_calc;
  options[AZ_max_iter] = 1550;
  options[AZ_poly_ord] = 5;
  options[AZ_overlap]  = AZ_none;
  options[AZ_kspace]   = 60;
  options[AZ_aux_vec]  = AZ_resid;

  params[AZ_tol]       = 4.00e-9;
  params[AZ_drop]      = 0.0;
  params[AZ_ilut_fill] = 1.5;
  params[AZ_omega]     = 1.;

} /* init_options */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void init_guess_and_rhs(int update_index[], int update[], double *x[],double
                        *ax[],int data_org[], double val[], int indx[], int
                        bindx[], int rpntr[], int cpntr[], int bpntr[], int
                        proc_config[])

/*
 * Set the initial guess and the right hand side where the right hand side
 * is obtained by doing a matrix-vector multiplication.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date :  3/15/95
 *
 * Parameters
 *
 *    update_index   ==      On input, ordering of update and external
 *                           locally on this processor. For example
 *                           'update_index[i]' gives the index location
 *                           of the block which has the global index
 *                           'update[i]'.
 *    update         ==      On input, list of pts to be updated on this node
 *    data_org       ==      On input, indicates how data is set on this node.
 *                           For example, data_org[] contains information on
 *                           how many unknowns are internal, external and
 *                           border unknowns as well as which points need
 *                           to be communicated. See User's Guide for more
 *                           details.
 *    val, indx,     ==      On input, holds matrix nonzeros. See User's Guide
 *    bindx, rpntr,          for more details.
 *    cpntr, bpntr
 *    x              ==      On output, 'x' is allocated and set to all zeros.
 *    ax             ==      On output, 'ax' is allocated and is set to the
 *                           result of a matrix-vector product.
 */

{

  int i,j;
  int temp,num;
  double sum = 0.0;
  AZ_MATRIX *Amat;

  temp = data_org[AZ_N_int_blk]  + data_org[AZ_N_bord_blk];
  num  = data_org[AZ_N_internal] + data_org[AZ_N_border];

  /* allocate vectors */

  i       = num + data_org[AZ_N_external];
  *x      = (double *) AZ_allocate((i+1)*sizeof(double));
  *ax     = (double *) AZ_allocate((i+1)*sizeof(double));
  if (*ax == NULL) {
    (void) fprintf(stderr, "Not enough space in init_guess_and_rhs() for ax\n");
    exit(1);
  }
  for (j = 0 ; j < i ; j++ ) (*x)[j] = 0.0;
  for (j = 0 ; j < i ; j++ ) (*ax)[j] = 0.0;

  /* initialize 'x' to a function which will be used in matrix-vector product */

  if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
    for (i = 0; i < temp; i++) {
      for (j = rpntr[i]; j < rpntr[i+1]; j++) {
        (*x)[j] = (double) (update[i]) + (double)(j-rpntr[i]) /
          (double)(num_PDE_eqns);
      }
    }
  }
  else {
    for (i = 0; i < temp; i++) {
      (*x)[i] = (double) (update[i]) / (double) (num_PDE_eqns);
    }
  }

  /* Reorder 'x' so that it conforms to the transformed matrix */
 
  AZ_reorder_vec(*x,data_org,update_index,rpntr);

  if (application == 2) {

    /* take out the constant vector. Used for the */
    /* finite element problem because it is singular */

    sum = AZ_gsum_double(sum, proc_config);
    i   = AZ_gsum_int(num, proc_config);
    if (i != 0) sum = sum / ((double) i);
    for (i = 0; i < num; i++) (*x)[i] -= sum;
  }
  Amat = AZ_matrix_create(num);

  if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX)
     AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);
  else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
     AZ_set_VBR(Amat, rpntr,cpntr, bpntr, indx,bindx, val, data_org, 0, NULL,AZ_LOCAL);


  Amat->matvec(*x, *ax, Amat, proc_config);
  AZ_matrix_destroy( &Amat );

  for (i = 0; i < num; i++) (*x)[i] = 0.0;

} /* init_guess_and_rhs */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_matrix_vector_structures(int proc_config[], int *update_index[], int
                                   *update[], int  *data_org[], int *external[],
                                   int *extern_index[], int input_option, double
                                   *val[], int *bindx[], int *indx[], int
                                   *bpntr[], int *rpntr[], int *cpntr[])


/*
 * Read in the points to be updated on this processor, create the global
 * distributed form of the application matrix, and then convert it to a
 * local distributed form for AZTEC kernels. Along the way, initialize the
 * following quantities:
 *     update_index[], update[], data_org[], a[], bindx[], bpntr[], cpntr[],
 *     rpntr[], indx[], external[], extern_index[].
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date:   3/15/95
 *
 * Parameters
 *
 *    proc_config    ==      On input, processor information:
 *                              proc_config[AZ_node] = name of this processor
 *                              proc_config[AZ_N_procs] = # of processors used
 *    update         ==      On output, list of pts to be updated on this node
 *    val,bindx      ==      On output, local distributed form of arrays
 *                           holding matrix values
 *    external       ==      On output, list of external vector elements
 *    update_index   ==      On output, ordering of update and external
 *    extern_index   ==      locally on this processor. For example
 *                           'update_index[i]' gives the index location
 *                           of the block which has the global index
 *                           'update[i]'.
 *    data_org       ==      On output, indicates how the data is set out on
 *                           this node. For example, data_org[] contains
 *                           information on how many unknowns are internal,
 *                           external, and border unknowns as well as which
 *                           points need to be communicated. See User's Guide
 *                           for more details.
 *    input_option   ==      Indicates how update[] will be initialized.
 *                           = 0, linear decomposition
 *                           = 1, points read from file 'update'.
 *                           = 2, box decomposition
 *                           See AZ_read_update() comments for more details.
 *
 *      The default finite difference MSR problem corresponds to a setting up
 *      a series of uncoupled 3D Poisson equations on a cube.
 *      To solve other problems, the call 'add_row_3D(...)' in
 *      'create_msr_matrix()' can be changed to 'add_row_5pt()' or
 *      'add_row_9pt()'.
 */

{

  int    N_update;            /* Number of pts updated on this processor     */
  int    MSRorVBR;
  int    chunks;
int blk_size, num_blk_cols,num_blk_rows,size,kk, convert_to_vbr = 0;
double *val2;
int    *bindx2;


  MSRorVBR = AZ_MSR_MATRIX;
  if (application == 1) MSRorVBR = AZ_VBR_MATRIX;

  chunks = num_PDE_eqns;
  if (MSRorVBR == AZ_VBR_MATRIX) chunks = 1;

  /* initialize the list of global indices. NOTE: the list of global */
  /* indices must be in ascending order so that subsequent calls to  */
  /* AZ_find_index() will function properly. */

  AZ_read_update(&N_update, update, proc_config, N_grid_pts, chunks,
                 input_option);

  /* create the matrix: each processor creates only the      */
  /* rows appearing in update[] ... however this row is      */
  /* created as if it were on a serial machine (i.e. using   */
  /* the global column numbers)                              */

  if (application == 1)
    create_vbr_matrix(*update, val, indx, N_update, rpntr, bpntr, bindx);
  else {
    *indx = NULL; *bpntr = NULL; *rpntr = NULL; *cpntr = NULL;

    if (application == 0) create_msr_matrix(*update, val, bindx, N_update);
    if (application == 2) create_fe_matrix(*update, proc_config[AZ_node],
                                           bindx, val, N_update);
    if (application == 3) { 
        AZ_read_msr_matrix(*update, val, bindx, N_update, proc_config);
    }
  }

  /* convert matrix to a distributed parallel matrix */

  AZ_transform(proc_config, external, *bindx, *val, *update, update_index,
               extern_index, data_org, N_update, *indx, *bpntr, *rpntr, cpntr,
               MSRorVBR);

  if ( (convert_to_vbr == 1) && (application == 3) ) {
     if (proc_config[AZ_node] == 0 ) {
	 printf("enter the block size\n");
	 scanf("%d",&blk_size);
     }
     AZ_broadcast((char *) &blk_size,  sizeof(int), proc_config, AZ_PACK);
     AZ_broadcast((char *) NULL         , 0          , proc_config, AZ_SEND);

     if ( N_update%blk_size != 0 ) {
        (void) fprintf(stderr," The block size must be a multiple of the number of rows per processor.\n");
        exit(-1);
     }

     num_blk_rows = N_update/blk_size;
     num_blk_cols = ( (*data_org)[AZ_N_external] + N_update)/blk_size;
     *cpntr = (int *) AZ_allocate( (num_blk_cols+2)*sizeof(int));
     *rpntr = (int *) AZ_allocate( (num_blk_cols+2)*sizeof(int));
     *bpntr = (int *) AZ_allocate( (num_blk_cols+2)*sizeof(int));
     size   = 20*(num_blk_cols+2);
     *indx  =  (int *) AZ_allocate(size*sizeof(int));
     bindx2 = *bindx;
     val2   = *val;
     *bindx = (int *) AZ_allocate(size*sizeof(int));
     *val   =  (double *) AZ_allocate(size*blk_size*blk_size*sizeof(double));

     for (kk = 0 ; kk < num_blk_cols ; kk++ ) (*cpntr)[kk] = blk_size;
     AZ_msr2vbr(*val,*indx,*rpntr,*cpntr,*bpntr,*bindx,bindx2,val2,
		num_blk_rows,num_blk_cols,size,size*blk_size*blk_size,blk_size);
     MSRorVBR = AZ_VBR_MATRIX;
     N_update /= blk_size;
     num_PDE_eqns = blk_size; 
     for (kk = 0 ; kk < N_update ; kk++ )
           (*update)[kk] = (*update)[blk_size*kk]/blk_size;
     for (kk = 0 ; kk < (*data_org)[AZ_N_external] ; kk++ ) 
           (*external)[kk] = (*external)[blk_size*kk]/blk_size;

     (*data_org)[AZ_matrix_type] = AZ_VBR_MATRIX;
     (*data_org)[AZ_N_int_blk ] /= blk_size;
     (*data_org)[AZ_N_bord_blk] /= blk_size;
     (*data_org)[AZ_N_ext_blk ] /= blk_size;
     AZ_free(bindx2);  AZ_free(val2);
  }


} /* init_matrix_vector_structures */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void print_global_element(int element,int update[],int data_org[],
           int update_index[], int rpntr[], double vector[],int proc_config[])
{
/*
 * Print out the vector element corresponding to the global numbering
 * 'element'. Note: if the VBR format is used, this routine will print
 * out all the vector elements corresponding to this block.
 *
 * Author: Ray Tuminaro, Div 1422, SNL
 * Date:   6/15/96
 *
 * Parameters
 *
 *    element        ==      On input, global number of vector element that
 *                           will be printed.
 *    update         ==      On input, list of pts updated on this node
 *    data_org       ==      On input, indicates how the data is set out on
 *                           this node. For example, data_org[] contains
 *                           information on how many unknowns are internal,
 *                           external, and border unknowns as well as which
 *                           points need to be communicated. See User's Guide
 *                           for more details.
 *    update_index   ==      On input, ordering of update locally on this
 *                           processor. For example, 'update_index[i]' gives 
 *                           the index location of the block which has the 
 *                           global index 'update[i]'.
 *    rpntr          ==      On input, rpntr[i+1]-rpntr[i] gives the block
 *                           size of the ith local block.
 *    vector         ==      On input, vector to be printed (just one element).
 *    proc_config    ==      On input, processor information:
 *                              proc_config[AZ_node] = name of this processor
 *                              proc_config[AZ_N_procs] = # of processors used
 */
   int i,k;

   /* synchronize things */

   fflush(stdout);
   i = AZ_gsum_int(1,proc_config);
 


   i = AZ_find_index(element,update,
                     data_org[AZ_N_int_blk]+data_org[AZ_N_bord_blk]);
   if (i !=-1) {
      i = update_index[i];
      if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
         fprintf(stdout,"(%d) = %e\n",element,vector[i]);
      else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX) {
        for (k = rpntr[i]; k < rpntr[i+1]; k++ ) 
           fprintf(stdout,"(%d,%d) = %e\n",element,k-rpntr[i],vector[k]);
      }
      fflush(stdout);
   }

   /* synchronize things */
   i = AZ_gsum_int(i,proc_config);

}

