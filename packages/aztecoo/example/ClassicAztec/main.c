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

/*******************************************************************************
 * Sample driver for AZTEC package.  The software is tested by setting up a
 * system of equations and a right hand side and then solving the system of
 * equations using AZTECs iterative solvers.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"

int N_grid_pts;                  /* number of grid points in application      */

int num_PDE_eqns;                /* number of equations (usually PDEs)        */
                                 /* associated with each grid point.          */

int application;                 /* problem that is being run:                */
                                 /*  0: Finite Difference Poisson (MSR format)*/
                                 /*  1: Finite Difference Poisson (VBR format)*/
                                 /*  2: Finite Element Poisson (MSR format)   */
                                 /*  3: Read matrix from file using           */
                                 /*     AZ_read_msr_matrix() (MSR format)     */

/* -------------  external function declarations ------------------------- */

extern void init_matrix_vector_structures(int proc_config[],
                                          int *update_index[], int *update[],
                                          int  *data_org[], int *external[],
                                          int *extern_index[], int
                                          input_option, double *a[],
                                          int *bindx[], int *indx[],
                                          int *bpntr[], int *rpntr[],
                                          int *cpntr[]);

extern void init_options(int options[], double params[]);

extern void init_guess_and_rhs(int update_index[], int update[], double *x[],
                               double *ax[], int data_org[], double val[],
                               int indx[], int bindx[], int rpntr[],
                               int cpntr[], int bpntr[], int proc_config[]);

extern void fill_fe_matrix(double val[],int bindx[], int update[],
                           int update_index[], int external[],
                           int extern_index[], int data_org[]);

/* ------------------------------------------------------------------------*/

int main(int argc, char *argv[])

/* Set up and solve a test problem defined in the subroutine
   init_matrix_vector_structures().

   Author:   Ray Tuminaro, Div 1422, Sandia National Labs
   date:     11/10/94

 ******************************************************************************/

{

  double *ax,*x;                   /* ax is the right hand side for the test
                                      problem. x is the approximate solution
                                      obtained using AZTEC.                  */

  int    i,input_option;


  /* See Aztec User's Guide for more information   */
  /* on the variables that follow.                 */

  int    proc_config[AZ_PROC_SIZE];/* Processor information:                 */
  /*  proc_config[AZ_node] = node name      */
  /*  proc_config[AZ_N_procs] = # of nodes  */

  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.   */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.        */
  int    *data_org;                /* Array to specify data layout */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve()
                                      indicating success or failure.         */

  int    *update,                  /* vector elements (global index) updated
                                      on this processor.                     */
    *external;                     /* vector elements needed by this node.   */

  int    *update_index;            /* ordering of update[] and external[]    */
  int    *extern_index;            /* locally on this processor. For example
                                      update_index[i] gives the index
                                      location of the vector element which
                                      has the global index 'update[i]'.      */

                                   /* Sparse matrix to be solved is stored
                                      in these arrays.                       */
  int    *rpntr,*cpntr,*indx, *bpntr, *bindx;
  double *val;


  /* ----------------------- execution begins --------------------------------*/

  /* Put the # of processors, the node id,    */
  /* and an MPI communicator into proc_config */

#ifdef AZTEC_MPI
  MPI_Init(&argc,&argv);
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config, AZ_NOT_MPI );
#endif

  /*
   * Read and broadcast: problem choice, problem size, equations per grid point
   * and how we wish to initialize 'update'.
   */

  if (proc_config[AZ_node] == 0) {
    (void) printf("enter the application problem number\n");
    (void) printf("  = 0: Finite Difference MSR Poisson on n x n x n grid.\n");
    (void) printf("  = 1: Finite Difference VBR Poisson on n x n x n grid.\n");
    (void) printf("  = 2: Finite Element MSR Poisson\n");
    (void) printf("  = 3: Use AZ_read_msr_matrix() to read file '.data'\n");
    scanf("%d",&application);
    if ((application < 0) || (application > 3)){
      (void) fprintf(stderr, "Error: Invalid application (%d) selected\n",
                     application);
      exit(1);
    }

    if (application == 0) {
      (void) printf("\nNote: To try other problems, change add_row_3D()");
      (void) printf("\n      in create_msr_matrix() to add_row_5pt() or");
      (void) printf("\n      add_row_9pt().\n\n");
    }
    if (application == 2) {
      (void) printf("\nNote: Input files are provided for 1 finite element ");
      (void) printf("\n      problem. This problem can be run on either 1  ");
      (void) printf("\n      or 4 processors. To run on 1 processor, copy  ");
      (void) printf("\n      the file fe_1proc_grid_0 to fe_grid_0. To run on");
      (void) printf("\n      4 processors, copy the files fe_4proc_grid_k to ");
      (void) printf("\n      fe_grid_k (k = 0,1,2,3). In both cases enter 197");
      (void) printf("\n      when prompted for the number of grid points and ");
      (void) printf("\n      linear when prompted for the partitioning!!!\n\n");
    }

    if (application == 3)
      (void) printf("enter the total number of matrix rows\n");
    else (void) printf("enter the total number of grid points\n");
    scanf("%d", &N_grid_pts);

    num_PDE_eqns = 1;
    if (application < 2) {
      (void) printf("enter the number of equations per grid point\n");
      scanf("%d", &num_PDE_eqns);
    }

    (void) printf("partition option \n");
    (void) printf("     = %d: linear\n", AZ_linear);
    (void) printf("     = %d: update pts from file '.update'\n", AZ_file);
    if (application < 2)
      (void) printf("     = %d: box decomposition\n", AZ_box);
    scanf("%d", &input_option);
  }
  AZ_broadcast((char *) &N_grid_pts  , sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &num_PDE_eqns, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &input_option, sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) &application , sizeof(int), proc_config, AZ_PACK);
  AZ_broadcast((char *) NULL         , 0          , proc_config, AZ_SEND);

  /* create an application matrix for AZTEC */

  init_matrix_vector_structures(proc_config, &update_index, &update, &data_org,
                                &external, &extern_index, input_option, &val,
                                &bindx, &indx, &bpntr, &rpntr, &cpntr);

  /* initialize AZTEC options */

  init_options(options,params);

  if ( (i = AZ_check_input(data_org, options, params, proc_config) ) < 0) {
    AZ_print_error(i);
    exit(-1);
  }

  /* Matrix fill for finite element example (see Aztec User's Guide). */

  if (application == 2)
    fill_fe_matrix(val, bindx, update, update_index, external, extern_index,
                   data_org);

  /* Initialize right hand side and initial guess */
  /* NOTE: STORAGE ALLOCATED FOR 'x' IS GREATER THAN THE NUMBER */
  /*       OF MATRIX ROWS PER PROCESSOR. 'x' INCLUDES SPACE FOR */
  /*       EXTERNAL (GHOST) ELEMENTS. THUS, THE SIZE OF 'x' IS  */
  /*       'data_org[AZ_N_internal] + data_org[AZ_N_border] +   */
  /*       data_org[AZ_N_external]'.                            */

  init_guess_and_rhs(update_index, update, &x, &ax, data_org, val, indx, bindx,
                     rpntr, cpntr, bpntr, proc_config);

  /* update[], update_index[], external[], extern_index[] are used to map
   * between Aztec's ordering of the equations and the user's ordering
   * (see the User's guide for more details). If these mapping arrays 
   * are not needed by the user, they can be deallocated as they are not 
   * used by AZ_solve().
   */

  free((void *) update);   free((void *) update_index);
  free((void *) external); free((void *) extern_index);

  /* solve the system of equations using ax as the right hand side */

  AZ_solve(x,ax, options, params, indx, bindx, rpntr, cpntr, bpntr, val,
           data_org, status, proc_config);

  /* Free allocated memory */

  free((void *) x);        free((void *) ax);           free((void *) indx);
  free((void *) bindx);    free((void *) rpntr);        free((void *) cpntr);
  free((void *) bpntr);    free((void *) val);          free((void *) data_org);

#ifdef AZTEC_MPI
  MPI_Finalize();
#endif
  return(1);

}
