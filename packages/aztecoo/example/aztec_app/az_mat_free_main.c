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
static char rcsid[] = "$Id$"
;
#endif

/*******************************************************************************
 * Sample matrix-free driver for AZTEC package.  The software is tested by 
 * solving a simple system (Poisson equation) via Aztec's iterative solvers.
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "example_specific.h"

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int main(int argc, char *argv[])

/*
 * Set up and solve a linear system of equations using Aztec's 
 * matrix-free interface. The matrix system corresponds to 
 * Poisson's equation on a regular 2-dimensional mesh of size
 * n x n. The associated matrix has n^2 equations (one for each 
 * mesh grid point). The kth equation corresponds to grid point 
 * (i,j) where    
 *            k = n*j + i       and  i = 0, ... , n-1
 *                                   j = 0, ... , n-1
 * and the equation is of the form
 *
 *    4*u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1) = f(i,j)
 *
 * where u(i,j) corresponds to the unknown located at the (i,j)th 
 * grid point and f(i,j) is the right-hand side vector corresponding
 * to the (i,j)th grid point.
 * 
 * For example a 4x4 grid would correspond to a 16 equation matrix
 * with the following grid point numbering:
 *
 *                     12 13 14 15 
 *                      8  9 10 11 
 *                      4  5  6  7 
 *                      0  1  2  3 
 * 
 * The equations are solved using p^2 processors where the processors
 * are distributed as a 2-D array. For example with 4 processors, the
 * above equations would be distributed as follows:
 *      processor 0 would contain equations  0, 1, 4, 5
 *      processor 1 would contain equations  2, 3, 6, 7
 *      processor 2 would contain equations  8, 9,12,13
 *      processor 3 would contain equations 10,11,14,15
 */

{
/*****************************************************************************/
/* Generic Aztec variables common to most matrix-free applications           */
/*****************************************************************************/
                     /*                                                      */
int                  /*                                                      */
   proc_config[      /* Processor information.                               */
     AZ_PROC_SIZE],  /*                                                      */
   options[          /* Used to pass parameters to Aztec.                    */
     AZ_OPTIONS_SIZE];/*                                                     */
double               /*                                                      */
   params[           /* Used to pass in parameters to Aztec                  */
     AZ_PARAMS_SIZE],/*                                                      */
   status[           /* Returned from AZ_iterate() indicating success or     */
     AZ_STATUS_SIZE];/* failure of the iterative solvers                     */ 
AZ_MATRIX *Amat=NULL,/* Structure representing matrix to be solved.          */
          *Pmat=NULL;/* Structure representing preconditioning matrix.       */
                     /*                                                      */
AZ_PRECOND *Prec;    /* Structure representing preconditioner.               */
                     /*                                                      */
double *x, *rhs;     /* Initial guess and right-hand side to linear system.  */

/*****************************************************************************/
/* Example specific variables                                                */
/*****************************************************************************/
int N_equations;     /* Number of equations in the linear system residing on */
                     /* this processor. In our example, this is 'nx*ny'.     */
                     /*                                                      */
int N_coord;         /* Total number of points in each coordinate direction. */
                     /* That is, the matrix dimension is of size N_coord^2   */
                     /*                                                      */
int  Nproc_coord;    /* Number of processor in each coordinate direction     */
                     /* That is, the linear system is solved by a processor  */
                     /* array of size Nproc_coord x Nproc_coord.             */
                     /*                                                      */
int nx, ny;          /* Number of grid points in each coordinate direction   */
                     /* within the local grid owned by this processor.       */
                     /*                                                      */
struct pass_data     /* Used to pass data through Aztec into the USER's      */
       pass_data;    /* matrix-vector product routine.                       */
struct exchange      /* Used to exchange data within the USER's matrix-vector*/
      comm_structure;/* product routine.                                     */

int i, N_ghost, precond_choice;
/* ----------------------- execution begins --------------------------------*/

  /* get number of processors and the name of this processor */

#ifdef AZTEC_MPI
  MPI_Init(&argc,&argv);
   AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
#else
   AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  /* initialize AZTEC options */

   AZ_defaults(options, params);


   precond_choice = example_specific_read_input(options, proc_config,
                                                &N_coord,&Nproc_coord);

  /* Set up the example specific data structure 'comm_structure' which  */
  /* handles communication. In addition, set (nx, ny) to the number of  */
  /* grid points in the x and y direction for this processor.           */
  
   example_specific_comm_setup(N_coord, &nx, &ny, proc_config[AZ_node], 
                               Nproc_coord, &comm_structure);
   N_equations = nx*ny;

  /* Set up an EXAMPLE SPECIFIC data structure to pass information      */
  /* through Aztec to the USER's routines. NOTE: while specific to      */
  /* this example, most codes will do something similar.                */

   pass_data.nx    = nx;                      
   pass_data.ny    = ny;                      
   pass_data.nproc = (int) sqrt( ((double) proc_config[AZ_N_procs]) + .01);
   pass_data.py    = proc_config[AZ_node]/pass_data.nproc;
   pass_data.px    = proc_config[AZ_node]%pass_data.nproc;
   pass_data.comm_structure  = &comm_structure;

  /* initialize Aztec matrix to be solved */

   Amat = AZ_matrix_create(N_equations);
   AZ_set_MATFREE(Amat, &pass_data, example_specific_matvec);

   Prec = NULL;    /* When preconditioning structure is NULL, AZ_iterate()  */
                   /* applies Aztec's preconditioners to application matrix */
                   /* (i.e. user does not supply a preconditioning routine  */
                   /* or an additional matrix for preconditioning.          */

   switch (options[AZ_precond]) {

   
   case AZ_none:           /***** NO preconditioning. ***********************/
     break;

   case AZ_ls:             /***** Polynomial preconditioning (least-squares */
   case AZ_Neumann:        /***** or Neumann). ******************************/
     AZ_set_MATFREE_matrix_norm(Amat, 8.0); /* Aztec needs upper bound for  */
     break;                                 /* matrix norm                  */

   case AZ_Jacobi:         /***** Jacobi preconditioning. *******************/
      /*                                                                    */
      /* Matrix diagonal is needed and can be supplied via:                 */
      /*     1) a getrow() function associated with Amat                    */
      /*     2) a new matrix (MSR, VBR, or USER with getrow function)       */
      /*        associated with preconditioner.                             */
      /* We choose the second option in this example.                       */


      Pmat = AZ_matrix_create(N_equations);
      AZ_set_MATFREE_getrow(Pmat, NULL, example_specific_diagonal_getrow, 
                            NULL, 0, proc_config);
      Prec = AZ_precond_create(Pmat, AZ_precondition, NULL);
      break;

   case AZ_dom_decomp:     /***** Domain Decomposition preconditioning. *****/
      /*                                                                    */
      /* Local matrix must be supplied via a getrow associated with Amat or */
      /* via a new matrix associated with preconditioner.                   */

      if (precond_choice == NONOVERLAPDD_PRECOND) {

         /* Since no overlapping is done, communication information and     */
         /* matrix data corresponding to ghost unknowns is not needed.      */

         Pmat = AZ_matrix_create(N_equations);
         AZ_set_MATFREE_getrow(Pmat, &pass_data,simple_example_specific_getrow, 
			       NULL, 0, proc_config);
         Prec = AZ_precond_create(Pmat, AZ_precondition, NULL);
         options[AZ_overlap] = AZ_none;
         options[AZ_subdomain_solve]  = AZ_lu; /* use LU solver on  */
                                               /* each subdomain.   */
      }
      else {

         /* Since overlapping is done, communication information and matrix */
         /* data corresponding to ghost unknowns is supplied.               */

         N_ghost = 0;
         if (pass_data.px !=    0   ) N_ghost += ny;
         if (pass_data.py !=    0   ) N_ghost += nx;
         if (pass_data.px != pass_data.nproc - 1) N_ghost += ny;
         if (pass_data.py != pass_data.nproc - 1) N_ghost += nx;

         AZ_set_MATFREE_getrow(Amat, &pass_data, example_specific_getrow, 
			     example_specific_comm_wrapper,N_ghost,proc_config);
         options[AZ_subdomain_solve]  = AZ_lu; /* use LU solver on  */
                                               /* each subdomain.   */
         options[AZ_overlap] = 2;
      }
      break;

   case AZ_user_precond:   /***** User supplied preconditioner.         *****/
      /*                                                                    */
      /* A local Gauss-Seidel preconditioner.                               */
      /*                                                                    */
     Prec = AZ_precond_create(Amat, example_specific_precond, &pass_data);
     AZ_set_precond_print_string(Prec,"My Gauss Seidel");
     break;
  }

   /* Set initialize guess to all zeros and rhs to a delta       */
   /* function (ie. zero everywhere except at one point).        */

   x   = (double *) malloc(N_equations*sizeof(double));
   rhs = (double *) malloc(N_equations*sizeof(double));
   if ( rhs == NULL) { 
      printf("Not enough space for rhs\n"); exit(1); 
   }
   for (i = 0 ; i < N_equations; i++) x[i]   = 0.0;
   for (i = 0 ; i < N_equations; i++) rhs[i] = 0.0;
   if (proc_config[AZ_node] == 0) rhs[0] = 1.; 

   /* solve linear system using Aztec. */

   AZ_iterate(x, rhs, options, params, status,proc_config,Amat,Prec,NULL);

   /* Free allocated memory */

   free (x);      free (rhs);
   example_specific_frees(&comm_structure);

   AZ_matrix_destroy(&Amat);
   if (Pmat != NULL) AZ_matrix_destroy(&Pmat);
   if (Prec != NULL) AZ_precond_destroy(&Prec);

#ifdef AZTEC_MPI
  MPI_Finalize();
#endif
  return(1);

} /* main */


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

int example_specific_read_input(int *options, int *proc_config,int *N_coord,
				 int *Nproc_coord){
/*
 * Read the user's input and fill the array 'options' and the variables
 * N_coord (number of grid points in each coordinate direction) and 
 * Nproc_coord (number of nodes in processor array in each coordinate
 * direction).
 */
  int N_grid_pts, precond_choice;

  /* Read and broadcast the size of the linear system that will    */
  /* be solved. Since this linear system will correspond to the    */
  /* discretization of a PDE on a square grid (N_coord x N_coord), */
  /* the sqrt() of the linear system size must be an integer.      */
  /* Also read the preconditioning choice.                         */

   if (proc_config[AZ_node] == 0) {
     printf("enter the total number of grid points.\n");
     printf("Must be a square(i.e. x^2).\n");
     scanf("%d", &N_grid_pts);
     printf("Enter the type of preconditioning:\n");
     printf("    %1d:   no preconditioning\n",NO_PRECOND);
     printf("    %1d:   polynomial/least-squares\n",LS_PRECOND); 
     printf("    %1d:   polynomial/Neumman series\n",NS_PRECOND); 
     printf("    %1d:   3-step Jacobi\n",JAC_PRECOND); 
     printf("    %1d:   Nonoverlapping Domain decomposition\n",
			NONOVERLAPDD_PRECOND); 
     printf("    %1d:   Overlapping Domain decomposition\n",OVERLAPDD_PRECOND); 
     printf("    %1d:   User supplied\n",USER_PRECOND); 
     scanf("%d", &precond_choice);
   }
   AZ_broadcast((char *) &N_grid_pts    ,sizeof(int), proc_config, AZ_PACK);
   AZ_broadcast((char *) &precond_choice,sizeof(int), proc_config, AZ_PACK);
   AZ_broadcast((char *) NULL           ,0          , proc_config, AZ_SEND);

   if      (precond_choice == NO_PRECOND) options[AZ_precond] = AZ_none;
   else if (precond_choice == LS_PRECOND) options[AZ_precond] = AZ_ls;
   else if (precond_choice == NS_PRECOND) options[AZ_precond] = AZ_Neumann;
   else if (precond_choice == JAC_PRECOND) options[AZ_precond] = AZ_Jacobi;
   else if (precond_choice == NONOVERLAPDD_PRECOND) 
	options[AZ_precond] = AZ_dom_decomp;
   else if (precond_choice == OVERLAPDD_PRECOND) 
	options[AZ_precond] = AZ_dom_decomp;
   else if (precond_choice == USER_PRECOND) options[AZ_precond]=AZ_user_precond;
   else {
      printf("Invalid preconditoner\n");
      exit(1);
   }


   /* figure out the number of grid points in each coordinate     */
   /* direction. Also figure out the number of processors in      */
   /* each coordinate direct. That is, the processors are         */
   /* configured as a processor array (Nproc_coord x Nproc_coord) */
   /* and so the sqrt() of the total number of processors must    */
   /* also be an integer.                                         */

   *N_coord     = (int) sqrt(((double) N_grid_pts) + .01);
   *Nproc_coord = (int) sqrt(((double) proc_config[AZ_N_procs]) + .01); 

   if ((*N_coord)*(*N_coord) != N_grid_pts) {
      if (proc_config[AZ_node] == 0)
         printf("Number of grid points must be a square (i.e. n^2).\n");
      exit(1);
   }
   if ((*Nproc_coord)*(*Nproc_coord) != proc_config[AZ_N_procs]) {
      if (proc_config[AZ_node] == 0)
         printf("Number of processors must be a square (i.e. p^2).\n");
      exit(1);
   }
   return(precond_choice);
}

/******************************************************************************/
/******************************************************************************/

void example_specific_precond(double x[], int input_options[],
        int proc_config[], double input_params[], AZ_MATRIX *Amat,
        AZ_PRECOND *prec)


/******************************************************************************/
/*
 * A LOCAL Gauss-Seidel routine subroutine.
 *
 * Parameters:
 * =========
 * x                         On input, a vector. On output, x[] is
 *                           smoothed by taking averages with neighbors.
 *
 * prec                      On input, AZ_get_precond_data(prec) points to
 *                           that data_structure 'pass_data' which contains
 *                           the local grid size (nx,ny) on this processor.
 *
 * Amat, input_options,      Not used.
 * proc_config, Amat,
 * input_params
 *
 */


{
   int i, j, nx, ny, me;

   struct pass_data *pass_data;     /* Data passing structure. This user-     */
                                    /* defined data structure is used to pass */
                                    /* information through Aztec and back into*/
                                    /* the user's subroutines.                */
   double *rhs, temp;

  /*-------------------------------------------------------------------------*/

    pass_data      = (struct pass_data *) AZ_get_precond_data(prec);
    nx = pass_data->nx;
    ny = pass_data->ny;
    rhs            = (double *) malloc(nx*ny*sizeof(double));
    if (rhs == NULL) {
       printf("Out of memory in preconditioning routine\n");
       exit(1);
    }
    for (i = 0; i < nx*ny; i++) rhs[i] = x[i];
    for (i = 0; i < nx*ny; i++) x[i]   = 0.0;

    for (i = 0 ; i < nx; i++ ) {
       for (j = 0 ; j < ny; j++ ) {
          me = j*nx + i;
          temp = rhs[me];
          if (i !=   0 ) temp += x[me-1];
          if (i != nx-1) temp += x[me+1];
          if (j !=   0 ) temp += x[me-nx];
          if (j != ny-1) temp += x[me+nx];
          x[me] = temp/4.;
       }
    }
    free(rhs);
}
