/*@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef AZTEC_MPI
#include <mpi.h>
#endif
#include "az_aztec.h"
#include "azk_komplex.h"

int main(int argc, char *argv[])
{
  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */

  int    *bindx_real;              /* index and values arrays for MSR matrices */
  double *val_real, *val_imag;

  int * update;                    /* List of global eqs owned by the processor */
  double *x_real, *b_real;         /* initial guess/solution, RHS  */
  double *x_imag, *b_imag;

  int    N_local;                  /* Number of equations on this node */
  double residual;                 /* Used for computing residual */

  double *xx_real, *xx_imag, *xx; /* Known exact solution */
  int myPID, nprocs;

  AZ_MATRIX *Amat_real;             /* Real matrix structure */
  AZ_MATRIX  *Amat;                 /* Komplex matrix to be solved. */
  AZ_PRECOND *Prec;                 /* Komplex preconditioner */
  double *x, *b;                    /* Komplex Initial guess and RHS */

  int i;

  /******************************/
  /* First executable statement */
  /******************************/

#ifdef AZTEC_MPI
  MPI_Init(&argc,&argv);
#endif

  /* Get number of processors and the name of this processor */

#ifdef AZTEC_MPI
  AZ_set_proc_config(proc_config,MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config,0);
#endif

  nprocs = proc_config[AZ_N_procs];
  myPID  = proc_config[AZ_node];

  printf("proc %d of %d is alive\n",myPID, nprocs);

  /* Define two real diagonal matrices. Will use as real and imaginary parts */ 

  /* Get the number of local equations from the command line */
  if (argc!=2)
   {
    if (myPID==0) printf("Usage: %s number_of_local_equations\n",argv[0]);
    exit(1);
   }
  N_local = atoi(argv[1]);

 /* Need N_local+1 elements for val/bindx arrays */
  val_real = malloc((N_local+1)*sizeof(double));
  val_imag = malloc((N_local+1)*sizeof(double));

  /* bindx_imag is not needed since real/imag have same pattern  */
  bindx_real = malloc((N_local+1)*sizeof(int));

  update = malloc(N_local*sizeof(int)); /* Malloc equation update list */

  b_real = malloc(N_local*sizeof(double)); /* Malloc x and b arrays */
  b_imag = malloc(N_local*sizeof(double));
  x_real = malloc(N_local*sizeof(double));
  x_imag = malloc(N_local*sizeof(double));
  xx_real = malloc(N_local*sizeof(double));
  xx_imag = malloc(N_local*sizeof(double));

  for (i=0; i<N_local; i++) 
    {
      val_real[i] = 10 + i/(N_local/10); /* Some very fake diagonals */
      val_imag[i] = 10 - i/(N_local/10); /* Should take exactly 20 GMRES steps */

      x_real[i] = 0.0;         /* Zero initial guess */
      x_imag[i] = 0.0; 

      xx_real[i] = 1.0;        /* Let exact solution = 1 */
      xx_imag[i] = 0.0;

      /* Generate RHS to match exact solution */
      b_real[i] = val_real[i]*xx_real[i] - val_imag[i]*xx_imag[i];
      b_imag[i] = val_imag[i]*xx_real[i] + val_real[i]*xx_imag[i];

      /* All bindx[i] have same value since no off-diag terms */
      bindx_real[i] = N_local + 1;

      /* each processor owns equations 
	 myPID*N_local through myPID*N_local + N_local - 1 */
      update[i] = myPID*N_local + i; 
      
    }

  bindx_real[N_local] = N_local+1; /* Need this last index */

  /* Register Aztec Matrix for Real Part, only imaginary values are needed*/

  Amat_real = AZ_matrix_create(N_local);

  AZ_set_MSR(Amat_real, bindx_real, val_real, NULL, N_local, update, AZ_GLOBAL);

  /* initialize AZTEC options */
 
  AZ_defaults(options, params);
  options[AZ_solver]  = AZ_gmres; /* Use CG with no preconditioning */
  options[AZ_precond] = AZ_none;
  options[AZ_kspace] = 21;
  options[AZ_max_iter] = 21;
  params[AZ_tol] = 1.e-14;


    /**************************************************************/
    /* Construct linear system.  Form depends on input parameters */
    /**************************************************************/

       /**************************************************************/
       /* Method 1:  Construct A, x, and b in one call.              */
       /* Useful if using A,x,b only one time. Equivalent to Method 2*/
       /**************************************************************/
      
     AZK_create_linsys_ri2k (x_real,  x_imag,  b_real,  b_imag,
                    options,  params, proc_config,
                    Amat_real, val_imag, &x, &b, &Amat);
      
       /**************************************************************/
       /* Method 2:  Construct A, x, and b in separate calls.        */
       /* Useful for having more control over the construction.      */
       /* Note that the matrix must be constructed first.            */
       /**************************************************************/

     /* Uncomment these three calls and comment out the above call

     AZK_create_matrix_ri2k (options,  params, proc_config,
			     Amat_real, val_imag, &Amat);

      AZK_create_vector_ri2k(options,  params, proc_config, Amat,
			     x_real, x_imag, &x);

      AZK_create_vector_ri2k(options,  params, proc_config, Amat,
			     b_real, b_imag, &b);
     */

    /**************************************************************/
    /* Build exact solution vector.                               */
    /* Check residual of init guess and exact solution            */
    /**************************************************************/

      AZK_create_vector_ri2k(options,  params, proc_config, Amat,
			     xx_real, xx_imag, &xx);

      residual = AZK_residual_norm(x, b, options, params, proc_config, Amat);
    if (proc_config[AZ_node]==0)
      printf("\n\n\nNorm of residual using initial guess = %12.4g\n",residual);
 
      residual = AZK_residual_norm(xx, b, options, params, proc_config, Amat);
      AZK_destroy_vector(options,  params, proc_config, Amat, &xx);
    if (proc_config[AZ_node]==0)
      printf("\n\n\nNorm of residual using exact solution = %12.4g\n",residual);

    /**************************************************************/
    /* Create preconditioner                                      */
    /**************************************************************/

   AZK_create_precon(options,  params, proc_config, x, b, Amat, &Prec);

    /**************************************************************/
    /* Solve linear system using Aztec.                           */
    /**************************************************************/

   AZ_iterate(x, b, options, params, status, proc_config, Amat, Prec, NULL);

    /**************************************************************/
    /* Extract solution.                                          */
    /**************************************************************/

     AZK_extract_solution_k2ri(options, params, proc_config, Amat, Prec, x, 
			       x_real,  x_imag); 
    /**************************************************************/
    /* Destroy Preconditioner.                                    */
    /**************************************************************/

      AZK_destroy_precon (options,  params, proc_config, Amat, &Prec);

    /**************************************************************/
    /* Destroy linear system.                                     */
    /**************************************************************/

      AZK_destroy_linsys (options,  params, proc_config, &x, &b, &Amat);
 
  if (proc_config[AZ_node]==0)
    {
      printf("True residual norm squared   = %22.16g\n",status[AZ_r]);
      printf("True scaled res norm squared = %22.16g\n",status[AZ_scaled_r]);
      printf("Computed res norm squared    = %22.16g\n",status[AZ_rec_r]);
    }

  /* Print comparison between known exact and computed solution */
  {double sum = 0.0;
  
  for (i=0; i<N_local; i++) sum += fabs(x_real[i]-xx_real[i]);
  for (i=0; i<N_local; i++) sum += fabs(x_imag[i]-xx_imag[i]);
  printf("Processor %d:  Difference between exact and computed solution = %12.4g\n",
	 proc_config[AZ_node],sum);
  }
  /*  Free memory allocated */

  free((void *) val_real );
  free((void *) bindx_real );
  free((void *) val_imag );
  free((void *) update );
  free((void *) b_real );
  free((void *) b_imag );
  free((void *) x_real );
  free((void *) x_imag );
  free((void *) xx_real );
  free((void *) xx_imag );
  AZ_matrix_destroy(&Amat_real);

#ifdef AZTEC_MPI
  MPI_Finalize();
#endif

return 0 ;
}
