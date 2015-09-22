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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#define AZTEC_MPI
#define AZ_MPI
#endif
#ifdef AZTEC_MPI
#include "mpi.h"
#endif
#include "az_aztec.h"
#include "prototypes.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define double_quote '"'

#ifdef __cplusplus
extern "C" void AZOO_iterate(double * xsolve, double * b, 
			     int * options, double * params, 
			     double * status, int *proc_config,
			     AZ_MATRIX * Amat);
#endif
int main(int argc, char *argv[])
{
  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
  int    *data_org;
                                   /* Array to specify data layout          */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
  int    *update;                  /* vector elements updated on this node. */
  int    *external;
                                   /* vector elements needed by this node.  */
  int    *update_index;
                                   /* ordering of update[] and external[]   */
  int    *extern_index;
                                   /* locally on this processor.            */
  int    *indx;   /* MSR format of real and imag parts */
  int    *bindx;
  int    *bpntr;
  int    *rpntr;
  int    *cpntr;
  AZ_MATRIX *Amat;
  double *val;
  double *x, *b, *xexact, *xsolve;
  int    n_nonzeros, n_blk_nonzeros;
  int    N_update;           /* # of block unknowns updated on this node    */
  int    N_local;
                                 /* Number scalar equations on this node */
  int    N_global, N_blk_global; /* Total number of equations */
  int    N_external, N_blk_eqns;

  double *val_msr;
  int *bindx_msr;
  
  double norm, d ;

  int matrix_type;

  int has_global_indices, option;
  int i, j, m, mp ;
  int ione = 1;

#ifdef AZTEC_MPI
  MPI_Init(&argc,&argv);
#endif

  /* get number of processors and the name of this processor */
 
#ifdef AZTEC_MPI
  AZ_set_proc_config(proc_config,MPI_COMM_WORLD);
#else
  AZ_set_proc_config(proc_config,0);
#endif

  printf("proc %d of %d is alive\n",
	 proc_config[AZ_node],proc_config[AZ_N_procs]) ;

#ifdef VBR
  if(argc != 3) 
    perror("error: enter name of data and partition file on command line") ; 
#else
  if(argc != 2) perror("error: enter name of data file on command line") ; 
#endif
  /* Set exact solution to NULL */
  xexact = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

#ifdef VBR
  read_hb(argv[1], proc_config, &N_global, &n_nonzeros, 
	  &val_msr,  &bindx_msr, &x, &b, &xexact);
  
  create_vbr(argv[2], proc_config, &N_global, &N_blk_global,
	     &n_nonzeros, &n_blk_nonzeros, &N_update, &update,
	     bindx_msr, val_msr, &val, &indx, 
	     &rpntr, &cpntr, &bpntr, &bindx);

  if(proc_config[AZ_node] == 0) 
    {
      free ((void *) val_msr);
      free ((void *) bindx_msr);
      free ((void *) cpntr);
    }
    matrix_type = AZ_VBR_MATRIX;


  distrib_vbr_matrix( proc_config, N_global, N_blk_global, 
		      &n_nonzeros, &n_blk_nonzeros,
		      &N_update, &update, 
		      &val, &indx, &rpntr, &cpntr, &bpntr, &bindx, 
		      &x, &b, &xexact);

#else
    read_hb(argv[1], proc_config, &N_global, &n_nonzeros,
             &val,  &bindx, &x, &b, &xexact);

  distrib_msr_matrix(proc_config, N_global, &n_nonzeros, &N_update,
		  &update, &val, &bindx, &x, &b, &xexact);

    matrix_type = AZ_MSR_MATRIX;
#endif
  /* convert matrix to a local distributed matrix */
  AZ_transform(proc_config, &external, bindx, val, update,
	       &update_index, &extern_index, &data_org, 
	       N_update, indx, bpntr, rpntr, &cpntr,
               matrix_type);

  printf("Processor %d: Completed AZ_transform\n",proc_config[AZ_node]) ;
      has_global_indices = 0;
      option = AZ_LOCAL;

#ifdef VBR
  N_local = rpntr[N_update];
#else
  N_local = N_update;
#endif
 
  Amat = AZ_matrix_create(N_local);

#ifdef VBR
  AZ_set_VBR(Amat, rpntr, cpntr, bpntr, indx, bindx, val, data_org,
          N_update, update, option);
#else
  AZ_set_MSR(Amat, bindx, val, data_org, N_update, update, option);
#endif


  printf("proc %d Completed AZ_create_matrix\n",proc_config[AZ_node]) ;


  if (has_global_indices)
    {
      N_external = 0;
    }
  else
    {
      N_external = data_org[AZ_N_external];
    }

  xsolve  = (double *) calloc(N_local + N_external, 
			   sizeof(double)) ;

  for (i=0; i<N_local; i++) xsolve[i] = x[i];

  /* Reorder rhs and xsolve to match matrix ordering from AZ_transform */
  if (!has_global_indices)
    {
      AZ_reorder_vec(b, data_org, update_index, rpntr) ;
      AZ_reorder_vec(xsolve, data_org, update_index, rpntr) ;
    }

#ifdef VBR
  AZ_check_vbr(N_update, data_org[AZ_N_ext_blk], AZ_LOCAL, 
	       bindx, bpntr, cpntr, rpntr, proc_config);
#else
  AZ_check_msr(bindx, N_update, N_external, AZ_LOCAL, proc_config);
#endif

  printf("Processor %d of %d N_local = %d N_external = %d NNZ = %d\n",
	 proc_config[AZ_node],proc_config[AZ_N_procs],N_local,N_external,
	 n_nonzeros);

  AZ_defaults(options, params);

  options[AZ_solver]  = AZ_gmres;
  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_ilut;
  params[AZ_ilut_fill] = 2.0;
 /*
  params[AZ_ilut_fill] = 2.0;
  params[AZ_drop] = 0.01;
  options[AZ_overlap] = 0;
  options[AZ_reorder] = 0;
  params[AZ_rthresh] = 1.0E-1;
  params[AZ_athresh] = 1.0E-1;
  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_bilu_ifp;
  options[AZ_reorder] = 0;
  options[AZ_graph_fill] = 0;
  params[AZ_rthresh] = 1.0E-7;
  params[AZ_athresh] = 1.0E-7;
 options[AZ_poly_ord] = 1;
 options[AZ_precond] = AZ_Jacobi;
  params[AZ_omega] = 1.0;
  options[AZ_precond] = AZ_none ;

  options[AZ_poly_ord] = 1;
  options[AZ_precond] = AZ_Jacobi ;
  options[AZ_scaling] = AZ_sym_row_sum ;
  options[AZ_scaling] = AZ_sym_diag;


  options[AZ_conv] = AZ_noscaled;
  options[AZ_scaling] = AZ_Jacobi ;

  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_icc ;
  options[AZ_subdomain_solve] = AZ_ilut ;
  params[AZ_omega] = 1.2;
  params[AZ_ilut_fill] = 2.0;
  params[AZ_drop] = 0.01;
  options[AZ_reorder] = 0;
  options[AZ_overlap] = 0;
  options[AZ_type_overlap] = AZ_symmetric;

  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_bilu ;
  options[AZ_graph_fill] = 0;
  options[AZ_overlap] = 0;

  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_bilu_ifp ;
  options[AZ_graph_fill] = 0;
  options[AZ_overlap] = 0;
  params[AZ_rthresh] = 1.0E-3;
  params[AZ_athresh] = 1.0E-3;

 options[AZ_poly_ord] = 1;
 options[AZ_precond] = AZ_Jacobi ;
 */


  options[AZ_kspace] = 600 ;

  options[AZ_max_iter] = 600 ;
  params[AZ_tol] = 1.0e-14;
  AZOO_iterate(xsolve, b, options, params, status, proc_config, Amat);

  if (proc_config[AZ_node]==0)
    {
      printf("True residual norm = %22.16g\n",status[AZ_r]);
      printf("True scaled res    = %22.16g\n",status[AZ_scaled_r]);
      printf("Computed res norm  = %22.16g\n",status[AZ_rec_r]);
    }

  /* Get solution back into original ordering */
  if (!has_global_indices)
      AZ_invorder_vec(xsolve, data_org, update_index, rpntr, x);
  else
    x = xsolve;


  if (xexact != NULL)
    {
      double sum = 0.0;
      double largest = 0.0;
      for (i=0; i<N_local; i++) sum += fabs(x[i]-xexact[i]);
 printf("Processor %d:  Difference between exact and computed solution = %12.4g\n",
	     proc_config[AZ_node],sum);
      for (i=0; i<N_local; i++) largest = AZ_MAX(largest,fabs(xexact[i]));
 printf("Processor %d:  Difference divided by max abs value of exact   = %12.4g\n",
	     proc_config[AZ_node],sum/largest);
    }

				       
  /*  NOTE:  This does not work because we need to use AZ_free for
      any arrays allocated by Aztec!!!

  free((void *) update);   free((void *) update_index);
  free((void *) external); free((void *) extern_index);

  free((void *) cpntr) ;
  free((void *) b) ;
  free((void *) x) ;
  free((void *) val) ;
  free((void *) rpntr) ;
    free((void *) bpntr) ;
    free((void *) bindx) ;
    free((void *) indx) ;

  */
#ifdef AZTEC_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return 0 ;
}
