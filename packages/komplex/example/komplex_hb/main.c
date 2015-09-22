/*
//@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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

#include "azk_komplex.h"
#include "az_aztec.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef AZTEC_MPI
#include "mpi.h"
#endif
#include "prototypes.h"

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }
#define double_quote '"'

int main(int argc, char *argv[])
{
  char  ri2k[]="ri2k";
  char  c2k[] ="c2k";
  char  g2k[] ="g2k";
  char  nc2k[]="nc2k";
  char  global[]="global";
  char  local[]="local";

  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    options[AZ_OPTIONS_SIZE]; /* Array used to select solver options.  */
  double params[AZ_PARAMS_SIZE];   /* User selected solver paramters.       */
  int    *data_org_mat0, *data_org_mat1;
                                   /* Array to specify data layout          */
  double status[AZ_STATUS_SIZE];   /* Information returned from AZ_solve(). */
  int    *update_mat0, *update_mat1;/* vector elements updated on this node. */
  int    *external_mat0, *external_mat1;
                                   /* vector elements needed by this node.  */
  int    *update_index_mat0, *update_index_mat1;
                                   /* ordering of update[] and external[]   */
  int    *extern_index_mat0, *extern_index_mat1;
                                   /* locally on this processor.            */
  int    *indx_mat0, *indx_mat1;   /* MSR/VBR format of real and imag parts */
  int    *bindx_mat0,*bindx_mat1;
  int    *bpntr_mat0,*bpntr_mat1;
  int    *rpntr_mat0,*rpntr_mat1;
  int    *cpntr_mat0,*cpntr_mat1;
  AZ_MATRIX *Amat_mat0, *Amat_mat1;
  double *val_mat0, *val_mat1;
  double *x_mat0, *b_mat0, *xexact_mat0, *xsolve_mat0;
  double *x_mat1, *b_mat1, *xexact_mat1, *xsolve_mat1;
  double *xsolve_complex, *b_complex, *val_complex, *xexact_complex;
  int    n_nonzeros_mat0, n_blk_nonzeros_mat0;
  int    n_nonzeros_mat1, n_blk_nonzeros_mat1;
  int    N_update_mat0, N_update_mat1;/* # of block unknowns updated on this node */
  int    N_local_mat0, N_local_mat1;
                                 /* Number scalar equations on this node */
  int    N_global_mat0, N_blk_global_mat0; /* Total number of equations */
  int    N_global_mat1, N_blk_global_mat1; /* Total number of equations */
  int    N_external_mat0;
  int    N_external_mat1;
  int    has_global_indices, option;
  int    is_VBR;

  double *val_msr_mat0, *val_msr_mat1;
  int   *bindx_msr_mat0, *bindx_msr_mat1;
  double residual, *xx_mat0, *xx_mat1, *xx;
  /*
  double c0r = 2.5;
  double c0i = 0.5;
  double c1r = 0.25;
  double c1i = - 0.3;
  */  
  double c0r = 1.0;
  double c0i = 0.0;
  double c1r = 0.0;
  double c1i = 1.0;
  AZ_MATRIX  *Amat;   /* Structure representing matrix to be solved.          */
  AZ_PRECOND *Prec;    /* Structure representing entire preconditioner.        */
                    /*                                                      */
  double *x, *b;      /* Initial guess and right-hand side to linear system.  */

  int matrix_type;
  int *ptmp0, *ptmp1;

  double * A0x0, * A0x1, * A1x0, * A1x1;

  int i, j;
  int ione = 1;
#ifdef AZTEC_MPI
  double MPI_Wtime(void) ;
#endif
  double time ;

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

#ifdef AZTEC_MPI
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

#ifdef VBRMATRIX
  if(argc != 6) 
#else
  if(argc != 5)
#endif
    {
      fprintf(stderr,
#ifdef VBRMATRIX
	     "Usage: %s c2k|g2k|ri2k local|global real_HB_filename imag_HB_filename partition_file\n",
#else
	     "Usage: %s c2k|g2k|ri2k local|global real_HB_filename imag_HB_filename\n",
#endif
	     argv[0]);
      fprintf(stderr,"where\n");
      fprintf(stderr,"Argument 1: Type of system to convert to komplex formulation: \n\n");

      fprintf(stderr,"  c2k -  Test program constructs a complex problem and calls\n");
      fprintf(stderr,"         appropriate komplex routines to solve it.\n\n");

      fprintf(stderr," ri2k -  Test program constructs a real/imaginary problem and calls\n");
      fprintf(stderr,"         appropriate komplex routines to solve it.\n\n");

      fprintf(stderr,"  g2k -  Test program constructs a general problem and calls\n");
      fprintf(stderr,"         appropriate komplex routines to solve it.\n\n");

      fprintf(stderr,"Argument 2: Type of indices:\n\n");

      fprintf(stderr,"  global - Construct input matrices using global index space.  The\n");
      fprintf(stderr,"           index values in bindx (the Aztec MSR/VBR index vector) refer to the\n");
      fprintf(stderr,"           global indexing of the matrix, not a tranformed index space used on\n");
      fprintf(stderr,"           each processor.  Global indices are what the user typically uses to \n");
      fprintf(stderr,"           construct their problem.\n");
      fprintf(stderr,"  local  - Construct input matrices using local index space.  The\n");
      fprintf(stderr,"           index values in bindx (the Aztec MSR/VBR index vector) refer to the \n");
      fprintf(stderr,"           local indexing of the matrix, typically the local index values are\n");
      fprintf(stderr,"           obtained by calling AZ_transform.\n\n");

      fprintf(stderr,"Argument 3: File name of a Harwell-Boeing format file containing the\n");
      fprintf(stderr,"            first matrix structure and values.  The exact use of this\n");
      fprintf(stderr,"            matrix will depend on the choice of Argument 1.\n\n");

      fprintf(stderr,"Argument 4: File name of a Harwell-Boeing format file containing the\n");
      fprintf(stderr,"            second matrix structure and values.The exact use of this\n\n");
      fprintf(stderr,"            matrix will depend on the choice of Argument 1.\n\n");

#ifdef VBRMATRIX
      fprintf(stderr,"Argument 5: The block size to use for partitioning the input matrices\n");
      fprintf(stderr,"            into variable block matrices.\n\n");
#endif

  exit(1);
    }




  /* Set some things to NULL */
  indx_mat0   = NULL;
  bpntr_mat0  = NULL;
  rpntr_mat0  = NULL;

  indx_mat1   = NULL;
  bpntr_mat1  = NULL;
  rpntr_mat1  = NULL;

  xexact_mat0 = NULL;
  xexact_mat1 = NULL;

  /* Read matrix file and distribute among processors.  
     Returns with this processor's set of rows */ 

#ifdef VBRMATRIX

  is_VBR = 1;

  /* Do real part first */

  read_hb(argv[3], proc_config, &N_global_mat0, &n_nonzeros_mat0, 
	  &val_msr_mat0,  &bindx_msr_mat0, &x_mat0, &b_mat0, &xexact_mat0);
  
  create_vbr(argv[5], proc_config, &N_global_mat0, &N_blk_global_mat0,
	     &n_nonzeros_mat0, &n_blk_nonzeros_mat0, 
	     &N_update_mat0, &update_mat0,
	     bindx_msr_mat0, val_msr_mat0, &val_mat0, &indx_mat0, 
	     &rpntr_mat0, &cpntr_mat0, &bpntr_mat0, &bindx_mat0);

  if(proc_config[AZ_node] == 0) 
    {
      free ((void *) val_msr_mat0);
      free ((void *) bindx_msr_mat0);
      free ((void *) cpntr_mat0);
    }
    matrix_type = AZ_VBR_MATRIX;

#ifdef AZTEC_MPI 
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

  distrib_vbr_matrix( proc_config, N_global_mat0, N_blk_global_mat0, 
		      &n_nonzeros_mat0, &n_blk_nonzeros_mat0,
		      &N_update_mat0, &update_mat0, 
		      &val_mat0, &indx_mat0, &rpntr_mat0, 
		      &cpntr_mat0, &bpntr_mat0, &bindx_mat0, 
		      &x_mat0, &b_mat0, &xexact_mat0);

  /* Now do Imaginary part */

  read_hb(argv[4], proc_config, &N_global_mat1, &n_nonzeros_mat1, 
	  &val_msr_mat1,  &bindx_msr_mat1, &x_mat1, &b_mat1, &xexact_mat1);
  
  create_vbr(argv[5], proc_config, &N_global_mat1, &N_blk_global_mat1,
	     &n_nonzeros_mat1, &n_blk_nonzeros_mat1, 
	     &N_update_mat1, &update_mat1,
	     bindx_msr_mat1, val_msr_mat1, &val_mat1, &indx_mat1, 
	     &rpntr_mat1, &cpntr_mat1, &bpntr_mat1, &bindx_mat1);

  if(proc_config[AZ_node] == 0) 
    {
      free ((void *) val_msr_mat1);
      free ((void *) bindx_msr_mat1);
      free ((void *) cpntr_mat1);
    }
    matrix_type = AZ_VBR_MATRIX;

#ifdef AZTEC_MPI
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

  distrib_vbr_matrix( proc_config, N_global_mat1, N_blk_global_mat1, 
		      &n_nonzeros_mat1, &n_blk_nonzeros_mat1,
		      &N_update_mat1, &update_mat1, 
		      &val_mat1, &indx_mat1, &rpntr_mat1, 
		      &cpntr_mat1, &bpntr_mat1, &bindx_mat1, 
		      &x_mat1, &b_mat1, &xexact_mat1);

#else

  is_VBR = 0;

  /* Do real part first */

    read_hb(argv[3], proc_config, &N_global_mat0, &n_nonzeros_mat0,
             &val_mat0, &bindx_mat0, &x_mat0, &b_mat0, &xexact_mat0);

#ifdef AZTEC_MPI
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

  distrib_msr_matrix(proc_config, N_global_mat0, &n_nonzeros_mat0, &N_update_mat0,
		     &update_mat0, &val_mat0, &bindx_mat0, 
		     &x_mat0, &b_mat0, &xexact_mat0);
#ifdef DEBUG
  for (i = 0; i<N_update_mat0; i++)
    if (val_mat0[i] == 0.0 ) printf("Zero diagonal at row %d\n",i);
#endif
    matrix_type = AZ_MSR_MATRIX;

    /* Do imaginary part next */

    read_hb(argv[4], proc_config, &N_global_mat1, &n_nonzeros_mat1,
             &val_mat1,  &bindx_mat1, &x_mat1, &b_mat1, &xexact_mat1);

#ifdef AZTEC_MPI
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

  distrib_msr_matrix(proc_config, N_global_mat1, &n_nonzeros_mat1, &N_update_mat1,
		     &update_mat1, &val_mat1, &bindx_mat1, 
		     &x_mat1, &b_mat1, &xexact_mat1);

#ifdef DEBUG
  for (i = 0; i<N_update_mat1; i++)
    if (val_mat1[i] == 0.0 ) printf("Zero diagonal at row %d\n",i);
#endif

  matrix_type = AZ_MSR_MATRIX;
#endif
  /* convert matrix to a local distributed matrix */
  if(!strcmp(argv[2],local))
    {
      AZ_transform(proc_config, &external_mat0, bindx_mat0, val_mat0, update_mat0,
		   &update_index_mat0, &extern_index_mat0, &data_org_mat0, 
		   N_update_mat0, indx_mat0, bpntr_mat0, rpntr_mat0, &cpntr_mat0,
		   matrix_type);
      
      AZ_transform(proc_config, &external_mat1, bindx_mat1, val_mat1, update_mat1,
		   &update_index_mat1, &extern_index_mat1, &data_org_mat1, 
		   N_update_mat1, indx_mat1, bpntr_mat1, rpntr_mat1, &cpntr_mat1,
		   matrix_type);
      
      printf("Processor %d: Completed AZ_transform\n",proc_config[AZ_node]) ;
      has_global_indices = 0;
      option = AZ_LOCAL;
    }
  else if(!strcmp(argv[2],global))
    {
      has_global_indices = 1;
      option = AZ_GLOBAL;
    }
  else
    perror("Must set global/local arg");
  
  if (!has_global_indices)
    {
      int ii,ndiff= 0;
  
      for (ii=0;ii<N_update_mat0; ii++)
	if (update_index_mat0[ii] != update_index_mat1[ii]) ndiff++;

      if (ndiff && strcmp(argv[1],g2k)) 
	{
	  fprintf(stderr,"Processor %d: Real and Imag equations ordered differently.\n",
		  proc_config[AZ_node]) ;
	  fprintf(stderr,"Number of difference is %d\n",ndiff) ;
	  perror("Not supported for this matrix form");
	}
    }

#ifdef VBRMATRIX
  N_local_mat0 = rpntr_mat0[N_update_mat0];
  N_local_mat1 = rpntr_mat1[N_update_mat1];
#else
  N_local_mat0 = N_update_mat0;
  N_local_mat1 = N_update_mat1;
#endif

  /* Real Part */

  Amat_mat0 = AZ_matrix_create(N_local_mat0);

#ifdef VBRMATRIX
  AZ_set_VBR(Amat_mat0, rpntr_mat0, cpntr_mat0, bpntr_mat0, 
             indx_mat0, bindx_mat0, val_mat0, data_org_mat0, 
	     N_update_mat0, update_mat0, option);
#else
  AZ_set_MSR(Amat_mat0, bindx_mat0, val_mat0, data_org_mat0, 
	     N_update_mat0, update_mat0, option);
#endif

  /* Imag Part */

  Amat_mat1 = AZ_matrix_create(N_local_mat1);

#ifdef VBRMATRIX
  AZ_set_VBR(Amat_mat1, rpntr_mat1, cpntr_mat1, bpntr_mat1, 
             indx_mat1, bindx_mat1, val_mat1, data_org_mat1, 
	     N_update_mat1, update_mat1, option);
#else
  AZ_set_MSR(Amat_mat1, bindx_mat1, val_mat1, data_org_mat1, 
	     N_update_mat1, update_mat1, option);
#endif

  printf("Processor %d Completed AZ_create_matrix\n",proc_config[AZ_node]) ;

#ifdef AZTEC_MPI
  MPI_Barrier(MPI_COMM_WORLD) ;
#endif

  /* initialize AZTEC options */
 
  AZ_defaults(options, params);
  options[AZ_solver]  = AZ_gmres;
  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_ilut;
  params[AZ_ilut_fill] = 8.0;
  options[AZ_graph_fill] = 1;
  options[AZ_overlap] = 1;
  params[AZ_drop] = 0.001;
  params[AZ_athresh] = 0.00003;
  params[AZ_rthresh] = 1.0;
  /*
  options[AZ_overlap] = 1;
  options[AZ_precond] = AZ_none ;

  options[AZ_poly_ord] = 1;
  options[AZ_precond] = AZ_Jacobi ;

 
  options[AZ_conv] = AZ_noscaled;
  options[AZ_scaling] = AZ_Jacobi ;

  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_ilut ;
  params[AZ_ilut_fill] = 2.0;
  params[AZ_drop] = 0.01;
  options[AZ_overlap] = 0;
  options[AZ_type_overlap] = AZ_symmetric;

  options[AZ_precond] = AZ_dom_decomp ;
  options[AZ_subdomain_solve] = AZ_bilu ;
  options[AZ_graph_fill] = 0;
  options[AZ_overlap] = 0;

 options[AZ_poly_ord] = 1;
 options[AZ_precond] = AZ_Jacobi ; */

  options[AZ_kspace] = 800 ;

  options[AZ_max_iter] = 800;
  params[AZ_tol] = 1.e-16 ;


  /* xsolve is a little longer vector needed to account for external 
     entries.  Make it and copy x (initial guess) into it.
  */
  if (has_global_indices)
    {
      N_external_mat0 = 0;
      N_external_mat1 = 0;
    }
  else
    {
      N_external_mat0 = data_org_mat0[AZ_N_external];
      N_external_mat1 = data_org_mat1[AZ_N_external];
    }

  xsolve_mat0  = (double *) calloc(N_local_mat0 + N_external_mat0, 
			   sizeof(double)) ;

  for (i=0; i<N_local_mat0; i++) xsolve_mat0[i] = x_mat0[i];

  xsolve_mat1  = (double *) calloc(N_local_mat1 + N_external_mat1, 
			   sizeof(double)) ;

  for (i=0; i<N_local_mat1; i++) xsolve_mat1[i] = x_mat1[i];

  /* Compute complex RHS such that, when solved by az_komplex, the solution
     matches the individual solution vectors for the real and imaginary parts.
     Precisely, we have xr and xi such that Ar*xr = br and Ai*xi = bi.  We got
     this data from the HB files.  Now we will use xr and xi to make a complex
     RHS bc such that the solution xc to (Ar + i*Ai) * xc = bc has
     real(xc) = xr and imag(xc) = xi.
  */

  if (!strcmp(argv[3],argv[4]))
    { double tmp1, tmp2;
    if(!strcmp(argv[1],g2k) || !strcmp(argv[1],ri2k)) 
      {
      for (i=0; i<N_local_mat0; i++)
	{
	  tmp1 = (c0r+c1r)*b_mat0[i] - (c0i+c1i)*b_mat1[i];
	  tmp2 = (c0i+c1i)*b_mat0[i] + (c0r+c1r)*b_mat1[i];
	  b_mat0[i] = tmp1;
	  b_mat1[i] = tmp2;
	}
    if(!strcmp(argv[1],ri2k))
	for (i=0; i<n_nonzeros_mat0+(1-is_VBR); i++)
	  {
	    val_mat0[i] *= (c0r+c1r);
	    val_mat1[i] *= (c0i+c1i);
	  }
      }
    else
      for (i=0; i<N_local_mat0; i++)
	{
	  tmp1 = b_mat0[i] - b_mat1[i];
	  tmp2 = b_mat0[i] + b_mat1[i];
	  b_mat0[i] = tmp1;
	  b_mat1[i] = tmp2;
	}
    }

  /* Copy exact solution for later use */

  if (xexact_mat0 != NULL)
    {
      xx_mat0  = (double *) calloc(N_local_mat0 + N_external_mat0, 
				   sizeof(double)) ;  
      for (i=0; i<N_local_mat0; i++) xx_mat0[i] = xexact_mat0[i];

      xx_mat1  = (double *) calloc(N_local_mat1 + N_external_mat1, 
				   sizeof(double)) ;
      for (i=0; i<N_local_mat1; i++) xx_mat1[i] = xexact_mat1[i];
    }
  

  /* Reorder rhs and xsolve and exact solution (if available) 
     to match matrix ordering from AZ_transform */

  if (!has_global_indices)
    {
      AZ_reorder_vec(b_mat0, data_org_mat0, update_index_mat0, rpntr_mat0) ;
      AZ_reorder_vec(xsolve_mat0, data_org_mat0, update_index_mat0, rpntr_mat0) ;
      
      AZ_reorder_vec(b_mat1, data_org_mat1, update_index_mat1, rpntr_mat1) ;
      AZ_reorder_vec(xsolve_mat1, data_org_mat1, update_index_mat1, rpntr_mat1) ;

   if (xexact_mat0 != NULL)
     {
       AZ_reorder_vec(xx_mat0, data_org_mat0, update_index_mat0, rpntr_mat0) ;
       AZ_reorder_vec(xx_mat1, data_org_mat1, update_index_mat1, rpntr_mat1) ;
     }

#ifdef VBRMATRIX
      AZ_check_vbr(N_update_mat0, data_org_mat0[AZ_N_ext_blk], AZ_LOCAL, 
		   bindx_mat0, bpntr_mat0, cpntr_mat0, rpntr_mat0, proc_config);

      AZ_check_vbr(N_update_mat1, data_org_mat1[AZ_N_ext_blk], AZ_LOCAL, 
		   bindx_mat1, bpntr_mat1, cpntr_mat1, rpntr_mat1, proc_config);
#else
      AZ_check_msr(bindx_mat0, N_update_mat0, N_external_mat0, AZ_LOCAL, 
		   proc_config);
      
      AZ_check_msr(bindx_mat1, N_update_mat1, N_external_mat1, AZ_LOCAL, 
		   proc_config);
#endif

      {int N_equations, N_external;

      N_equations = data_org_mat0[AZ_N_internal] + data_org_mat0[AZ_N_border];
      
      N_external = max(data_org_mat0[AZ_N_external],
		       data_org_mat1[AZ_N_external]);
      
      printf("Processor %d of %d N_equations = %d N_external = %d.\n",
	     proc_config[AZ_node],proc_config[AZ_N_procs],N_equations,N_external) ;
      
      }
    }

    /**************************************************************/
    /* Construct linear system.  Form depends on input parameters */
    /* Build exact solution vector also.                          */
    /**************************************************************/

    if (N_local_mat0 != N_local_mat1)
      perror("Error: Matrices must have same dimension\n");

#ifdef VBRMATRIX
	ptmp0 = bpntr_mat0;
	ptmp1 = bpntr_mat1;
#else
	ptmp0 = bindx_mat0;
	ptmp1 = bindx_mat1;
#endif

    if (!strcmp(argv[1],ri2k)) /* Real and imag have same structure */
      {

	for (i=0; i<N_update_mat0; i++)
	  {
	    if (ptmp0[i] != ptmp1[i])
	      {
		fprintf(stderr,"Processor %d:  pntr0[%d] = %d  pntr1[%d] = %d\n",
			proc_config[AZ_node],i,ptmp0[i],i,ptmp1[i]);
		perror("Error: ri2k must have matrices with same structure\n");
	      }
	    for (j=ptmp0[i]; j< ptmp0[i+1]; j++)
	      if (bindx_mat0[j] != bindx_mat1[j])
		{
		  fprintf(stderr,"Processor %d:  bindx0[%d] = %d  bindx1[%d] = %d\n",
			proc_config[AZ_node],j,bindx_mat0[j],j,bindx_mat1[j]);
		  perror("Error: ri2k must have matrices with same structure\n");
		}
	  }

	AZK_create_linsys_ri2k (xsolve_mat0,  xsolve_mat1,  b_mat0,  b_mat1, 
				options,  params, proc_config,  
				Amat_mat0, val_mat1, &x, &b, &Amat);

       if (xexact_mat0 != NULL)
	 AZK_create_vector_ri2k(options,  params, proc_config, Amat, 
				xx_mat0, xx_mat1, &xx);

      }
    else if (!strcmp(argv[1],c2k)) /* Complex matrix */
      {
	for (i=0; i<N_update_mat0; i++)
	  {
	    if (ptmp0[i] != ptmp1[i])
	      perror("Error: c2k must have matrices with same structure\n");
	    for (j=ptmp0[i]; j< ptmp0[i+1]; j++)
	      if (bindx_mat0[i] != bindx_mat1[i])
		perror("Error: c2k must have matrices with same structure\n");
	  }
	
	val_complex = (double *) calloc(2*(n_nonzeros_mat0+1), sizeof(double));
	xsolve_complex   = (double *) calloc(2*N_local_mat0, sizeof(double));
	b_complex   = (double *) calloc(2*N_local_mat0, sizeof(double));
	for (i=0; i<n_nonzeros_mat0+(1-is_VBR); i++)
	  {
	    val_complex[2*i  ] = val_mat0[i];
	    val_complex[2*i+1] = val_mat1[i];
	  }
	for (i=0; i<N_local_mat0; i++)
	  {
	    xsolve_complex[2*i  ] = xsolve_mat0[i];
	    xsolve_complex[2*i+1] = xsolve_mat1[i];
	    b_complex[2*i  ] = b_mat0[i];
	    b_complex[2*i+1] = b_mat1[i];
	  }

	Amat_mat0->val = val_complex;

	AZK_create_linsys_c2k (xsolve_complex,  b_complex, 
			       options,  params, proc_config,  
			       Amat_mat0, &x, &b, &Amat);

       if (xexact_mat0 != NULL)
	 {
	   xexact_complex   = (double *) calloc(2*N_local_mat0, sizeof(double));
	   for (i=0; i<N_local_mat0; i++)
	     {
	       xexact_complex[2*i  ] = xx_mat0[i];
	       xexact_complex[2*i+1] = xx_mat1[i];
	     }
	   
	   AZK_create_vector_c2k(options,  params, proc_config, Amat, 
			       xexact_complex, &xx);
	 }
      }
    else if(!strcmp(argv[1],g2k)) /* General matrix */
      {
	AZK_create_linsys_g2k (xsolve_mat0,  xsolve_mat1,  b_mat0,  b_mat1, 
			       options,  params, proc_config,  
			       c0r, c0i, Amat_mat0, 
			       c1r, c1i, Amat_mat1, &x, &b, &Amat);

       if (xexact_mat0 != NULL)
	 AZK_create_vector_ri2k(options,  params, proc_config, Amat, 
				xx_mat0, xx_mat1, &xx);

      }
   else if(!strcmp(argv[1],nc2k)) /* No copy matrix */
     {
       if (has_global_indices) perror("Global indices not supported for nc2k");

       AZK_create_linsys_no_copy (xsolve_mat0,  xsolve_mat1,  b_mat0,  b_mat1, 
				  options,  params,  
				  proc_config, Amat_mat0, Amat_mat1, &x, &b, &Amat);
     }
    else
      {
	fprintf(stderr, "Error: Matrix form set to: %s\n",argv[1]);
	perror("Error: This form not supported\n");
      }

    /**************************************************************/
    /* Check residual of init guess and exact solution            */
    /**************************************************************/

 
    if(!strcmp(argv[1],nc2k)) /* No copy matrix */
      residual = AZK_residual_norm_no_copy(xsolve_mat0, xsolve_mat1, 
					   b_mat0, b_mat1, 
					   options, params, proc_config,
					   Amat_mat0, Amat_mat1);
    else
      residual = AZK_residual_norm(x, b, options, params, proc_config, Amat);


    if (proc_config[AZ_node]==0)
      printf("\n\n\nExplicitly computed norm of residual using initial guess = %12.4g\n",residual);


    if (xexact_mat0 != NULL) {
  
      if(!strcmp(argv[1],nc2k)) /* No copy matrix */
	residual = AZK_residual_norm_no_copy(xx_mat0, xx_mat1, 
					     b_mat0, b_mat1, 
					     options, params, proc_config,
					     Amat_mat0, Amat_mat1);
      else {
	residual = AZK_residual_norm(xx, b, options, params, proc_config, Amat);
	AZK_destroy_vector(options, params, proc_config, Amat, &xx);
      }
    }
    if (proc_config[AZ_node]==0)
      printf("\n\n\nExplicitly computed norm of residual using exact solution = %12.4g\n",residual);

    /**************************************************************/
    /* Create preconditioner                                      */
    /**************************************************************/

   AZK_create_precon(options,  params, proc_config, x, b, Amat, &Prec);

    /**************************************************************/
    /* solve linear system using Aztec.                           */
    /**************************************************************/

   AZ_iterate(x, b, options, params, status, proc_config, Amat, Prec, NULL);

    /**************************************************************/
    /* Extract solution.                                          */
    /**************************************************************/

   if (!strcmp(argv[1],c2k)) /* Complex matrix */  
     AZK_extract_solution_k2c(options, params, proc_config, Amat, Prec, x, 
			       xsolve_complex); 
   else
     AZK_extract_solution_k2ri(options, params, proc_config, Amat, Prec, x, 
			       xsolve_mat0,  xsolve_mat1); 
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
  /* Get solution back into original ordering */
  if (!has_global_indices)
    {
      AZ_invorder_vec(xsolve_mat0, data_org_mat0, update_index_mat0, 
		      rpntr_mat0, x_mat0);
      AZ_invorder_vec(xsolve_mat1, data_org_mat1, update_index_mat1, 
		      rpntr_mat1, x_mat1);
    }
  else
    {
      free ((void *) x_mat0);
      free ((void *) x_mat1);
      x_mat0 = xsolve_mat0;
      x_mat1 = xsolve_mat1;
      xsolve_mat0 = NULL;
      xsolve_mat1 = NULL;
    }
  if (xexact_mat0 != NULL)
    {
      double sum = 0.0;
      if (!strcmp(argv[1],c2k)) /* Complex matrix */  
	for (i=0; i<N_local_mat0+N_local_mat1; i++) sum += fabs(xsolve_complex[i]-xexact_complex[i]);
      else
	{
	  for (i=0; i<N_local_mat0; i++) sum += fabs(x_mat0[i]-xexact_mat0[i]);
	  for (i=0; i<N_local_mat1; i++) sum += fabs(x_mat1[i]-xexact_mat1[i]);
	}
      printf("Processor %d:  Difference between exact and computed solution = %12.4g\n",
	     proc_config[AZ_node],sum);

      /* Test for null space if difference between computed and exact is large */
      /* Multiply the K matrix by the difference of the two solutions. */
      /* If the previous sum value is large and this one is small, the the matrix */
      /* is singular */
      if (strcmp(argv[1],c2k) && !strcmp(argv[2],local) && sum > 0.0e-5) {
   
	/* Compute difference between exact and computed solution */
	for (i=0; i<N_local_mat0; i++) x_mat0[i] = x_mat0[i]-xexact_mat0[i];
	for (i=0; i<N_local_mat1; i++) x_mat1[i] = x_mat1[i]-xexact_mat1[i];

	A0x0  = (double *) calloc(N_local_mat0, sizeof(double)) ;
	A0x1  = (double *) calloc(N_local_mat0, sizeof(double)) ;
	A1x0  = (double *) calloc(N_local_mat0, sizeof(double)) ;
	A1x1  = (double *) calloc(N_local_mat0, sizeof(double)) ;

	Amat_mat0->matvec(x_mat0, A0x0, Amat_mat0, proc_config);
	Amat_mat0->matvec(x_mat1, A0x1, Amat_mat0, proc_config);
	Amat_mat1->matvec(x_mat0, A1x0, Amat_mat1, proc_config);
	Amat_mat1->matvec(x_mat1, A1x1, Amat_mat1, proc_config);

	sum = 0.0;
	for (i=0; i<N_local_mat0; i++) sum += fabs(A0x0[i]-A1x1[i]);
	for (i=0; i<N_local_mat1; i++) sum += fabs(A1x0[i]+A0x1[i]);
	printf("Processor %d:  Difference between A* exact-computed solution = %12.4g\n",
	       proc_config[AZ_node],sum);

	

	free((void *) A0x0);
	free((void *) A0x1);
	free((void *) A1x0);
	free((void *) A1x1);
      }
    }

  /*  NOTE:  We have not free'd any memory.  Need to use AZ_free to
   *  free memory allocated by Aztec routines.
  */
  free ((void *) x_mat0);
  free ((void *) b_mat0);
  free ((void *) xexact_mat0);
  free ((void *) val_mat0);
  free ((void *) bindx_mat0);

  if(!strcmp(argv[2],local)) {
    AZ_free ((void *) external_mat0);
    AZ_free ((void *) update_index_mat0);
    AZ_free ((void *) extern_index_mat0);
    AZ_free ((void *) data_org_mat0);
  }
  if (xsolve_mat0!=NULL) free ((void *) xsolve_mat0);

  free ((void *) x_mat1);
  free ((void *) b_mat1);
  free ((void *) xexact_mat1);
  free ((void *) val_mat1);
  free ((void *) bindx_mat1);
  if(!strcmp(argv[2],local)) {
    AZ_free ((void *) external_mat1);
    AZ_free ((void *) update_index_mat1);
    AZ_free ((void *) extern_index_mat1);
    AZ_free ((void *) data_org_mat1);
  }
  if (xsolve_mat1!=NULL) free ((void *) xsolve_mat1);
  if (xexact_mat0 != NULL)
    {
      free ((void *) xx_mat0);  
      free ((void *) xx_mat1);
    }
  if (!strcmp(argv[1],c2k)) {
    free ((void *) val_complex);
    free ((void *) xsolve_complex);
    free ((void *) b_complex);
    if (xexact_complex != NULL) free ((void *) xexact_complex);
  }



#ifdef VBRMATRIX
  free ((void *)    indx_mat0);
  free ((void *)    rpntr_mat0);
  free ((void *)    bpntr_mat0);
  if(!strcmp(argv[2],local)) AZ_free ((void *) cpntr_mat0);

  free ((void *)    indx_mat1);
  free ((void *)    rpntr_mat1);
  free ((void *)    bpntr_mat1);
  if(!strcmp(argv[2],local)) AZ_free ((void *) cpntr_mat1);
#endif

  AZ_free((void *) update_mat0);
  AZ_free((void *) update_mat1);
  AZ_matrix_destroy(&Amat_mat0);
  AZ_matrix_destroy(&Amat_mat1);

#ifdef AZTEC_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return 0 ;
}
