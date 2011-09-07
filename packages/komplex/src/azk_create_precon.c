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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "azk_komplex.h"
/*! \file
\brief Creation routines for constructing a preconditioner for a
Komplex matrix.

KOMPLEX is an add-on module to AZTEC that allows users to solve
complex-valued linear systems.  As such, all Aztec preconditioners are
available.  To learn how to set preconditioner options, please see the
Aztec 2.1 User Guide.

*/


/*! \fn void AZK_create_precon(int *options, double *params,
		       int *proc_config,double *x, double *b,
		       AZ_MATRIX *Amat, AZ_PRECOND **Prec)
\brief Create a Preconditioner for a Komplex matrix.

Constructs a preconditioner for a Komplex matrix Amat.  All
preconditioning options available in Aztec are supported.


\param options (In)
       Determines specific preconditioner method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.

\param x (In/Out)
       Komplex version of initial guess and solution.  May be modified
       depending on preconditioner options.
\param b (In/Out)
       Komplex version of RHS.  May be modified
       depending on preconditioner options.

\param Amat (In)
       Komplex version of matrix stored as an AZ_MATRIX structure.

\param Prec (Out)
       Preconditioner for Amat stored as an AZ_PRECOND structure.

*/

void AZK_create_precon(int *options, double *params,
		       int *proc_config,double *x, double *b,
		       AZ_MATRIX *Amat, AZ_PRECOND **Prec)
{
  AZ_KOMPLEX *pass_data, *Prec_pass_data;
  AZ_MATRIX *Pmat, *Amat_real, *Amat_imag;
  int N_equations, N_real;
  int *data_org_real, *data_org_imag;
  double *val;
  int i, is_VBR;

  /* Extract necessary data from pass_data */

  pass_data      = (AZ_KOMPLEX *) Amat->aux_ptr;

  if (pass_data->Form_of_Equations != AZK_Komplex_No_Copy)
    (*Prec) = AZ_precond_create(Amat,AZ_precondition,NULL);
  else
    {
      Amat_real = pass_data->Amat_real; /* Real operator in AZ_MATRIX struct */
      Amat_imag = pass_data->Amat_imag; /* Imag operator in AZ_MATRIX struct */

      data_org_real = Amat_real->data_org;
      data_org_imag = Amat_imag->data_org;
      N_real = data_org_real[AZ_N_internal] + data_org_real[AZ_N_border];
  
      N_equations = 2 * N_real;
      if (data_org_real[AZ_matrix_type] == AZ_VBR_MATRIX &&
	  data_org_imag[AZ_matrix_type] == AZ_VBR_MATRIX )
	{
	  is_VBR = 1;
	}
      else if (data_org_real[AZ_matrix_type] == AZ_MSR_MATRIX &&
	       data_org_imag[AZ_matrix_type] == AZ_MSR_MATRIX)
	{
	  is_VBR = 0;
	}
      else
	{
	  printf("Unsupported Matrix types\n");
	  abort();
	}

      /* set the preconditioning structure 'Prec'. */

      Prec_pass_data = (AZ_KOMPLEX *) AZ_allocate(sizeof(AZ_KOMPLEX));
      if (Prec_pass_data == NULL)
      AZ_perror("AZK_create_precon: Out of memory.");


      switch (options[AZ_precond]) {

	/* NO preconditioning. There is nothing to do                 */
	/* We just need to give a valid matrix to the preconditioning */
      case AZ_none:
	(*Prec) = AZ_precond_create(Amat,AZ_precondition,NULL);
	break;

	/* Polynomial preconditioning (least-squares or Neumann).     */
	/* Here we must give Aztec an upper bound for the norm of the */
	/* matrix. In addition, we need to tell Aztec to use the      */
	/* USER's matrix-vector product when applying the polynomial  */
	/* preconditioner.                                            */

      case AZ_ls:                   
      case AZ_Neumann:
	Amat->matrix_norm = 8.0;   
	(*Prec) = AZ_precond_create(Amat,AZ_precondition,NULL);
	break;

	/* Jacobi preconditioning. In this case, Aztec needs the      */
	/* diagonal of the matrix. This can be passed in as an MSR    */
	/* matrix. However, when using Jacobi, it is important to note*/
	/* that the MSR 'bindx' array does not need to be specified.  */
	/* Only the diagonal needs to be placed in the 'val' array.   */ 

      case AZ_Jacobi:
	if (!is_VBR)
	  {
	    Pmat = AZ_create_matrix(N_equations, AZ_NO_EXTRA_SPACE,
				    AZ_MSR_MATRIX, N_equations,
				    AZ_NOT_USING_AZTEC_MATVEC);
	    val      = (double *) AZ_allocate(N_real * sizeof(double));
	    if (val == NULL) 
	      AZ_perror("AZK_create_precon: Out of memory");
	    for ( i = 0;i < N_real; i++)
	      val[i] = 1.0/(Amat_real->val[i]*Amat_real->val[i] +
			    Amat_imag->val[i]*Amat_imag->val[i]);
	    Pmat->val = val;
	    Pmat->bindx = NULL;
	    Pmat->indx = NULL;
	    Pmat->bpntr = NULL;
	    Pmat->rpntr = NULL;
	    Pmat->cpntr = NULL;
	    (*Prec) = AZ_precond_create(Pmat,AZK_precon,NULL);
	    options[AZ_precond] = AZ_user_precond;
	    Prec_pass_data->AZK_precond = AZ_Jacobi;
	    (*Prec)->Pmat->aux_ptr = (void *) Prec_pass_data;
	  }
	else
	  {
	    AZ_perror("Jacobi scaling is only supported for MSR matrices");
	  }
	break;

	/* Domain decomposition preconditioning. In this case, Aztec  */
	/* needs the local matrix defined within the processor. This  */
	/* can be passed in as an MSR array. Note: since we do not    */
	/* overlap the subdomains in this specific example, only the  */
	/* local columns associated with local equations need to be   */
	/* kept. That is, we drop all references to external variables*/

      case AZ_dom_decomp:

	AZ_perror("AZK_linsys_create_no_copy does not support dom_decomp");

	break;
      default:
	AZ_perror("AZK_linsys_create_no_copy does not support this option");

      }
    }
}
