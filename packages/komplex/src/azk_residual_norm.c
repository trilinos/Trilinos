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

/*
 * This code is a convenience function for computing the two-norm of the 
 * residual of a system 
 *    b - A * x.  
 * It takes the real and imaginary parts of a complex
 *  valued linear system and returns a double.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "az_aztec.h"
#include "azk_komplex.h"

double AZK_residual_norm_no_copy(double *xr, double *xi, double *br, double *bi, 
				int *options, double *params,
				int *proc_config,
				AZ_MATRIX *Amat_real, AZ_MATRIX *Amat_imag)

/*******************************************************************************

  Author:          Mike Heroux, SNL, 9222
  =======

  Return code:     double
  ============

  Parameter list:
  ===============

  xr,xi:           On input, contains the initial guess, real part in xr and
                   imaginary part in xi.

  br,bi:           Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

 Amat_real,
 Amat_imag:        The real and imaginary parts of the complex operator, each
                   stored separately as AZ_MATRIX structures.

  Overview
  ========

  AZK_residual_norm_no_copy computes the two norm of the residual ||r|| where
  r = b - A*x.  Specifically, writing in terms of real and imaginary parts, 
  we have

  (rr + i*ri) = (br + i*bi) - (Ar + i*Ai)*(xr + i*xi).

  The two-norm of the complex vector r is identical to the two-norm of the
  twice-length real vector formed by concatenating rr = real(r) and 
  ri = imag(r).
  


*******************************************************************************/


{

  AZ_MATRIX  *Amat;   /* Structure representing matrix to be solved.          */
  double *x, *b;      /* Solution  and right-hand side to linear system.      */
  int N_equations, i;
  double *y_tmp, result;

  /* Transform complex system into komplex system */
  
  AZK_create_linsys_no_copy (xr,  xi,  br,  bi, options,  params,  
			    proc_config, Amat_real, Amat_imag, &x, &b, &Amat);
  
  /* Allocate temp vector y */
  
  N_equations = Amat->data_org[AZ_N_internal] + Amat->data_org[AZ_N_border];
  
  y_tmp = (double *) AZ_allocate(N_equations*sizeof(double));
  if (y_tmp == NULL) 
  AZ_perror("AZK_residual_norm_no_copy: Out of memory.");
  
  /* Compute y = A*x. */
  Amat->matvec(x, y_tmp, Amat, proc_config);

  /* Compute r = b - A*x (put in y_tmp) */
  /*daxpy_(&N_equations, &neg_one, b, &ione, y_tmp, &ione);*/

	for (i=0; i<N_equations; i++) y_tmp[i] = y_tmp[i] - b[i];

  /* Use Aztec function to compute norm */

  result = AZ_gvector_norm(N_equations, 2, y_tmp, proc_config);

  /* Free memory space */

  AZK_destroy_linsys (options,  params, proc_config, &x, &b, &Amat);
 
  AZ_free((void *) y_tmp);

  result = sqrt(result);
  return(result);
/* AZK_residual_norm */
}



double AZK_residual_norm(double *xk, double *bk, 
			 int *options, double *params,
			 int *proc_config,
			 AZ_MATRIX *Amat_komplex)

/*******************************************************************************

  Author:          Mike Heroux, SNL, 9222
  =======

  Return code:     double
  ============

  Parameter list:
  ===============

  xk:              On input, contains the initial guess.

  bk:              Right hand side of linear system.

  options:         Determines specific solution method and other parameters.

  params:          Drop tolerance and convergence tolerance info.

  proc_config:     Machine configuration.  proc_config[AZ_node] is the node
                   number.  proc_config[AZ_N_procs] is the number of processors.

 Amat_komplex:     The komplex operator, stored as an AZ_MATRIX structure.


  Overview
  ========

  AZK_residual_norm computes the two norm of the residual ||r|| where
  r = b - A*x.  Specifically, writing in terms of real and imaginary parts, 
  we have

  (rr + i*ri) = (br + i*bi) - (Ar + i*Ai)*(xr + i*xi).

  The two-norm of the complex vector r is identical to the two-norm of the
  twice-length real vector formed by concatenating rr = real(r) and 
  ri = imag(r).
  


*******************************************************************************/


{
  int N_equations, i;
  double *y_tmp, result;

  /* Allocate temp vector y */
  
  N_equations = Amat_komplex->data_org[AZ_N_internal] + 
                Amat_komplex->data_org[AZ_N_border];
  
  y_tmp = (double *) AZ_allocate(N_equations*sizeof(double));
  if (y_tmp == NULL) 
  AZ_perror("AZK_residual_norm: Out of memory.");
  
  /* Compute y = A*x. */
  Amat_komplex->matvec(xk, y_tmp, Amat_komplex, proc_config);

  /* Compute r = b - A*x (put in y_tmp) */
  /*daxpy_(&N_equations, &neg_one, bk, &ione, y_tmp, &ione);*/
	for (i=0; i<N_equations; i++) y_tmp[i] = y_tmp[i] - bk[i];

  /* Use Aztec function to compute norm */

  result = AZ_gvector_norm(N_equations, 2, y_tmp, proc_config);

  /* Free memory space */

  AZ_free((void *) y_tmp);

  result = sqrt(result);
  return(result);
/* AZK_residual_norm */
}



