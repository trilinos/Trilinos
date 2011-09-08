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
/*#define DEBUG*/
#ifdef DEBUG
extern double dnrm2_(int *n1, double *v1, int *dum11);
#endif 
void AZK_matvec_no_copy(double *x, double *y, AZ_MATRIX *Amat, 
			int proc_config[])



/******************************************************************************/
/*
 * Perform matrix vector operation:  
 * 
 *                    y = A x 
 *
 */


{

  AZ_MATRIX *Amat_real, *Amat_imag;
  int *data_org_real, *data_org_imag;
  int *komplex_to_real, *komplex_to_imag;
  int i;
  int N_equations, N_external;
   double *x_tmp, *y_tmp;

   AZ_KOMPLEX *pass_data;           /* Data passing structure. This user-     */
                                    /* defined data structure is used to pass */
                                    /* information through Aztec and back into*/
                                    /* the user's matrix-vector product.      */

#ifdef DEBUG
      double sum1, sum2;
      int int_one=1;
#endif

  /*-------------------------------------------------------------------------*/

   /* Extract necessary data from pass_data */

    pass_data      = (AZ_KOMPLEX *) Amat->aux_ptr;
    Amat_real = pass_data->Amat_real; /* Real operator in AZ_MATRIX struct */
    Amat_imag = pass_data->Amat_imag; /* Imag operator in AZ_MATRIX struct */
    data_org_real = Amat_real->data_org;
    data_org_imag = Amat_imag->data_org;

    /* Integer maps from komplex vector format to real and imag parts.
       Maps internal and border entries properly and allows for flexibility
       in how the real/imag parts are ordered for the preconditioner. */

    komplex_to_real = pass_data->komplex_to_real;
    komplex_to_imag = pass_data->komplex_to_imag;

    /* Number of real equations (both internal and border) must equal number
       of imag equations */

    N_equations = data_org_real[AZ_N_internal] + data_org_real[AZ_N_border];

    N_external = AZ_MAX(data_org_real[AZ_N_external],
		     data_org_imag[AZ_N_external]);

    /* Allocate space for temp copies of real/imag vectors. */

    y_tmp = (double *) AZ_allocate(N_equations*sizeof(double));
    x_tmp = (double *) AZ_allocate((N_equations+N_external)*sizeof(double));
    if (x_tmp == NULL) 
      AZ_perror("AZK_matvec_no_copy: Out of memory.");
 
    /* Extract real part of input vector x */

    for (i = 0; i < N_equations; i++)
      x_tmp[i] = x[komplex_to_real[i]];

     /* Multiply A_real * x_real */

    Amat_real->matvec(x_tmp, y_tmp, Amat_real, proc_config);

#ifdef DEBUG
      sum1 = dnrm2_(&N_equations,x_tmp,&int_one);
      sum2 = dnrm2_(&N_equations,y_tmp,&int_one);
      printf("Processor %d of %d norm xr = %16.12e  norm Arxr = %16.12e.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],sum1, sum2);

#endif
    /* Save to real part of output vector */

    for (i = 0; i < N_equations; i++)
      y[komplex_to_real[i]] = y_tmp[i];

    /* Multiply A_imag * x_real */

    Amat_imag->matvec(x_tmp, y_tmp, Amat_imag, proc_config);

#ifdef DEBUG
      sum1 = dnrm2_(&N_equations,x_tmp,&int_one);
      sum2 = dnrm2_(&N_equations,y_tmp,&int_one);
      printf("Processor %d of %d norm xr = %16.12e  norm Aixr = %16.12e.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],sum1, sum2);

#endif
    /* Save to imag part of output vector */

    for (i = 0; i < N_equations; i++)
      y[komplex_to_imag[i]] = y_tmp[i];

    /* Extract imag part of input vector x */

    for (i = 0; i < N_equations; i++)
      x_tmp[i] = x[komplex_to_imag[i]];

    /* Multiply A_real * x_imag */

    Amat_real->matvec(x_tmp, y_tmp, Amat_real, proc_config);

#ifdef DEBUG
      sum1 = dnrm2_(&N_equations,x_tmp,&int_one);
      sum2 = dnrm2_(&N_equations,y_tmp,&int_one);
      printf("Processor %d of %d norm xi = %16.12e  norm Arxi = %16.12e.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],sum1, sum2);

#endif
    /* Update imag part of output vector */

    for (i = 0; i < N_equations; i++)
      y[komplex_to_imag[i]] += y_tmp[i];

    /* Multiply A_imag * x_imag */

    Amat_imag->matvec(x_tmp, y_tmp, Amat_imag, proc_config);

#ifdef DEBUG
      sum1 = dnrm2_(&N_equations,x_tmp,&int_one);
      sum2 = dnrm2_(&N_equations,y_tmp,&int_one);
      printf("Processor %d of %d norm xi = %16.12e  norm Aixi = %16.12e.\n",
	       proc_config[AZ_node],proc_config[AZ_N_procs],sum1, sum2);

#endif
    /* Update real part of output vector */

    for (i = 0; i < N_equations; i++)
      y[komplex_to_real[i]] -= y_tmp[i];

    AZ_free((void *) x_tmp);
    AZ_free((void *) y_tmp);
}
