
/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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

/*******************************************************************************
 * MATRIX FREE  matrix vector multiplication
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "az_ifpack.h"

void AZ_ifpack_precon(double x[], int options[], 
	int proc_config[], double params[], AZ_MATRIX *Amat,
        AZ_PRECOND *Prec)


/******************************************************************************/
/*
 * A bogus preconditioning subroutine which simply smooths x[]
 * in the interior of each local grid by taking averages with
 * neighboring grid points.
 *
 * Parameters:
 * =========
 * x                         On input, a vector. On output, x[] is 
 *                           smoothed by taking averages with neighbors.
 *
 * prec                      On input, prec->Pmat->aux_ptr points to 
 *                           that data_structure 'pass_data' which contains
 *                           the local grid size (nx,ny) on this processor.
 *
 * Amat, input_options,      Not used.
 * proc_config, Amat,
 * input_params
 *
 */


{
  int i, len;
  void *precon;
  AZ_IFPACK *Prec_pass_data;
  int nr, nc;
  double *input_vector;
                                    /* Data passing structure. This user-     */
                                    /* defined data structure is used to pass */
                                    /* information through Aztec and back into*/
                                    /* the user's subroutines.                */

  /*-------------------------------------------------------------------------*/
   /* Extract necessary data from pass_data */

  Prec_pass_data      = (AZ_IFPACK *) Prec->Pmat->aux_ptr;
  precon = (void *) Prec_pass_data->precon;
  nr = Prec_pass_data->nr;
  nc = Prec_pass_data->nc;
  if (nc != 1) abort();
  /* input_vector = (double *) Prec_pass_data->input_vector; */
  input_vector = (double *) malloc (nr * sizeof(double));
  len = nr*nc;
  /* dcopy_(&len, x, &ione, input_vector, &ione); */

  for (i=0; i<len; i++) input_vector[i] = x[i];

  ifp_apply(precon, nr, nc, input_vector, nr, x, nr);
  free((void *) input_vector);
}
