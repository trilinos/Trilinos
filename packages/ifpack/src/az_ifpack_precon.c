
/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
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
