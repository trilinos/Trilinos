

/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/


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
