

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

void AZ_ifpack_prec_destroy(int *options, double *params, int *proc_config,
		       AZ_MATRIX *Amat, AZ_PRECOND *Prec)
{
  AZ_IFPACK *Prec_pass_data;
  void *precon, *bmat;
  
  Prec_pass_data = (AZ_IFPACK *) Prec->Pmat->aux_ptr;

  precon = Prec_pass_data->precon;
  bmat = Prec_pass_data->bmat;
  options[AZ_precond] = Prec_pass_data->user_precon; /* Restore user prec*/



   /* Free allocated memory */


  ifp_freepreconditioner(precon);
  /*ifp_freeblockmatrix(bmat); Need to fix the destructor for BlockMat*/

  /* Must make sure to clean up everything!!!!*/

  AZ_free((void *) Prec_pass_data);

  }
