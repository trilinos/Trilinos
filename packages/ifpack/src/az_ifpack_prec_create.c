

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
void AZ_ifpack_prec_create(double *x, double *b,
			   int *options, double *params,
			   int *proc_config,
			   AZ_MATRIX *Amat, AZ_PRECOND **Prec)
{
  AZ_IFPACK  *Prec_pass_data;
  void *precon, *bmat ;
  int nr, nc, *data_org;
  double rthresh, athresh;
 
  Prec_pass_data = (AZ_IFPACK *) AZ_allocate(sizeof(AZ_IFPACK));
  az2ifp_blockmatrix(&bmat, Amat); /* Create IFPACK encapsulation of Amat */
  
   /* set the preconditioning structure 'Prec'. */

   if  (options[AZ_precond] == AZ_none)
     ifp_preconditioner(&precon, bmat, IFP_NONE, 
		      (double) options[AZ_graph_fill], 0.0,
		      IFP_INVERSE, 0.0, 0.0);
 
   else if (options[AZ_precond] == AZ_Jacobi)
     {
       rthresh = params[AZ_rthresh];
       athresh = params[AZ_athresh];
       ifp_preconditioner(&precon, bmat, IFP_BJACOBI, 0.0, 0.0,
			IFP_SVD, rthresh, athresh);
       /*IFP_INVERSE, 0.0, 0.0); */
     }

   else if (options[AZ_precond] == AZ_dom_decomp && 
	    options[AZ_subdomain_solve] == AZ_bilu_ifp)
     {
       rthresh = params[AZ_rthresh];
       athresh = params[AZ_athresh];
       ifp_preconditioner(&precon, bmat, 
			IFP_BILUK, (double) options[AZ_graph_fill], 0.0,
			IFP_SVD, rthresh, athresh);
       /*IFP_INVERSE, 0.0, 0.0); */
       
     }
   else
     {
       printf("Not a supported preconditioner in az_ifpack_prec_create\n");
       abort();
     }

    (*Prec) = AZ_precond_create(Amat,AZ_ifpack_precon,NULL);


   /* Store pointers to preconditioner and IFPACK encapsulation of Amat */
  Prec_pass_data->precon = precon;
  Prec_pass_data->bmat = bmat;  

  /* Construct auxiliary vector for use with apply function.
     NOTE:  We are assuming only one RHS at this time !!! */

  data_org = Amat->data_org;
  nr = data_org[AZ_N_internal] + data_org[AZ_N_border];
  nc = 1;
  /*input_vector = (double *) malloc (nr * sizeof(double));
    Prec_pass_data.input_vector = input_vector; */
  Prec_pass_data->nr = nr;
  Prec_pass_data->nc = nc;
  (*Prec)->Pmat         = Amat;
  Prec_pass_data->user_aux_ptr = (*Prec)->Pmat->aux_ptr; /* Save this to be able to restore*/
  (*Prec)->Pmat->aux_ptr = (void *) Prec_pass_data;
  (*Prec)->prec_function = AZ_ifpack_precon;
  Prec_pass_data->user_precon = options[AZ_precond]; /* Save this to be able to restore*/
  options[AZ_precond] = AZ_user_precond;
}
