/*******************************************************************************
 * Copyright 1995, Sandia Corporation.  The United States Government retains a *
 * nonexclusive license in this software as prescribed in AL 88-1 and AL 91-7. *
 * Export of this program may require a license from the United States         *
 * Government.                                                                 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "az_aztec.h"
#include "azk_komplex.h"

void AZK_precon(double x[], int options[], 
	int proc_config[], double params[], AZ_MATRIX *Amat,
        AZ_PRECOND *Prec)
{
   int i, j;
   double x_real, x_imag, *val, *val_real, *val_imag;
   AZ_MATRIX *Amat_real, *Amat_imag, *Pmat;
   int *data_org_real, N_real;
   int njacobi;

   AZ_KOMPLEX *pass_data, *Prec_pass_data;

   /* Extract necessary data from pass_data */


  Prec_pass_data      = (AZ_KOMPLEX *) Prec->Pmat->aux_ptr;


  pass_data      = (AZ_KOMPLEX *) Amat->aux_ptr;
  Amat_real = pass_data->Amat_real; /* Real operator in AZ_MATRIX struct */
  Amat_imag = pass_data->Amat_imag; /* Imag operator in AZ_MATRIX struct */
  val_real = Amat_real->val;
  val_imag = Amat_imag->val;

  data_org_real = Amat_real->data_org;
  N_real = data_org_real[AZ_N_internal] + data_org_real[AZ_N_border];
  
   /* set the preconditioning structure 'Prec'. */

   switch (Prec_pass_data->AZK_precond) {

   /* NO preconditioning. There is nothing to do                 */
   case AZ_none:
   case AZ_ls:                   
   case AZ_Neumann:
     break;

   /* Jacobi preconditioning. In this case, Aztec needs the      */
   /* diagonal of the matrix. This can be passed in as an MSR    */
   /* matrix. However, when using Jacobi, it is important to note*/
   /* that the MSR 'bindx' array does not need to be specified.  */
   /* Only the diagonal needs to be placed in the 'val' array.   */ 

   case AZ_Jacobi:
     /*case AZ_user_precond: */
    Pmat = Prec->Pmat;
     val = Pmat->val;
     njacobi = options[AZ_poly_ord];  /* Number of Jacobi steps */
     for (j = 0; j < njacobi; j++)
       for (i = 0; i<N_real; i++)
	 {
           x_real = x[2*i];
           x_imag = x[2*i+1];
           x[2*i]   = (val_real[i]*x_real + val_imag[i]*x_imag)*val[i];
	   x[2*i+1] = (val_real[i]*x_imag - val_imag[i]*x_real)*val[i];
	 }
     break;

   /* Domain decomposition preconditioning. In this case, Aztec  */
   /* needs the local matrix defined within the processor. This  */
   /* can be passed in as an MSR array. Note: since we do not    */
   /* overlap the subdomains in this specific example, only the  */
   /* local columns associated with local equations need to be   */
   /* kept. That is, we drop all references to external variables*/

  case AZ_dom_decomp:
      break;

  }
}
