/*@HEADER
// ***********************************************************************
// 
//                Komplex: Complex Linear Solver Package
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
