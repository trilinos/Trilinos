

/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
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
