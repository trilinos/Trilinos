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

/*! \file
\brief Destruction routine for deleting a preconditioner for a
Komplex matrix.

KOMPLEX is an add-on module to AZTEC that allows users to solve
complex-valued linear systems.  As such, all Aztec preconditioners are
available.  To learn how to set preconditioner options, please see the
Aztec 2.1 User Guide.

*/


/*! \fn void AZK_destroy_precon(int *options, double *params, int *proc_config,
		       AZ_MATRIX *Amat, AZ_PRECOND **Prec)

\brief Destroy a Komplex preconditioner.

Destroys a komplex preconditioner created by the
AZK_create_preconditioner function.  Deletes any memory allocated 
by creation routine.

\param options (In)
       Determines specific solution method and other parameters.
\param params (In)
       Drop tolerance and convergence tolerance info.
\param proc_config (In)
       Machine configuration.  proc_config[AZ_node] is the node
       number.  proc_config[AZ_N_procs] is the number of processors.
\param Amat (In)
       Komplex version of matrix stored as an AZ_MATRIX structure.
\param Prec (Out)
       Deleted komplex version of preconditioner stored as an
       AZ_PRECOND structure.

*/


void AZK_destroy_precon(int *options, double *params, int *proc_config,
		       AZ_MATRIX *Amat, AZ_PRECOND **Prec)
{
  AZ_KOMPLEX *Prec_pass_data;
  AZ_MATRIX *Pmat;
  double *val;
  int *bindx, *indx, *bpntr, *rpntr, *cpntr;

  Prec_pass_data = (AZ_KOMPLEX *) (*Prec)->Pmat->aux_ptr;
    
  if (Prec_pass_data->Form_of_Equations != AZK_Komplex_No_Copy)
    AZ_precond_destroy(Prec);
  else
    {
      Pmat  = (*Prec)->Pmat;
      val   = Pmat->val;
      bindx = Pmat->bindx;
      indx  = Pmat->indx;
      bpntr = Pmat->bpntr;
      rpntr = Pmat->rpntr;
      cpntr = Pmat->cpntr;
      
      /* Free allocated memory */
      if (Pmat != Amat) 
	{
	  if (val   != NULL) AZ_free ((void *) val);
	  if (bindx != NULL) AZ_free ((void *) bindx);
	  if (indx != NULL) AZ_free ((void *) indx);
	  if (bpntr != NULL) AZ_free ((void *) bpntr);
	  if (rpntr != NULL) AZ_free ((void *) rpntr);
	  
	  /* cpntr was created by Aztec, must be destroyed by Aztec */
	  if (cpntr != NULL) AZ_free ((void *) cpntr);
	  
	  AZ_free((void *) Prec_pass_data);
	  
	  /* Free data_org if Aztec doesn't do it */
	  if (!Pmat->must_free_data_org) AZ_free((void *) Pmat->data_org);
	  
	  AZ_matrix_destroy (&Pmat);
	}
      
      AZ_precond_destroy(Prec);
    }
}
