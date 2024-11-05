
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#ifndef _AZOO_ITERATE_H_
#define _AZOO_ITERATE_H_

#if defined(AztecOO_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The AztecOO package is deprecated"
#endif
#endif

#ifndef __cplusplus
#define __cplusplus
#endif

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "AztecOO.h"
#include "Aztec2Petra.h"

/*! \file 
\brief AZOO_iterate:  A function built around AztecOO that mimics the Aztec funciton AZ_iterate.

    AZOO_iterate is intended to facilitate the smooth transition from Aztec to Trilinos/AztecOO.
    The AZOO_iterate interface is essentially identical to the AZ_iterate interface and should be,
    for most uses a simple replacement in functionality.  

    However, because AztecOO uses Petra for
    distributed matrix and vector services (instead of AZ_MATRIX as defined by Aztec), there will 
    be some differences.  Some known differences are:

    <ol>
    <li> AZOO_iterate does not support Aztec's matrix-free version of AZ_MATRIX at this time.  
         Note that AztecOO has its own version of matrix-free implementation via the pure virtual
	 class Epetra_RowMatrix.
    <li> Scaling is performed independently in AztecOO.  All of the Aztec scaling options 
         (options[AZ_scaling]) are recognized, but block Jacobi scaling is implemented as point
	 Jacobi scaling.
    <li> Block entry preconditioners are not supported in AztecOO.  This functionality will be
         provided by IFPACK in a future release.

    </ol>

\warning {This function is intended as a temporary bridge for users migrating from Aztec to 
          Trilinos/AztecOO.  As such, it is not optimal in terms of memory use.  Further 
          flexibility can be had by copying AZOO_iterate() and customizing it to your own needs.
	  Ultimately, users will be best served by making a complete transition to the Trilinos/AztecOO
	  framework, building problems using Petra classes.
}
*/
/*! \fn void AZOO_iterate(double * xsolve, double * b,
                          int * options, double * params,
			     double * status, int *proc_config,
			     AZ_MATRIX * Amat,
			     AZ_PRECOND *precond, struct AZ_SCALING *scaling)

\brief Provides essentially equivalent functionality as the AZ_iterate function in Aztec 2.1.

*/


#ifdef __cplusplus
extern "C" void AZOO_iterate(double * xsolve, double * b,
			     int * options, double * params,
			     double * status, int *proc_config,
			     AZ_MATRIX * Amat,
			     AZ_PRECOND *precond, struct AZ_SCALING *scaling);
#endif
#endif /* _AZOO_ITERATE_H_ */
