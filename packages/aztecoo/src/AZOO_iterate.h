
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _AZOO_ITERATE_H_
#define _AZOO_ITERATE_H_

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
