#ifndef _AZOO_ITERATE_H_
#define _AZOO_ITERATE_H_

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#ifdef PETRA_MPI
#include "mpi.h"
#define AZ_MPI
#define AZTEC_MPI
#endif
#include "az_aztec.h"
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Petra_Comm.h"
#include "Petra_Map.h"
#include "Petra_BlockMap.h"
#include "Petra_RDP_Vector.h"
#include "Petra_RDP_VBR_Matrix.h"
#include "Petra_RDP_CRS_Matrix.h"
#include "Petra_RDP_RowMatrix.h"
#include "Petra_RDP_LinearProblem.h"
#include "Aztec_OO.h"
#include "Aztec2Petra.h"

/*! \file 
\brief AZOO_iterate:  A function built around Aztec_OO that mimics the Aztec funciton AZ_iterate.

    AZOO_iterate is intended to facilitate the smooth transition from Aztec to Trilinos/Aztec_OO.
    The AZOO_iterate interface is essentially identical to the AZ_iterate interface and should be,
    for most uses a simple replacement in functionality.  

    However, because Aztec_OO uses Petra for
    distributed matrix and vector services (instead of AZ_MATRIX as defined by Aztec), there will 
    be some differences.  Some known differences are:

    <ol>
    <li> AZOO_iterate does not support Aztec's matrix-free version of AZ_MATRIX at this time.  
         Note that Aztec_OO has its own version of matrix-free implementation via the pure virtual
	 class Petra_RDP_RowMatrix.
    <li> Scaling is performed independently in Aztec_OO.  All of the Aztec scaling options 
         (options[AZ_scaling]) are recognized, but block Jacobi scaling is implemented as point
	 Jacobi scaling.
    <li> Block entry preconditioners are not supported in Aztec_OO.  This functionality will be
         provided by IFPACK in a future release.

    </ol>

\warning {This function is intended as a temporary bridge for users migrating from Aztec to 
          Trilinos/Aztec_OO.  As such, it is not optimal in terms of memory use.  Further 
          flexibility can be had by copying AZOO_iterate() and customizing it to your own needs.
	  Ultimately, users will be best served by making a complete transition to the Trilinos/Aztec_OO
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
