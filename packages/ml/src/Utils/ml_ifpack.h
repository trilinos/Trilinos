/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef __ML_IFPACK_H__
#define __ML_IFPACK_H__

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"


#define ML_IFPACK_ICT                     0
#define ML_IFPACK_RICK                    1
#define ML_IFPACK_RILUK                   2
#define ML_IFPACK_JACOBI                 10
/* options for ML_Gen_Smoother_Ifpack */
#define ML_IFPACK_OPTIONS_SIZE           10
#define ML_IFPACK_OVERLAP                 0
#define ML_IFPACK_LEVEL_OF_FILL           1
/* params for ML_Gen_Smoother_Ifpack */
#define ML_IFPACK_PARAMS_SIZE            20
#define ML_IFPACK_RELAX_VALUE             0

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

extern int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

int ML_Gen_Smoother_Ifpack(ML *ml, int nl, int choice, int *options,
			   double *params);

#ifndef ML_CPP
#ifdef __cplusplus
   }
#endif
#endif

#endif /* #ifndef __ML_IFPACK_H__ */
