/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef ML_IFPACK_H
#define ML_IFPACK_H

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"

#define ML_IFPACK_OPTIONS_SIZE 10
#define ML_IFPACK_PARAMS_SIZE 5

#define ML_IFPACK_TYPE           0
#define ML_IFPACK_OVERLAP        1
#define ML_IFPACK_LOCAL_PARTS    2
#define ML_IFPACK_SWEEPS         3
#define ML_IFPACK_BLOCK_OVERLAP  4
#define ML_IFPACK_LEVEL_OF_FILL  5

#define ML_IFPACK_DAMPING_FACTOR 0

#define ML_IFPACK_AMESOS               0
#define ML_IFPACK_JACOBI               1
#define ML_IFPACK_GS                   2     
#define ML_IFPACK_SGS                  3
#define ML_IFPACK_BLOCK_JACOBI         4
#define ML_IFPACK_BLOCK_GS             5
#define ML_IFPACK_BLOCK_SGS            6
#define ML_IFPACK_BLOCK_JACOBI_AMESOS  7
#define ML_IFPACK_BLOCK_GS_AMESOS      8
#define ML_IFPACK_BLOCK_SGS_AMESOS     9
#define ML_IFPACK_ICT                 10
#define ML_IFPACK_RILUK               11

/* sets default parameters for IFPACK smoothers */

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

extern int ML_Ifpack_Defaults(int options[], double params[]);

extern int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

int ML_Gen_Smoother_Ifpack(ML *ml, int nl, int pre_or_post, int *options,
			   double *params);

#ifndef ML_CPP
#ifdef __cplusplus
   }
#endif
#endif

#endif /* #ifndef __ML_IFPACK_H__ */
