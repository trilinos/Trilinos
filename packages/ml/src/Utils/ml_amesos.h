/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef __ML_AMESOS_H__
#define __ML_AMESOS_H__

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"

#define ML_AMESOS_KLU 0
#define ML_AMESOS_UMFPACK 1
#define ML_AMESOS_SUPERLUDIST 2


extern int ML_Smoother_Amesos(void *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Amesos(void *Amesos_Handle);

int ML_Gen_Smoother_Amesos(ML *ml, int nl, int AmesosSolver, int);


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

extern int ML_Smoother_Ifpack(void *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

int ML_Gen_Smoother_Ifpack(ML *ml, int nl, int choice, int *options, double *params);

#endif /* #ifndef __ML_AMESOS_H__ */
