/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef __ML_AMESOS_H__
#define __ML_AMESOS_H__

#include <stdlib.h>
#include <assert.h>
#include "ml_struct.h"

#define ML_AMESOS_KLU            0
#define ML_AMESOS_UMFPACK        1
#define ML_AMESOS_SUPERLUDIST    2
#define ML_AMESOS_MUMPS          3
#define ML_AMESOS_SCALAPACK      4


#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

extern int ML_Smoother_Amesos(void *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Amesos(void *Amesos_Handle);

int ML_Gen_Smoother_Amesos(ML *ml, int nl, int AmesosSolver, int);


#ifndef ML_CPP
#ifdef __cplusplus
   }
#endif
#endif

#endif /* #ifndef __ML_AMESOS_H__ */
