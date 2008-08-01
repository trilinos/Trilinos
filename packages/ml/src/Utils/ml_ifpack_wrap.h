/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

#ifndef ML_IFPACK_WRAP
#define ML_IFPACK_WRAP

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif


int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

/** Generates the ifpack smoother */
int ML_Gen_Smoother_Ifpack(ML *ml, const char* Type, int Overlap,
                           int nl, int pre_or_post, 
                           void *List,
                           void *Comm);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
#endif
