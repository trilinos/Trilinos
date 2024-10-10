/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef ML_SELF_H
#define ML_SELF_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS)

int ML_Smoother_Self(ML_Smoother *sm,int inlen,double x[],int outlen,
                     double rhs[]);

void ML_Smoother_Clean_Self(void * Self_Handle);

extern int ML_Self_Solve(void * Self_Handle, double * x, double * rhs);

extern void ML_Self_Destroy(void * Self_Handle);

#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_SELF_H */
