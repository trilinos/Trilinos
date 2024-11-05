/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_SELF_WRAP
#define ML_SELF_WRAP

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

namespace Teuchos {
  class ParameterList;
}
class Epetra_Comm;

int ML_Smoother_Self(ML_Smoother *sm,int inlen,double x[],int outlen,
                     double rhs[]);

void ML_Smoother_Clean_Self(void * Self_Handle);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

int ML_Gen_Smoother_Self(ML *ml, int Overlap, int nl, int pre_or_post,
                         int niters, Teuchos::ParameterList& List,
                         const Epetra_Comm& Comm);

#endif
#endif
