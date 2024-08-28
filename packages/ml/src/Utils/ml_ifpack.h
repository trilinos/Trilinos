/* ************************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact person,   */
/* and disclaimer.                                                           */
/* ************************************************************************* */

#ifndef ML_IFPACK_H
#define ML_IFPACK_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_include.h"

typedef struct Ifpack_Handle_Struct Ifpack_Handle_Type;

struct Ifpack_Handle_Struct {
  void  *A_Base;        /* really Ifpack_Preconditioner pointer */
  int   freeMpiComm;  /*0 = false, 1 = true */
};

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS)

int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

/** Solves using IFPACK */
int ML_Ifpack_Solve(void * Ifpack_Handle, double * x, double * rhs);

/** Destroy all data associated to the IFPACK smoother. */
void ML_Ifpack_Destroy(void * Ifpack_Handle);

#endif

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_IFPACK_H */
