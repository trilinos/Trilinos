/*!
 * \file ml_ifpack_wrap.h
 *
 * \brief Interface to the Trilinos package Ifpack.
 *
 * The ML/Amesos interface allows ML users to apply Ifpack iterative methods
 * as smoothers.
 *
 *
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_IFPACK_WRAP
#define ML_IFPACK_WRAP

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

/** Apply the Ifpack smoother. */
int ML_Smoother_Ifpack(ML_Smoother *sm,int inlen,double x[],int outlen,
			      double rhs[]);

/** Destroy the Ifpack smoother. */
void ML_Smoother_Clean_Ifpack(void * Ifpack_Handle);

/** Generate the Ifpack smoother */
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
