/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the CG solver                                          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#ifndef __MLCG__
#define __MLCG__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_krylov.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int ML_CG_Solve(ML_Krylov *, int, double *, double *);
extern int ML_CG_ComputeEigenvalues(ML_Krylov *data, int length, int);
extern int ML_Power_ComputeEigenvalues(ML_Krylov *data, int length, int);
extern int ML_SubspaceIteration_ComputeEigenvalues(ML_Krylov *data,
                int length, int scale_by_diag);



#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif
