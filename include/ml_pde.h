/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Declaration of the ML_PDE functions                                  */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : February, 2000                                       */
/* ******************************************************************** */

#ifndef _MLPDE__
#define _MLPDE__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "mli_solver.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int ML_PDE_GenMat(MLI_Solver*,int);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
