/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/* ******************************************************************** */
/* Miscellaneous functions for efficient searching and sorting          */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL) and Raymond Tuminaro (SNL)       */
/* Date          : September, 1998                                      */
/* ******************************************************************** */

#ifndef __MLRBMH__
#define __MLRBMH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif


extern int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[], double rbm[], int Ndof, int NscalarDof);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif
