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

#include "ml_common.h"

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif


   int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[],
                    double rbm[], int Ndof);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

