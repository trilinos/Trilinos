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

#ifndef __MLUTILH__
#define __MLUTILH__

#ifdef __cplusplus
   extern "C" {
#endif

   int ML_Coord2RBM(int Nnodes, double x[], double y[], double z[],
                    double rbm[], int Ndof);


#ifdef __cplusplus
}
#endif

#endif

