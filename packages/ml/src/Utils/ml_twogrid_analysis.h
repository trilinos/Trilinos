/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Some tools for two grid analysis.                                    */
/* ******************************************************************** */
/* Author        : Jonathan Hu (SNL)                                    */
/* Date          : September, 2001                                      */
/* ******************************************************************** */

#ifndef __MLTWOGRID__
#define __MLTWOGRID__

#include "ml_common.h"
#include "ml_defs.h"
#include "ml_struct.h"

#ifndef ML_CPP
#ifdef __cplusplus
   extern "C" {
#endif
#endif

extern double ML_gdot_H0(ML_Operator *Amat, double *vec1, double *vec2);
extern double ML_gdot_H1(ML_Operator *Amat, double *vec1, double *vec2);
extern double ML_gdot_H2(ML_Operator *Amat, double *vec1, double *vec2);
extern double ML_GetCoarseGridConst(ML_Operator *Amat, ML_Operator *Rmat,
                                    ML_Operator *Pmat, double *err_h);
extern double ML_GetSmoothingConst(ML_Operator *Amat, double *err_h,
                                    ML_Smoother *sm);
/*
extern double ML_GetTwoLevelConvergenceFactor(ML_Operator *Amat,
                                    ML_Operator *Rmat, ML_Operator *Pmat,
                                    ML_Smoother *sm,
                                    double *approx_soln, double *exact_soln);
*/
double ML_GetTwoLevelConvergenceFactor(ML *ml, double *approx_soln,
									   double *exact_soln);


#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /*ifdef __MLTWOGRID__*/

