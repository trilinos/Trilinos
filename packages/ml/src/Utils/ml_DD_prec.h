/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_DD_PREC_H
#define ML_DD_PREC_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"
#include "ml_struct.h"
#include "ml_smoother.h"
#include "ml_defs.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

double ML_DD_OneLevel(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
double ML_DD_Additive(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
double ML_DD_Hybrid(ML_1Level *curr, double *sol, double *rhs,
			   int approx_all_zeros, ML_Comm *comm,
			   int res_norm_or_not, ML *ml);
double ML_DD_Hybrid_2(ML_1Level *curr, double *sol, double *rhs,
			     int approx_all_zeros, ML_Comm *comm,
			     int res_norm_or_not, ML *ml);
int ML_Aggregate_Stats_CleanUp_Info( ML *ml);
int ML_Aggregate_Stats_ComputeCoordinates( ML *ml, ML_Aggregate *ag,
					  double *x, double *y, double *z);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* ifdef ML_DD_PREC_H */
