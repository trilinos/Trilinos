/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef __MLCHECK__
#define __MLCHECK__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

#include "ml_common.h"

extern void ML_check_it(double sol[], double rhs[], ML *ml);

extern void ML_interp_check(ML *ml, int, int);

extern int  ML_Check(ML *ml);

extern int ML_Reitzinger_Check_Hierarchy(ML *ml, ML_Operator **Tmat_array, int incr_or_decr);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif
#endif
