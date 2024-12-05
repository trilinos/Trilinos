/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLSUPERLULAPACKH__
#define __MLSUPERLULAPACKH__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(SUPERLU)

#include "ml_common.h"

#ifdef SGI
#define ML_DLAMCH_FUNC
#define ML_DLAMC1_FUNC
#define ML_DLAMC2_FUNC
#define ML_DLAMC3_FUNC
#define ML_DLAMC4_FUNC
#define ML_DLAMC5_FUNC
#endif
#endif
#endif
