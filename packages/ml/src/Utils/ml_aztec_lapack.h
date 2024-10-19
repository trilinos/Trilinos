/********************************************************************* */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

/********************************************************************* */
/*          BLAS/LAPACK Utilities for Aztec/ML users                   */
/********************************************************************* */

#ifndef __MLAZTECLAPACK__
#define __MLAZTECLAPACK__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif



#ifdef AZTEC
#include "ml_common.h"
#include "ml_defs.h"
#define ML_IDAMAX_FUNC
#define ML_DSWAP_FUNC
#define ML_DSCAL_FUNC
#define ML_DAXPY_FUNC
#define ML_DASUM_FUNC
#define ML_DDOT_FUNC
#define ML_DNRM2_FUNC
#define ML_DCOPY_FUNC
#define ML_DGEMM_FUNC
#define ML_DTRSM_FUNC
#define ML_DTRMM_FUNC
#define ML_DGETRS_FUNC
#define ML_LSAME_FUNC
#define ML_XERBLA_FUNC
#define ML_DLASWP_FUNC
#define ML_DGEMV_FUNC
#define ML_DGETRF_FUNC
#define ML_DGER_FUNC
#define ML_DTRMV_FUNC
#define ML_DTRSV_FUNC

#ifndef FSUB_TYPE
#  if defined(ncube)
#     define  FSUB_TYPE void
#  elif defined(paragon)
#     define  FSUB_TYPE void
#  elif defined(hp)
#     define  FSUB_TYPE void
#  else
#     define  FSUB_TYPE int
#  endif
#endif


#endif
#endif
