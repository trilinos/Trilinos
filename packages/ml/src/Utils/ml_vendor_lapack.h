#ifndef __MLVENDORLAPACK__
#define __MLVENDORLAPACK__





#ifdef USE_VENDOR_BLAS
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
#endif

#ifdef USE_VENDOR_LAPACK

#ifndef ML_DTRSV_FUNC
#define ML_DTRSV_FUNC
#endif
#ifndef ML_DGETRF_FUNC
#define ML_DGETRF_FUNC
#endif
#ifndef ML_DGETRS_FUNC
#define ML_DGETRS_FUNC
#endif
#ifndef ML_XERBLA_FUNC
#define ML_XERBLA_FUNC
#endif
#ifndef ML_LSAME_FUNC
#define ML_LSAME_FUNC
#endif
#ifndef ML_DLASWP_FUNC
#define ML_DLASWP_FUNC
#endif
#ifndef ML_DGETF2_FUNC
#define ML_DGETF2_FUNC
#endif
#ifndef ML_DGER_FUNC
#define ML_DGER_FUNC
#endif
#ifndef ML_ILAENV_FUNC
#define ML_ILAENV_FUNC
#endif
#ifndef ML_DGEMV_FUNC
#define ML_DGEMV_FUNC
#endif
#ifndef ML_DTRMV_FUNC
#define ML_DTRMV_FUNC
#endif
#ifndef ML_DLAMCH_FUNC
#define ML_DLAMCH_FUNC
#endif
#ifndef ML_DLAMC1_FUNC
#define ML_DLAMC1_FUNC
#endif
#ifndef ML_DLAMC2_FUNC
#define ML_DLAMC2_FUNC
#endif
#ifndef ML_DLAMC3_FUNC
#define ML_DLAMC3_FUNC
#endif
#ifndef ML_DLAMC4_FUNC
#define ML_DLAMC4_FUNC
#endif
#ifndef ML_DLAMC5_FUNC
#define ML_DLAMC5_FUNC
#endif

#ifndef TFLOP
#define ML_DPOTRS_FUNC
#define ML_DGELS_FUNC
#define ML_DGELQ2_FUNC
#define ML_DGELQF_FUNC
#define ML_DLASCL_FUNC
#define ML_DLASET_FUNC
#define ML_DLASSQ_FUNC
#define ML_DLANGE_FUNC
#define ML_DLABAD_FUNC
#define ML_DORMQR_FUNC
#define ML_DORMLQ_FUNC
#define ML_DORM2R_FUNC
#define ML_DORML2_FUNC
#endif

#endif
#endif
