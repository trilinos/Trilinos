/*      LAPACK++ (V. 1.0 Beta)						*/
/*      (C) 1992-1994 All Rights Reserved.				*/
/*
              LAPACK++ 1.0: Linear Algebra Package 1.0
               University of Tennessee, Knoxvilee, TN.
            Oak Ridge National Laboratory, Oak Ridge, TN.
        Authors: J. J. Dongarra, E. Greaser, R. Pozo, D. Walker
                 (C) 1992-1993 All Rights Reserved

                             NOTICE

 Permission to use, copy, modify, and distribute this software and
 its documentation for any purpose and without fee is hereby granted
 provided that the above copyright notice appear in all copies and
 that both the copyright notice and this permission notice appear in
 supporting documentation.

 Neither the Institutions (University of Tennessee, and Oak Ridge National
 Laboratory) nor the Authors make any representations about the suitability
 of this software for any purpose.  This software is provided ``as is''
 without express or implied warranty.

 LAPACK++ was funded in part by the U.S. Department of Energy, the
 National Science Foundation and the State of Tennessee.
*/

//		Double precision Lapack routines

#ifndef _IFP_DLAPACK_H_
#define _IFP_DLAPACK_H_

#include "Ifpack_config.h"

#include "ifp_arch.h"

extern "C"
{
// *************************** Utility Routines **********************


    IFPACK_DEPRECATED double F77NAME(dlamch)(char *t);



//******************** Linear Equation Solvers *************************
    IFPACK_DEPRECATED void F77NAME(dgesv)(integer *n, integer *k, doublereal *A, integer *lda, integer *ipiv,
            doublereal *X, integer *ldx, integer *info);

    IFPACK_DEPRECATED void F77NAME(dposv)(char *uplo, integer *m, integer *k , doublereal *A, integer *lda,
        doublereal *X, integer *ldx, integer *info);

	IFPACK_DEPRECATED void F77NAME(dsysv)(const char *uplo, integer *n, integer *nrhs, doublereal *A, 
			integer *lda, integer *ipv, integer *lidb, doublereal *work, integer *lwork, integer *info);

//******************** Lapack Utility Routines ************************

    IFPACK_DEPRECATED void F77NAME(dgels)(char *trans, integer *m, integer *n, integer *nrhs, doublereal *A,
        integer *lda, doublereal *B, integer *ldb, doublereal *work, integer *lwork, integer *info);

    IFPACK_DEPRECATED void F77NAME(dtimmg)(integer *iflag, integer *m, integer *n, doublereal *A, integer *lda,
                integer *kl, integer *ku);

    IFPACK_DEPRECATED void F77NAME(dlaswp)(integer *n, doublereal *A, integer *lda, integer *k1, integer *k2,
                integer *ipiv, integer *incx);

    doublereal F77NAME(dopla)(char *subname, integer *m, integer *n, integer *kl, integer *ku,
            integer *nb);

// ******************* LU Factorization Routines **********************

    IFPACK_DEPRECATED void F77NAME(dgetrf)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                integer *info);

    IFPACK_DEPRECATED void F77NAME(dgetf2)(integer *m, integer *n, doublereal *A, integer *lda, integer *ipiv,
                integer *info);

    IFPACK_DEPRECATED void F77NAME(dgbtrf)(integer *m, integer *n, integer *KL, integer *KU, doublereal *BM,
                integer *LDBM, integer *ipiv, integer *info);

    IFPACK_DEPRECATED void F77NAME(dgttrf)(integer *N, doublereal *DL, doublereal *D, doublereal *DU,
                doublereal *DU2, integer *ipiv, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpotrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                integer *info);

    IFPACK_DEPRECATED void F77NAME(dsytrf)(char *UPLO, integer *N, doublereal *SM, integer *LDSM,
                integer *ipiv, doublereal *WORK, integer *LWORK, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpbtrf)(char *UPLO, integer *N, integer *KD, doublereal *SBM,
                integer *LDSM, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpttrf)(integer *N, doublereal *D, doublereal *E, integer *info);

// ********************* LU Solve Routines ***************************

    IFPACK_DEPRECATED void F77NAME(dgetrs)(char *trans, integer *N, integer *nrhs, doublereal *A, integer *lda, 
            integer * ipiv, doublereal *b, integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dgbtrs)(char *trans, integer *N, integer *kl, integer *ku, integer *nrhs,
            doublereal *AB, integer *ldab, integer *ipiv, doublereal *b, integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dsytrs)(char *uplo, integer *N, integer *nrhs, doublereal *A, integer *lda, 
            integer *ipiv, doublereal *b, integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dgttrs)(char *trans, integer *N, integer *nrhs, doublereal *DL, 
                doublereal *D, doublereal *DU, doublereal *DU2, integer *ipiv, doublereal *b, 
                integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpotrs)(char *UPLO, integer *N, integer *nrhs, doublereal *A, integer *LDA,
                doublereal *b, integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpttrs)(integer *N, integer *nrhs, doublereal *D, doublereal *E, 
                doublereal *b, integer *ldb, integer *info);

    IFPACK_DEPRECATED void F77NAME(dpbtrs)(char *UPLO, integer *N, integer *KD, integer *nrhs, doublereal *AB,
                integer *LDAB, doublereal *b, integer *ldb, integer *info);

// ********************* Eigen Solve Routines ***************************

	IFPACK_DEPRECATED void F77NAME(dsyev)(char *jobz, char *uplo, integer *N, doublereal *S,
	integer *lda, doublereal *eig, doublereal *work, integer *lwork, integer *info);

// *******************************
}

#endif 
