
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

#include "Ifpack_config.h"

#include "ifp_arch.h"

extern "C"
{
    IFPACK_DEPRECATED void F77NAME(dgemm)(char *transa, char *transb, integer *m, integer *n, integer *k,
        double *alpha, const double *a, integer *lda, const double *b, 
        integer *ldb, double *beta, double *c, integer *ldc);
    
    IFPACK_DEPRECATED void F77NAME(dtrsm)(char *side, char *uplo, char *transa, char *diag,
        integer *m, integer *n, double *alpha, const double *A, integer *lda,
        const double *B, integer *ldb);

    IFPACK_DEPRECATED void F77NAME(dtrmm)(char *side, char *uplo, char *transa, char *diag,
        integer *m, integer *n, double *alpha, const double *A, integer *lda,
        const double *B, integer *ldb);

    IFPACK_DEPRECATED void F77NAME(dsymm)(char *side, char *uplo, integer *m, integer *n, 
        double *alpha, const double *A, integer *lda, const double *B, 
        integer *ldb, double *beta, double *C, integer *ldc);

    IFPACK_DEPRECATED void F77NAME(dsyrk)(char *uplo, char *transa, integer *n, integer *k, 
        double *alpha, double *A, integer *lda, double *beta, double *C, 
        integer *ldc);

    IFPACK_DEPRECATED void F77NAME(dsyr2k)(char *uplo, char *transa, integer *n, integer *k, 
        double *alpha, double *A, integer *lda, double *B, integer *ldb,
        double *beta, double *C, integer *ldc);
}

