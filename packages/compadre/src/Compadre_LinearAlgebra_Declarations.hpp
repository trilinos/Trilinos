// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#ifndef _COMPADRE_LINEAR_ALGEBRA_DECLARATIONS_HPP_
#define _COMPADRE_LINEAR_ALGEBRA_DECLARATIONS_HPP_

#include "Compadre_Config.h"
#include "Compadre_Typedefs.hpp"
#include "Compadre_ParallelManager.hpp"

namespace Compadre {

namespace GMLS_LinearAlgebra {

    /*! \brief Calculates two eigenvectors corresponding to two dominant eigenvalues
        \param teamMember    [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param V            [out] - dimension * dimension Kokkos View
        \param PtP           [in] - dimension * dimension Kokkos View
        \param dimensions    [in] - dimension of PtP
    */
    KOKKOS_INLINE_FUNCTION
    void largestTwoEigenvectorsThreeByThreeSymmetric(const member_type& teamMember, scratch_matrix_right_type V, scratch_matrix_right_type PtP, const int dimensions, pool_type& random_number_pool);

    /*! \brief Solves a batch of problems with QR+Pivoting
 
         ~ Note: Very strong assumption on B. ~

         A contains num_matrices * lda * nda data which is num_matrices different (lda x nda) matrices with valid entries of size (M x N), and
         B contains num_matrices * ldb * ndb data which is num_matrices different (ldb x ndb) right hand sides.
         B is assumed to have one of two forms:
         - \f$M==N\f$, the valid entries are the first (N X NRHS)
         - \f$M>N\f$, the valid entries are the first (NRHS)
           i.e. for this case, B is intended to store non-zero entries from a diagonal matrix (as a vector). For the \f$k^{th}\f$ matrix, the (m,m) entry of a diagonal matrix would here be stored in the \f$m^{th}\f$ position.

 

        \param pm                   [in] - manager class for team and thread parallelism
        \param A                [in/out] - matrix A (in), meaningless workspace output (out)
        \param lda                  [in] - row dimension of each matrix in A
        \param nda                  [in] - columns dimension of each matrix in A
        \param B                [in/out] - right hand sides (in), solution (out)
        \param ldb                  [in] - row dimension of each matrix in B
        \param ndb                  [in] - column dimension of each matrix in B
        \param M                    [in] - number of rows containing data (maximum rows over everything in batch) in A
        \param N                    [in] - number of columns containing data in each matrix in A
        \param NRHS                 [in] - number of columns containing data in each matrix in B
        \param num_matrices         [in] - number of problems
        \param implicit_RHS         [in] - determines whether RHS will be stored implicitly. If true, instead of RHS storing the full sqrt(W) explicitly, only the diagonal entries of sqrt(W) will be stored as a 1D array beginning at entry with matrix coordinate (0,0).
    */
    template <typename A_layout=layout_right, typename B_layout=layout_right, typename X_layout=layout_right>
    void batchQRPivotingSolve(ParallelManager pm, double *A, int lda, int nda, double *B, int ldb, int ndb, int M, int N, int NRHS, const int num_matrices, const bool implicit_RHS = true);

} // GMLS_LinearAlgebra
} // Compadre

#endif
