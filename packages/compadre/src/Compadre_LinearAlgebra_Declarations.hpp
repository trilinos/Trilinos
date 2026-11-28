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
#include "Compadre_Operators.hpp"

namespace Compadre {

namespace GMLS_LinearAlgebra {

    /*! \brief Calculates two eigenvectors corresponding to two dominant eigenvalues
        \param teamMember    [in] - Kokkos::TeamPolicy member type (created by parallel_for)
        \param V            [out] - dimension * dimension Kokkos View
        \param PtP           [in] - dimension * dimension Kokkos View
        \param dimensions    [in] - dimension of PtP
    */
    KOKKOS_INLINE_FUNCTION
    void largestTwoEigenvectorsThreeByThreeSymmetric(const member_type& teamMember, device_unmanaged_matrix_right_type V, scratch_matrix_right_type PtP, const int dimensions, pool_type& random_number_pool);

    /*! \brief Solves a batch of linear equations
 
         ~ Note: Very strong assumption on B. ~

         A contains num_matrices * lda * nda data which is num_matrices different (lda x nda) matrices with valid entries of size (M x N), and
         B contains num_matrices * ldb * ndb data which is num_matrices different (ldb x ndb) right hand sides.
         B is assumed to have one of two forms:
         - \f$M==N\f$, the valid entries are the first (N X NRHS)
         - \f$M>N\f$, the valid entries are the first (NRHS)
           i.e. for this case, B is intended to store non-zero entries from a diagonal matrix (as a vector). For the \f$k^{th}\f$ matrix, the (m,m) entry of a diagonal matrix would here be stored in the \f$m^{th}\f$ position.

        With QR with pivoting, we can solve AX=B, where A=sqrt(W)*P and B=sqrt(W)*Identity or we can solve AX=B where A=P^T*W*P and B=P^T*W. 
        Cholesky, on the other hand, is only set up to solve AX=B where A=P^T*W*P and B=P^T*W.
        QR can therefore either take a full RHS matrix B or the beginning of B is a vector of sqrt(W), in which case the RHS is implicit (implicit_RHS).

        \param pm                   [in] - manager class for team and thread parallelism
        \param dst                  [in] - Compadre::DenseSolverType
        \param A                [in/out] - matrix A (rank 3 A that is num_matrices * lda * nda, in), meaningless workspace output (out)
        \param lda                  [in] - row dimension of each matrix in A
        \param nda                  [in] - columns dimension of each matrix in A
        \param transpose_A          [in] - should A be transposed before use?
        \param B                [in/out] - right hand sides (rank 3 B that is num_matrices * ldb * ndb, in), solution (out)
        \param ldb                  [in] - row dimension of each matrix in B
        \param ndb                  [in] - column dimension of each matrix in B
        \param transpose_B          [in] - should B be transposed before use?
        \param M                    [in] - number of rows containing data (maximum rows over everything in batch) in A (after transpose, if specified)
        \param N                    [in] - number of columns containing data in each matrix in A (after transpose, if specified)
        \param NRHS                 [in] - number of columns containing data in each matrix in B (after transpose, if specified)
        \param num_matrices         [in] - number of problems
        \param implicit_RHS         [in] - determines whether RHS will be stored implicitly. If true, instead of RHS storing the full sqrt(W) explicitly, only the diagonal entries of sqrt(W) will be stored as a 1D array beginning at entry with matrix coordinate (0,0).
        \param rank_full            [in] - whether this A*X=B has an A of full rank
    */
    void batchSolve(ParallelManager pm, DenseSolverType dst, double *A, int lda, int nda, bool transpose_A, double *B, int ldb, int ndb, bool transpose_B, int M, int N, int NRHS, const int num_matrices, const bool implicit_RHS = true, const bool rank_full = true);

} // GMLS_LinearAlgebra
} // Compadre

#endif
