// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Tsqr_Random_GlobalMatrix_hpp
#define __Tsqr_Random_GlobalMatrix_hpp

#include "Tsqr_Matrix.hpp"
#include "Tsqr_Random_MatrixGenerator.hpp"
#include "Tsqr_RMessenger.hpp"
#include "Tsqr_Impl_SystemBlas.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace TSQR {
  namespace Random {

    template<class MatrixViewType>
    static void
    scaleMatrix (MatrixViewType& A,
                 const typename MatrixViewType::non_const_value_type& denom)
    {
      using LO = typename MatrixViewType::ordinal_type;

      const LO nrows = A.extent(0);
      const LO ncols = A.extent(1);
      const LO lda = A.stride(1);
      if (nrows == lda) { // A is stored contiguously.
        const LO nelts = nrows * ncols;
        auto A_ptr = A.data ();
        for (LO k = 0; k < nelts; ++k) {
          A_ptr[k] /= denom;
        }
      }
      else { // Each column of A is stored contiguously.
        for (LO j = 0; j < ncols; ++j) {
          auto A_j = &A(0,j);
          for (LO i = 0; i < nrows; ++i) {
            A_j[i] /= denom;
          }
        }
      }
    }

    template<class MatrixViewType, class Generator>
    void
    randomGlobalMatrix (Generator* const pGenerator,
                        MatrixViewType& A_local,
                        const typename Teuchos::ScalarTraits<typename MatrixViewType::non_const_value_type>::magnitudeType singular_values[],
                        MessengerBase< typename MatrixViewType::ordinal_type>* const ordinalMessenger,
                        MessengerBase< typename MatrixViewType::non_const_value_type>* const scalarMessenger)
    {
      using Teuchos::NO_TRANS;
      using ordinal_type = typename MatrixViewType::ordinal_type;
      using scalar_type = typename MatrixViewType::non_const_value_type;
      using STS = Teuchos::ScalarTraits<scalar_type>;

      const int rootProc = 0;
      const int nprocs = ordinalMessenger->size();
      const int myRank = ordinalMessenger->rank();
      Impl::SystemBlas<scalar_type> blas;

      const ordinal_type nrowsLocal = A_local.extent(0);
      const ordinal_type ncols = A_local.extent(1);

      // Theory: Suppose there are P processors.  Proc q wants an m_q by n
      // component of the matrix A, which we write as A_q.  On Proc 0, we
      // generate random m_q by n orthogonal matrices Q_q (in explicit
      // form), and send Q_q to Proc q.  The m by n matrix [Q_0; Q_1; ...;
      // Q_{P-1}] is not itself orthogonal.  However, the m by n matrix
      // Q = [Q_0 / P; Q_1 / P; ...; Q_{P-1} / P] is orthogonal:
      //
      // \sum_{q = 0}^{P-1} (Q_q^T * Q_q) / P = I.

      if (myRank == rootProc) {
        using matgen_type = Random::MatrixGenerator<ordinal_type,
          scalar_type, Generator>;
        matgen_type matGen (*pGenerator);

        // Generate a random ncols by ncols upper triangular matrix R
        // with the given singular values.
        Matrix<ordinal_type, scalar_type> R (ncols, ncols, scalar_type {});
        matGen.fill_random_R (ncols, R.data(), R.stride(1), singular_values);

        // Broadcast R to all the processors.
        scalarMessenger->broadcast (R.data(), ncols*ncols, rootProc);

        // Generate (for myself) a random nrowsLocal x ncols
        // orthogonal matrix, stored in explicit form.
        Matrix<ordinal_type, scalar_type> Q_local (nrowsLocal, ncols);
        matGen.explicit_Q (nrowsLocal, ncols, Q_local.data(), Q_local.stride(1));

        // Scale the (local) orthogonal matrix by the number of
        // processors P, to make the columns of the global matrix Q
        // orthogonal.  (Otherwise the norm of each column will be P
        // instead of 1.)
        const scalar_type P (static_cast<double> (nprocs));
        // Do overflow check.  If casting P back to scalar_type
        // doesn't produce the same value as nprocs, the cast
        // overflowed.  We take the real part, because scalar_type
        // might be complex.
        if (nprocs != static_cast<int> (STS::real (P))) {
          throw std::runtime_error ("Casting nprocs to Scalar failed");
        }

        scaleMatrix (Q_local, P);

        // A_local := Q_local * R
        blas.GEMM (NO_TRANS, NO_TRANS, nrowsLocal, ncols, ncols,
                   scalar_type(1), Q_local.data(), Q_local.stride(1),
                   R.data(), R.stride(1),
                   scalar_type(0), A_local.data(), A_local.stride(1));

        for (int recvProc = 1; recvProc < nprocs; ++recvProc) {
          // Ask the receiving processor how big (i.e., how many rows)
          // its local component of the matrix is.
          ordinal_type nrowsRemote = 0;
          ordinalMessenger->recv (&nrowsRemote, 1, recvProc, 0);

          // Make sure Q_local is big enough to hold the data for
          // the current receiver proc.
          Q_local.reshape (nrowsRemote, ncols);

          // Compute a random nrowsRemote * ncols orthogonal
          // matrix Q_local, for the current receiving processor.
          matGen.explicit_Q (nrowsRemote, ncols, Q_local.data(), Q_local.stride(1));

          // Send Q_local to the current receiving processor.
          scalarMessenger->send (Q_local.data(), nrowsRemote*ncols, recvProc, 0);
        }
      }
      else {
        // Receive the R factor from Proc 0.  There's only 1 R
        // factor for all the processes.
        Matrix<ordinal_type, scalar_type> R (ncols, ncols, scalar_type {});
        scalarMessenger->broadcast (R.data(), ncols*ncols, rootProc);

        // Q_local (nrows_local by ncols, random orthogonal matrix)
        // will be received from Proc 0, where it was generated.
        const ordinal_type recvSize = nrowsLocal * ncols;
        Matrix<ordinal_type, scalar_type> Q_local (nrowsLocal, ncols);

        // Tell Proc 0 how many rows there are in the random orthogonal
        // matrix I want to receive from Proc 0.
        ordinalMessenger->send (&nrowsLocal, 1, rootProc, 0);

        // Receive the orthogonal matrix from Proc 0.
        scalarMessenger->recv (Q_local.data(), recvSize, rootProc, 0);

        // Scale the (local) orthogonal matrix by the number of
        // processors, to make the global matrix Q orthogonal.
        const scalar_type P (static_cast<double> (nprocs));
        // Do overflow check.  If casting P back to scalar_type
        // doesn't produce the same value as nprocs, the cast
        // overflowed.  We take the real part, because scalar_type
        // might be complex.
        if (nprocs != static_cast<int> (STS::real (P))) {
          throw std::runtime_error ("Casting nprocs to Scalar failed");
        }
        scaleMatrix (Q_local, P);

        // A_local := Q_local * R
        blas.GEMM (NO_TRANS, NO_TRANS, nrowsLocal, ncols, ncols,
                   scalar_type(1), Q_local.data(), Q_local.stride(1),
                   R.data(), R.stride(1),
                   scalar_type(0), A_local.data(), A_local.stride(1));
      }
    }
  } // namespace Random
} // namespace TSQR

#endif // __Tsqr_Random_GlobalMatrix_hpp
