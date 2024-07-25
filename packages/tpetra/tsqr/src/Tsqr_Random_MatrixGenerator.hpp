// @HEADER
// *****************************************************************************
//          Kokkos: Node API and Parallel Node Kernels
//
// Copyright 2008 NTESS and the Kokkos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __TSQR_Random_MatrixGenerator_hpp
#define __TSQR_Random_MatrixGenerator_hpp

#include "Tsqr_Matrix.hpp"
#include "Tsqr_Impl_Lapack.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <algorithm>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace TSQR {
  namespace Random {

    template< class Ordinal, class Scalar, class Generator >
    class MatrixGenerator {
    private:
      typedef Teuchos::ScalarTraits<Scalar> STS;

    public:
      typedef Ordinal ordinal_type;
      typedef Scalar scalar_type;
      typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
      typedef Generator generator_type;

      MatrixGenerator (Generator& generator) : gen_ (generator) {}

      void
      fill_random (const Ordinal nrows,
                   const Ordinal ncols,
                   Scalar A[],
                   const Ordinal lda)
      {
        for (Ordinal j = 0; j < ncols; ++j) {
          Scalar* const A_j = &A[j*lda];
          for (Ordinal i = 0; i < nrows; ++i) {
            A_j[i] = gen_ ();
          }
        }
      }

      /// Fill the nrows by ncols matrix Q (in column-major order, with
      /// leading dimension ldq >= nrows) with a random orthogonal
      /// matrix.
      ///
      /// \note If you want the resulting Q factor to be from a Haar
      /// distribution (the usual one for random orthogonal matrices),
      /// Generator must have a normal distribution.
      void
      explicit_Q (const Ordinal nrows,
                  const Ordinal ncols,
                  Scalar Q[],
                  const Ordinal ldq)
      {
        // Fill Q with random numbers
        this->fill_random (nrows, ncols, Q, ldq);

        // Get ready for QR factorization
        Impl::Lapack<Scalar> lapack;
        std::vector<Scalar> tau (std::min(nrows, ncols));

        // Workspace query
        const int lwork1 =
          lapack.compute_QR_lwork (nrows, ncols, Q, ldq);
        const int lwork2 =
          lapack.compute_explicit_Q_lwork (nrows, ncols, ncols,
                                           Q, ldq, tau.data ());
        const Ordinal lwork = std::max (lwork1, lwork2);
        std::vector<Scalar> work (lwork);

        lapack.compute_QR (nrows, ncols, Q, ldq, tau.data(),
                           work.data(), lwork);

        // Compute explicit Q factor in place
        lapack.compute_explicit_Q (nrows, ncols, ncols,
                                   Q, ldq, tau.data(),
                                   work.data(), lwork);
      }

      /// Fill the nrows by ncols matrix Q (in column-major order, with
      /// leading dimension ldq >= nrows) with a random orthogonal
      /// matrix, stored implicitly.  tau (of length min(nrows,ncols))
      /// is part of this storage.
      ///
      /// \note If you want the resulting Q factor to be from a Haar
      /// distribution (the usual one for random orthogonal matrices),
      /// Generator must have a normal distribution.
      void
      implicit_Q (const Ordinal nrows,
                  const Ordinal ncols,
                  Scalar Q[],
                  const Ordinal ldq,
                  Scalar tau[])
      {
        // Fill Q with random numbers
        this->fill_random (nrows, ncols, Q, ldq);

        Impl::Lapack<Scalar> lapack;
        const int lwork =
          lapack.compute_QR_lwork (nrows, ncols, Q, ldq);
        std::vector<Scalar> work (lwork);
        lapack.compute_QR (nrows, ncols, Q, ldq, tau,
                           work.data (), lwork);
      }

      template< class MatrixViewType >
      void
      implicit_Q (MatrixViewType& Q,
                  typename MatrixViewType::non_const_value_type tau[])
      {
        implicit_Q (Q.extent(0), Q.extent(1), Q.data(), Q.stride(1), tau);
      }

      void
      fill_random_svd (const Ordinal nrows,
                       const Ordinal ncols,
                       Scalar A[],
                       const Ordinal lda,
                       const magnitude_type singular_values[])
      {
        typedef Matrix<Ordinal, Scalar> matrix_type;
        typedef MatView<Ordinal, Scalar> mat_view_type;

        matrix_type U (nrows, ncols, Scalar(0));
        matrix_type V (ncols, ncols, Scalar(0));
        std::vector<Scalar> tau_U (std::min (nrows, ncols));
        std::vector<Scalar> tau_V (ncols);

        // Fill A with zeros, and then make its diagonal the given set
        // of singular values.
        mat_view_type A_view (nrows, ncols, A, lda);
        deep_copy (A_view, Scalar {});
        for (Ordinal j = 0; j < ncols; ++j) {
          A_view(j,j) = Scalar (singular_values[j]);
        }

        // Generate random orthogonal U (nrows by ncols) and V (ncols by
        // ncols).  Keep them stored implicitly.
        implicit_Q (U, tau_U.data());
        implicit_Q (V, tau_V.data());

        // Workspace query for ORMQR.
        Impl::Lapack<Scalar> lapack;
        const int lwork1 =
          lapack.apply_Q_factor_lwork ('L', 'N', nrows, ncols, ncols,
                                       U.data (), U.stride (1),
                                       tau_U.data (), A, lda);
        int lwork2 = 0;
        if (STS::isComplex) {
          lwork2 =
            lapack.apply_Q_factor_lwork ('R', 'C',
                                         nrows, ncols, ncols,
                                         V.data (), V.stride (1),
                                         tau_V.data (), A, lda);
        }
        else {
          lwork2 =
            lapack.apply_Q_factor_lwork ('R', 'T',
                                         nrows, ncols, ncols,
                                         V.data (), V.stride (1),
                                         tau_V.data (), A, lda);
        }

        // Allocate workspace.
        Ordinal lwork (std::max (lwork1, lwork2));
        std::vector<Scalar> work (lwork);

        // Apply U to the left side of A, and V^H to the right side of A.
        lapack.apply_Q_factor ('L', 'N', nrows, ncols, ncols,
                               U.data(), U.stride(1), tau_U.data(),
                               A, lda, work.data(), lwork);
        if (STS::isComplex) {
          lapack.apply_Q_factor ('R', 'C', nrows, ncols, ncols,
                                 V.data(), V.stride(1), tau_V.data(),
                                 A, lda, work.data(), lwork);
        }
        else {
          lapack.apply_Q_factor ('R', 'T', nrows, ncols, ncols,
                                 V.data(), V.stride(1), tau_V.data(),
                                 A, lda, work.data(), lwork);
        }
      }

      /// \brief Fill in a random upper triangular matrix
      ///
      /// Fill R with a random n by n upper triangular matrix, with
      /// the specified singular values.
      ///
      /// \param n [in] Dimension of R (n by n)
      ///
      /// \param R [out] On output: n by n upper triangular matrix
      ///   with stride ldr, stored in non-packed, column-major
      ///   format.
      ///
      /// \param ldr [in] Stride (a.k.a. leading dimension) of R
      ///
      /// \param singular_values [in] n nonnegative real values, to be
      ///   used as the singular values of the output R.
      void
      fill_random_R (const Ordinal n,
                     Scalar R[],
                     const Ordinal ldr,
                     const magnitude_type singular_values[])
      {
        // Fill R with an n x n (not upper triangular) random matrix
        // having the given singular values.
        fill_random_svd (n, n, R, ldr, singular_values);

        // Compute the QR factorization in place of R (which isn't upper triangular yet).
        std::vector<Scalar> tau (n);

        // Workspace size query for QR factorization.
        Impl::Lapack<Scalar> lapack;
        const int lwork = lapack.compute_QR_lwork (n, n, R, ldr);

        // Compute QR factorization (implicit representation in place).
        std::vector<Scalar> work (lwork);
        lapack.compute_QR (n, n, R, ldr, tau.data (),
                           work.data (), lwork);

        // Zero out the stuff below the diagonal of R, leaving just the R factor.
        for (Ordinal j = 0; j < n; ++j) {
          for (Ordinal i = j+1; i < n; ++i) {
            R[i + j*ldr] = Scalar {};
          }
        }
      }

    private:
      static Ordinal
      checkedCast (const magnitude_type& x)
      {
        if (x < std::numeric_limits<Ordinal>::min() || x > std::numeric_limits<Ordinal>::max()) {
          throw std::range_error("Scalar input cannot be safely cast to an Ordinal");
        }
        else if (std::numeric_limits< magnitude_type >::is_signed &&
                 x < magnitude_type(0) &&
                 ! std::numeric_limits<Ordinal>::is_signed) {
          throw std::range_error("Scalar input is negative, but Ordinal is unsigned");
        }
        else {
          return static_cast<Ordinal> (x);
        }
      }

      Generator& gen_;
    };
  } // namespace Random
} // namespace TSQR

#endif // __TSQR_Random_MatrixGenerator_hpp
